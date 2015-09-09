# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150807
#
# Functions to load meshes, warp images via piecewise affine transforms, and
# display images with meshes.

"""
`MERGE_IMAGES` - Place images in global reference image

    merged_img, [bb.i, bb.j] = merge_images(imgs, offsets)

* `imgs`: 1D array, images (2D arrays)
* `offsets`: 1D array, 2-element array positions of corresponding image 
  in 'imgs` in global space

""" 
function merge_images(imgs, offsets)
    bbs = []
    for (img, offset) in zip(imgs, offsets)
        push!(bbs, BoundingBox(offset..., size(img)...))
    end
    global_ref = sum(bbs)
    merged_img = zeros(global_ref.h, global_ref.w)
    for (idx, (img, bb)) in enumerate(zip(imgs, bbs))
        println("Merging tile ", idx)
        i = bb.i - global_ref.i+1
        j = bb.j - global_ref.j+1
        w = bb.w-1
        h = bb.h-1
        merged_img[i:i+h, j:j+w] = max(merged_img[i:i+h, j:j+w], img)
        imgs[idx] = 0
        gc(); gc();
    end
    return merged_img, [global_ref.i, global_ref.j]
end

"""
`MERGE_IMAGES` - Using array of Mesh objects
"""
function merge_images(meshes)
    warps = map(meshwarp, meshes)
    return merge_images([[x[i] for x in warps] for i=1:2]...)
end

"""
`MERGE_IMAGES_PARALLEL` - Using array of Mesh objects, in parallel
"""
function merge_images_parallel(meshes)
    warps = pmap(meshwarp, meshes)
    return merge_images([[x[i] for x in warps] for i=1:2]...)
end   

"""
Overlay two images on top of each other using their offsets. Colors one
image red, the other green, and the overlap yellow.
Uses rounded interpolation.

Args:

* A: image A (2D array)
* offset_A: 2-element array for global position of image A
* B: image B (2D array)
* offset_B: 2-element array for global position of image A

Returns:

* O: Image object combining both image A & B

    imfuse(A, offset_A, B, offset_B)
"""
function imfuse(A, offset_A, B, offset_B)
    # pad to common origin
    BB_C = offset_B - offset_A
    if BB_C[1] > 0
        B = padimage(B, 0, BB_C[1], 0, 0)
    elseif BB_C[1] < 0
        A = padimage(A, 0, -BB_C[1], 0, 0)
    end 
    if BB_C[2] > 0
        B = padimage(B, BB_C[2], 0, 0, 0)
    elseif BB_C[2] < 0
        A = padimage(A, -BB_C[2], 0, 0, 0)
    end 
    # pad to match sizes
    szA = collect(size(A))
    szB = collect(size(B))
    szC = szB - szA
    if szC[1] > 0
        A = padimage(A, 0, 0, 0, szC[1])
    elseif szC[1] < 0
        B = padimage(B, 0, 0, 0, -szC[1])
    end 
    if szC[2] > 0
        A = padimage(A, 0, 0, szC[2], 0)
    elseif szC[2] < 0
        B = padimage(B, 0, 0, -szC[2], 0)
    end
    O = Overlay((A,B), (RGB(1,0,0), RGB(0,1,0)))
    BB_O = min(offset_A, offset_B)
    return O, BB_O
end

"""
`PADIMAGE` - Specify image padding in each of four directions
    
    new_img = padimage(img, xlow, ylow, xhigh, yhigh)

     _________________________________  
    |                                 |  
    |             ylow                |  
    |         ______________          |  
    |        |              |         |  
    |  xlow  |     img      |  xhigh  |  
    |        |              |         |  
    |        |______________|         |  
    |                                 |  
    |             yhigh               |  
    |_________________________________|  

Args:

* img: 2D or 3D array
* xlow: amount to pad in x prior to the image
* ylow: amount to pad in y prior to the image
* xhigh: amount to pad in x after to the image
* yhigh: amount to pad in y after to the image

Returns:

* new_img: original img, extended with rows and columns of zeros
"""
function padimage(img, xlow::Int64, ylow::Int64, xhigh::Int64, yhigh::Int64)
    h = ylow + size(img, 1) + yhigh
    w = xlow + size(img, 2) + xhigh
    z = zeros(h, w)
    z[ylow+1:ylow+size(img,1), xlow+1:xlow+size(img,2)] = img
    return z
end

"""
`RESCOPE` - Crop/pad an image to fill a bounding box
    
    new_img = rescope(img, offset, boundingbox)

Args:

* img: 2D or 3D array
* offset: 2-element array, specifying i & j offset from global origin
* bb: bounding box object in the global reference space

Returns:

* new_img: original img, cropped &/or extended with rows and columns of zeros
"""
function rescopeimage(img, offset, bb)
  z = zeros(bb.h+1, bb.w+1)
  imgbb = BoundingBox(offset..., size(img,1)-1, size(img,2)-1)
  xbb = imgbb - bb
  if !isnan(xbb.i) || !isnan(xbb.j) || !isnan(xbb.h) || !isnan(xbb.h)
    crop_img = xbb.i-offset[1]+1 : xbb.i-offset[1]+1+xbb.h, 
                  xbb.j-offset[2]+1 : xbb.j-offset[2]+1+xbb.w
    crop_z = xbb.i-bb.i+1:xbb.i-bb.i+1+xbb.h, xbb.j-bb.j+1:xbb.j-bb.j+1+xbb.w
    z[crop_z...] = img[crop_img...]
  end
  return z
end

function padimages(imgA, imgB)
  szA = collect(size(imgA))
  szB = collect(size(imgB))
  szC = max(szA, szB)
  imgA = padimage(imgA, 0, 0, reverse(szC-szA)...)
  imgB = padimage(imgB, 0, 0, reverse(szC-szB)...)
  return imgA, imgB
end

function sort_dir(dir, file_extension="jld")
    files_in_dir = filter(x -> x[end-2:end] == file_extension, readdir(dir))
    return sort(files_in_dir, by=x->parseName(x))
end

function get_global_bb(meshset)
    bbs = []
    println("Calculating global bounding box")
    for mesh in meshset.meshes
        nodes = hcat(mesh.nodes_t...)'
        push!(bbs, snap_bb(find_mesh_bb(nodes)))
    end
    global_bb = sum(bbs)
    global_bb.h += 1
    global_bb.w += 1
    println(global_bb)
    return global_bb
end    

function render_montaged(section_range::UnitRange{Int64})
  filenames = sort_dir(MONTAGED_DIR)[section_range]
  for fn in filenames
    println("Rendering ", fn[1:end-4])
    meshset = JLD.load(joinpath(MONTAGED_DIR, fn))["MeshSet"]
    img, offset = merge_images_parallel(meshset.meshes)
    img = grayim(img)
    img["spatialorder"] = ["y", "x"]
    println("Writing ", fn[1:end-4])
    @time imwrite(img, joinpath(MONTAGED_DIR, string(fn[1:end-4], ".tif")))

    imfuse_section(meshset)
  end
end

function warp_pad_write(mesh)
    println("Warping ", mesh.name)
    @time img, offset = meshwarp(mesh)
    println(offset)
    println(size(img))
    img = rescopeimage(img, offset, global_bb)
    println(size(img))
    println("Writing ", mesh.name)
    @time imwrite(img, joinpath(ALIGNED_DIR, string(mesh.name, ".tif")))
end

function render_aligned(section_range::UnitRange{Int64})
  scale = 0.0625
  s = [scale 0 0; 0 scale 0; 0 0 1]

  # Log file for image offsets
  log_path = joinpath(ALIGNED_DIR, "aligned_offsets.txt")
  if !isfile(log_path)
    f = open(log_path, "w")
    close(f)
  end

  filenames = sort_dir(ALIGNED_DIR)[section_range]
  for filename in filenames
    println("Rendering meshes in ", filename)
    meshset = JLD.load(joinpath(ALIGNED_DIR, filename))["MeshSet"]
    images = Dict()
    
    # Check images dict for thumbnail, otherwise render, save, & resize it
    function retrieve_image(mesh)
      index = (mesh.index[1:2]..., -4, -4)
      if index in keys(images)
        img = images[index]
      else
        path = getPath(index)
        if isfile(path)
          img, _ = imwarp(getUfixed8Image(index), s)
        else
          println("Warping ", mesh.name)
          @time img, offset = meshwarp(mesh)
          @time img = rescopeimage(img, offset, GLOBAL_BB)
          println("Writing ", mesh.name)
          new_fn = string(join(mesh.index[1:2], ","), "_aligned.tif")
          @time imwrite(img, joinpath(ALIGNED_DIR, new_fn))
          img, _ = imwarp(img, s)

          # Log image offsets
          log_file = open(log_path, "a")
          log_line = join((new_fn, offset[1], offset[2], 
                              size(img,1), size(img,2)), " ")
          write(log_file, log_line, "\n")
          close(log_file)
        end
        images[index] = img
      end
      return img
    end

    # map(warp_pad_write, meshset.meshes)
    for (k, matches) in enumerate(meshset.matches[1:1])
      src_index = matches.src_index
      dst_index = matches.dst_index
      src_mesh = meshset.meshes[findIndex(meshset, src_index)]
      dst_mesh = meshset.meshes[findIndex(meshset, dst_index)]

      src_nodes, dst_nodes = get_matched_points_t(meshset, k)
      src_index = (src_index[1:2]..., src_index[3]-1, src_index[4]-1)
      dst_index = (dst_index[1:2]..., dst_index[3]-1, dst_index[4]-1)
      src_offset = [GLOBAL_BB.i, GLOBAL_BB.j]
      dst_offset = [GLOBAL_BB.i, GLOBAL_BB.j]

      src_img = retrieve_image(src_mesh)
      dst_img = retrieve_image(dst_mesh)

      src_offset *= scale
      dst_offset *= scale

      O, O_bb = imfuse(src_img, src_offset, dst_img, dst_offset)

      src_nodes = (hcat(src_nodes...)[1:2, :] .- src_offset)*scale
      dst_nodes = (hcat(dst_nodes...)[1:2, :] .- dst_offset)*scale

      imgc, img2 = view(O, pixelspacing=[1,1])
      vectors = [src_nodes; dst_nodes]
      an_pts, an_vectors = draw_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1))
      c = draw_indices(imgc, img2, src_nodes)
      # an_src_pts = draw_points(imgc, img2, src_nodes, RGB(1,1,1), 2.0, 'x')
      # an_dst_pts = draw_points(imgc, img2, dst_nodes, RGB(0,0,0), 2.0, '+')

      thumbnail_fn = string(join(dst_index[1:2], ","), "-", join(src_index[1:2], ","), "_aligned_thumbnail.png")
      println("Writing ", thumbnail_fn)
      Cairo.write_to_png(imgc.c.back, joinpath(ALIGNED_DIR, "review", thumbnail_fn))
      destroy(toplevel(imgc))
    end
  end
end

function write_alignment_blockmatches(section_range::UnitRange{Int64})
  filenames = sort_dir(ALIGNED_DIR)[section_range]
  for filename in filenames
    println("Rendering meshes in ", filename)
    meshset = JLD.load(joinpath(ALIGNED_DIR, filename))["MeshSet"]
    save_blockmatch_imgs(meshset, k, [], joinpath(ALIGNED_DIR, "blockmatches"))
  end
end  