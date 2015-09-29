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
    end
    return merged_img, [global_ref.i, global_ref.j]
end

"""
`IMFUSE` - Overlay two images on top of each other using their offsets. Colors 
one image red, the other green, and the overlap yellow.
Uses rounded interpolation.

Args:

* A: image A (2D array)
* offset_A: 2-element array for global position of image A
* B: image B (2D array)
* offset_B: 2-element array for global position of image A

Returns:

* O: Image object combining both image A & B

    `imfuse(A, offset_A, B, offset_B)`
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
    #O = Overlay((A,B), (RGB(1.,0.,0.), RGB(0.,0.6,0.)))
    BB_O = min(offset_A, offset_B)
    return O, BB_O
end

"""
`PADIMAGE` - Specify image padding in each of four directions
    
    `new_img = padimage(img, xlow, ylow, xhigh, yhigh)`

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
  z = zeros(Ufixed8, bb.h+1, bb.w+1)
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

"""
`MATCH_PADDING` - For each dimension, pad images to match the larger of the two

  `imgA, imgB = match_padding(imgA, imgB)`
"""
function match_padding(imgA, imgB)
  szA = collect(size(imgA))
  szB = collect(size(imgB))
  szC = max(szA, szB)
  imgA = padimage(imgA, 0, 0, reverse(szC-szA)...)
  imgB = padimage(imgB, 0, 0, reverse(szC-szB)...)
  return imgA, imgB
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

"""
Cycle through JLD files in montaged directory and render montage
"""
function render_montaged(section_range::Array{Int64})
  log_path = joinpath(MONTAGED_DIR, "montaged_offsets.txt")
  filenames = sort_dir(MONTAGED_DIR)[section_range]
  for filename in filenames
    println("Rendering ", filename[1:end-4])
    meshset = JLD.load(joinpath(MONTAGED_DIR, filename))["MeshSet"]
    warps = pmap(meshwarp, meshset.meshes)
    img, offset = merge_images([[x[i] for x in warps] for i=1:2]...)
    img = grayim(img)
    img["spatialorder"] = ["y","x"]
    println("Writing ", filename[1:end-4])
    new_filename = string(filename[1:end-4], ".tif")
    @time imwrite(img, joinpath(MONTAGED_DIR, new_filename))
    update_offset_log!(log_path, new_filename, [0,0], size(img))

    imfuse_section(meshset)
  end
end

"""
Write thumbnail image with vectors and match indices overlayed
"""
function write_thumbnail(path, img, vectors, factor)
  imgc, img2 = view(img, pixelspacing=[1,1])
  a, b = draw_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1), factor)
  c = draw_indices(imgc, img2, vectors[1:2,:])
  println("Writing ", path)
  Cairo.write_to_png(imgc.c.back, path)
  destroy(toplevel(imgc))
end

"""
Calculate prealignment transforms from first section through section_num
"""
function calculate_global_tform(index, dir=PREALIGNED_DIR)
  global_tform = eye(3)
  if index != (1,1,-2,-2)
    index_pairs = create_sequential_index_pairs((1,1,-2,-2), index)
    for (indexA, indexB) in index_pairs
      meshset = load(indexA, indexB)
      # tform = affine_approximate(meshset)
      tform = regularized_approximate(meshset, lambda=0.9)
      global_tform *= tform
    end
  end
  return global_tform
end

"""
INCOMPLETE

Lazy function to generate list of indices between indexA and indexB

Fixes to include:
* switch between wafers
* check indexA[3:4] against indexB[3:4]
"""
function create_range(indexA, indexB)
  return [(indexA[1], i, indexA[3:4]...) for i in indexA[2]:indexB[2]]
end

function create_sequential_index_pairs(indexA, indexB)
  indices = create_range(indexA, indexB)
  return zip(indices[1:end-1], indices[2:end])
end

"""
Copy a section from one process step to the next
"""
function copy_section_through(index)
  println("copy_section_through INCOMPLETE")
  return
end

"""
Test if an index is the first section in the stack
"""
function is_first_section(index)
  return index[1] == 1 && index[2] == 1
end

"""
Find first row of the offset file that matches index and return cumulative sum 
of previous offset arrays
"""
function find_offset(offset_file, index)
  if findfirst(offset_file[:,2], index) != 0
    return collect(sum(offset_file[1:findfirst(offset_file[:,2], index), 3:4], 1))
  else
    return [0,0]
  end
end

"""
Find appropriate offset file and pull out the offset array for the index
"""
function load_offset(index)
  if is_montaged(index)
    return find_offset(MONTAGED_OFFSETS, index)
  elseif is_aligned(index)
    return find_offset(PREALIGNED_OFFSETS, index)
  else
    return [0,0]
  end
end

"""
Return Dictionary of staged image to remove redundancy in loading
"""
function stage_image(mesh, tform, scale=0.05)
  s = [scale 0 0; 0 scale 0; 0 0 1]
  stage = Dict()
  stage["index"] = mesh.index
  img = get_ufixed8_image(mesh)
  montaged_offset = load_offset(mesh.index)
  println("tform: ", tform)
  println("montaged_offset: ", montaged_offset)
  println("Warping ", mesh.name)
  stage["img"], stage["offset"] = imwarp(img, tform, montaged_offset)
  println("Warping thumbnail for ", mesh.name)
  stage["thumb"], stage["thumb_offset"] = imwarp(img, tform*s, montaged_offset)
  stage["scale"] = scale
  return stage
end

"""
Prealignment where offsets are global
"""
function render_prealigned(indexA, indexB)
  dir = PREALIGNED_DIR
  fixed = Dict()

  global_tform = calculate_global_tform(indexA)
  log_path = joinpath(dir, "prealigned_offsets.txt")

  function save_image(stage, dir, log_path)
    fn = string(join(stage["index"][1:2], ","), "_prealigned.tif")
    update_offset_log!(log_path, fn, stage["offset"], size(stage["img"]))
    println("Writing ", fn)
    @time imwrite(stage["img"], joinpath(dir, fn))
  end

  function save_thumbnails(A, B)
    fn = string(join(B["index"][1:2], ","), "_prealigned_thumbnail.png")
    println("Saving thumbnail ", fn)
    path = joinpath(dir, "review", fn)
    O, O_bb = imfuse(A["thumb"], A["thumb_offset"], B["thumb"], B["thumb_offset"])
    moving_nodes = A["nodes"][:,1:2]'*A["scale"]
    fixed_nodes = B["nodes"][:,1:2]'*B["scale"]
    moving_nodes .-= O_bb
    fixed_nodes .-= O_bb
    vectors = [moving_nodes; fixed_nodes]
    write_thumbnail(path, O, vectors, 1.0)
  end

  index_pairs = create_sequential_index_pairs(indexA, indexB)
  for (k, (indexA, indexB)) in enumerate(index_pairs)
    println("\nPrealigning ", indexA, " & ", indexB)
    meshset = load(indexA, indexB)
    if k==1
      fixed = stage_image(meshset.meshes[1], global_tform)
      if is_first_section(indexA)
        save_image(fixed, dir, log_path)
      end
    end
    tform = regularized_approximate(meshset, lambda=0.9)
    moving = stage_image(meshset.meshes[2], global_tform*tform)
    save_image(moving, dir, log_path)
    moving_nodes, fixed_nodes = get_matched_points(meshset, 1)
    fixed["nodes"] = points_to_Nx3_matrix(fixed_nodes)*global_tform
    moving["nodes"] = points_to_Nx3_matrix(moving_nodes)*global_tform*tform
    save_thumbnails(fixed, moving)
    fixed = 0
    fixed = moving
    moving = 0
    global_tform *= tform
  end
end

"""
Cycle through JLD files in aligned directory and render alignment
"""
function render_aligned(file_index)
  dir = ALIGNED_DIR
  scale = 0.0625
  s = [scale 0 0; 0 scale 0; 0 0 1]

  # Log file for image offsets
  log_path = joinpath(dir, "aligned_offsets.txt")

  filename = sort_dir(dir)[file_index]
  println("Rendering meshes in ", filename)
  meshset = JLD.load(joinpath(dir, filename))["MeshSet"]
  images = Dict()
  
  # Check images dict for thumbnail, otherwise render, save, & resize it
  function retrieve_image(mesh)
    index = (mesh.index[1:2]..., -4, -4)
    if index in keys(images)
      img = images[index]
    else
      path = get_path(index)
      if isfile(path)
        img, _ = imwarp(get_ufixed8_image(index), s)
      else
        println("Warping ", mesh.name)
        @time img, offset = meshwarp(mesh)
        @time img = rescopeimage(img, offset, GLOBAL_BB)
        println("Writing ", mesh.name)
        new_fn = string(join(mesh.index[1:2], ","), "_aligned.tif")
        @time imwrite(img, joinpath(dir, new_fn))
        img, _ = imwarp(img, s)

        # Log image offsets
        update_offset_log!(log_path, new_fn, offset, size(img))
      end
      images[index] = img
    end
    return img
  end

  # map(warp_pad_write, meshset.meshes)
  for (k, matches) in enumerate(meshset.matches)
    src_index = matches.src_index
    dst_index = matches.dst_index
    src_mesh = meshset.meshes[find_index(meshset, src_index)]
    dst_mesh = meshset.meshes[find_index(meshset, dst_index)]

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

    src_nodes = hcat(src_nodes...)[1:2, :]*scale .- src_offset
    dst_nodes = hcat(dst_nodes...)[1:2, :]*scale .- dst_offset

    imgc, img2 = view(O, pixelspacing=[1,1])
    vectors = [src_nodes; dst_nodes]
    an_pts, an_vectors = draw_vectors(imgc, img2, vectors, RGB(0,0,1), RGB(1,0,1))
    c = draw_indices(imgc, img2, src_nodes)
    # an_src_pts = draw_points(imgc, img2, src_nodes, RGB(1,1,1), 2.0, 'x')
    # an_dst_pts = draw_points(imgc, img2, dst_nodes, RGB(0,0,0), 2.0, '+')

    thumbnail_fn = string(join(dst_index[1:2], ","), "-", join(src_index[1:2], ","), "_aligned_thumbnail.png")
    println("Writing ", thumbnail_fn)
    Cairo.write_to_png(imgc.c.back, joinpath(dir, "review", thumbnail_fn))
    destroy(toplevel(imgc))
  end
end

function write_alignment_blockmatches(section_range::Array{Int64})
  filenames = sort_dir(ALIGNED_DIR)[section_range]
  for filename in filenames
    println("Rendering meshes in ", filename)
    meshset = JLD.load(joinpath(ALIGNED_DIR, filename))["MeshSet"]
    save_blockmatch_imgs(meshset, k, [], joinpath(ALIGNED_DIR, "blockmatches"))
  end
end  
