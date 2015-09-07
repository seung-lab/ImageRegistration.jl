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
  i_start = bb.i - offset[1]
  j_start = bb.j - offset[2]
  i_end = (offset[1] + size(img,1)) - (bb.i + bb.h)
  j_end = (offset[2] + size(img,2)) - (bb.j + bb.w)
  i_start = i_start > 0 ? i_start : 1
  j_start = j_start > 0 ? j_start : 1
  i_end = i_end > 0 ? bb.i + bb.h : size(img,1)
  j_end = j_end > 0 ? bb.j + bb.w : size(img,2)
  img = img[i_start:i_end, j_start:j_end]
  i_low = offset[1] - bb.i + 1
  j_low = offset[2] - bb.j + 1
  i_high = offset[1] - bb.i + size(img,1)
  j_high = offset[2] - bb.j + size(img,2)
  z = zeros(bb.h, bb.w)
  z[i_low:i_high, j_low:j_high] = img
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
  for fn in filenames[4:end]
    println("Rendering ", fn[1:end-4])
    meshset = JLD.load(joinpath(MONTAGED_DIR, fn))["MeshSet"]
    img, offset = merge_images(meshset.meshes)
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
  filenames = sort_dir(ALIGNED_DIR)[section_range]
  for filename in filenames
    println("Rendering meshes in ", filename)
    meshset = JLD.load(joinpath(ALIGNED_DIR, filename))["MeshSet"]
    global_bb = BoundingBox(0,0,36000,36000)
    # map(warp_pad_write, meshset.meshes)
    for mesh in meshset.meshes[5:5]
      println("Warping ", mesh.name)
      @time img, offset = meshwarp(mesh)
      println(global_bb)
      println(offset)
      println(size(img))
      println("Writing ", mesh.name)
      new_fn = string(join(mesh.index[1:2], ","), "_aligned.tif")
      @time imwrite(rescopeimage(img, offset, global_bb), 
                                                joinpath(ALIGNED_DIR, new_fn))
    end
  end
end
