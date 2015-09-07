#=
include("convolve.jl")


using Images
include("imwarp.jl")
include("visualize.jl")
#include("render.jl") # cyclic inclusion
=#

function test()
  getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)))
  moving_section = getimage("../sections/S2-W001_fixed_section4_0.175.tif")
  #fixed_section = getimage("./sections/S2-W001_fixed_section4_0.175.tif")
  fixed_section = getimage("../sections/S2-W001_fixed_section5_0.175.tif")
  println(size(moving_section))
  println(size(fixed_section))
  trans, moving_points, fixed_points, res1, res2 = affine_align_sections(moving_section, fixed_section, PARAMS_PREALIGNMENT; return_points=true)
  moving_points = points_to_3xN_matrix(moving_points)
  fixed_points = points_to_3xN_matrix(fixed_points)

  println(trans)
  #out_img, offset = imwarp(fixed_section, trans)
  #imwrite(moving_section, joinpath(".","test_outputs", string("moving_section", ".tif")))
  #imwrite(fixed_section, joinpath(".","test_outputs", string("fixed_section_", offset[1], "_", offset[2], ".tif")))
  #imwrite(out_img, joinpath(".","test_outputs", string("warped_", offset[1], "_", offset[2], ".tif")))
  p22 = fixed_points + res2
  p11 = moving_points + res1
  moving_points = moving_points[1:2,:]
  fixed_points = fixed_points[1:2,:]
  draw_vectors(fixed_section, vcat(fixed_points[1:2,:], p22[1:2,:]))
  #ccp.write_image_from_points(moving_points[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
  draw_points(moving_section, moving_points)
  draw_points(fixed_section, fixed_points)
  #draw_points(out_img, moving_points)
end


function test2()
  getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)))
  moving_section = getimage("../output_images/(1,1)_montage.tif")
  fixed_section = getimage("../output_images/(1,2)_montage.tif")
  println(size(moving_section))
  println(size(fixed_section))
  trans, moving_points, fixed_points, res1, res2 = affine_align_sections(moving_section, fixed_section, PARAMS_PREALIGNMENT; return_points=true)
  moving_points = points_to_3xN_matrix(moving_points)
  fixed_points = points_to_3xN_matrix(fixed_points)

  println(trans)
  downsample = 4
  moving_section = moving_section[1:downsample:end, 1:downsample:end]
  fixed_section = fixed_section[1:downsample:end, 1:downsample:end]

  p11 = moving_points + res1
  p22 = fixed_points + res2

  p11 = p11[1:2,:]
  p22 = p22[1:2,:]
  moving_points = moving_points[1:2,:]
  fixed_points = fixed_points[1:2,:]

  p11 = p11/downsample
  p22 = p22/downsample
  moving_points = moving_points/downsample
  fixed_points = fixed_points/downsample
  
  trans = adjust_affine_for_scaling(trans.', downsample).'
  out_img, offset = imwarp(fixed_section, inv(trans))
  println(offset)
  fused, fused_offset = imfuse(moving_section, [0,0], out_img, offset)
  println(fused_offset)

  block_radius = 150
  scalebar = [1; 1; 2*block_radius/downsample; 2*block_radius/downsample]
  draw_vectors(make_isotropic(moving_section), hcat(vcat(moving_points, p11), scalebar))
  draw_vectors(make_isotropic(fixed_section), hcat(vcat(fixed_points, p22), scalebar+100))
  draw_vectors(make_isotropic(fused), hcat(vcat(moving_points.-fused_offset, p11.-fused_offset), scalebar)) # todo: check offset
    #ccp.write_image_from_points(moving_points[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
end

function test3()
  moving_section = "../output_images/(1,1)_montage.tif"
  fixed_section = "../output_images/(1,2)_montage.tif"
  A, meshset = affine_align_sections(moving_section, fixed_section)
  return meshset
end

function points_to_Nx3_matrix(points)
  points = hcat(points...)
  if size(points,1)==2
    points = [points; ones(eltype(points), 1, size(points,2))]'
  end
  return points
end

"""
`AFFINE_ALIGN_IMAGES`

Returns the affine transformation tform, such that p_in_B = p_in_A * tform, where p is 
point coordinates in row vector.
"""
function affine_align_images(moving_img::Array{}, fixed_img::Array{}, 
                                params=PARAMS_PREALIGNMENT; return_points=false)
  # points here are in column vector convention
  if params.scaling_factor != 1.0
    println("Finding points at ", params.scaling_factor, "x")
    s = [params.scaling_factor 0 0; 0 params.scaling_factor 0; 0 0 1]
    moving_img, moving_offset = imwarp(moving_img, s)
    fixed_img, fixed_offset = imwarp(fixed_img, s)
  end
  points = generate_match_points(moving_img, fixed_img, params)
  moving_pointslist, fixed_pointslist = get_block_matches(moving_img, fixed_img, points, params)

  n_matches = length(moving_pointslist)
  println("n_matches: ", n_matches)

  moving_points = points_to_Nx3_matrix(moving_pointslist)
  fixed_points = points_to_Nx3_matrix(fixed_pointslist)

  tform = find_affine(moving_points, fixed_points)
  residualIn2 = moving_points*tform - fixed_points
  rmsIn2 = mean( sum(residualIn2.^2, 1) )^0.5
  fixed_pointsin1 = fixed_points*inv(tform)
  residualIn1 = fixed_pointsin1 - moving_points
  rmsIn1 = mean( sum(residualIn1.^2, 1) )^0.5

  tolerance_ratio = 0.001
  rmsThres = tolerance_ratio * mean([size(moving_img)..., size(fixed_img)...])
  rmsTotal = mean([rmsIn1, rmsIn2].^2)^0.5

  if n_matches < 0.5 * length(points)
    println("WARNING [affine_align_sections]: # of matches is small. n_matches: ", n_matches)
  end
  if  rmsTotal > rmsThres
    println("WARNING [affine_align_sections]: high residual. RMS error: ", rmsTotal)
  end

  tform = adjust_affine_for_scaling(tform, params.scaling_factor)

  if return_points
    return tform, moving_pointslist, fixed_pointslist, residualIn1, residualIn2, rmsIn1, rmsIn2, rmsTotal
  else
    return tform
  end
end


function generate_match_points(moving_img::Array{}, fixed_img::Array{}, 
                                                  params=PARAMS_PREALIGNMENT)
  border_ratio = 0.1
  mesh_length = params.mesh_length * params.scaling_factor
  grid_size = minimum(floor(Int, collect(size(moving_img))/mesh_length))
  block_radius = params.block_size * params.scaling_factor
  search_radius = params.search_r * params.scaling_factor
  overlap = [min(size(moving_img), size(fixed_img))...]
  #overlap = size(moving_img)

  border = round(Int, border_ratio * minimum(overlap))
  combined_radius = search_radius+block_radius
  border = border > combined_radius ? border : combined_radius
  println("border excluded: ", border)

  step = floor(Int, (overlap - 2*border - 1) ./ (grid_size-1))
  points = [[x; y] for x = 1+border:step[1]:overlap[1]-border, 
                                y = 1+border:step[2]:overlap[2]-border]
  points = points[:]
  println("n points: ", size(points))
  return points
end

function get_block_matches(moving_img::Array{}, fixed_img::Array{}, points, 
                          params=PARAMS_PREALIGNMENT)
  moving_points = []
  fixed_points = []
  for pt = points
    offset, r, xc = block_match_at_point(moving_img, pt, fixed_img, pt, 
                                            params.block_size, params.search_r)
    println(pt, offset, r)
    if r >= params.min_r
      push!(moving_points, collect(pt))
      push!(fixed_points, collect(pt+offset))
    end
  end
  return moving_points, fixed_points
end


function adjust_affine_for_scaling(tform, scale)
  s = [scale 0 0; 0 scale 0; 0 0 1]
  return s * tform * inv(s)
end

"""
`FIND_AFFINE` - Calculate left-hand affine transform from two Nx3 point sets

Inputs:
    2D arrays with each column being a point in homogeneous coords.
Returns:
    Affine transformation that produces moving_points*tform = fixed_points


Affine transform of a global position:

* homogeneous coordinates [x, y, 1]  [ax + by + c, dx + ey + f, 1]
* or equivalently [x, y, 1]  [x, y, 1] * tform

        where `tform` = [a d 0;  
                         b e 0;  
                         c f 1]
"""
function find_affine(moving_points, fixed_points)
  return moving_points \ fixed_points
end

function block_match_at_point(A, pointInA, B, pointInB, block_radius, search_r)
# A is the template
# Returns:
# offset: offset from pointInB to pointInA's actual match in B.
# r_max:  correlation value at the match point
# xc:   raw cross-correlation map

  b = block_radius
  B_radius = block_radius + search_r

  A_lower = pointInA - b
  A_upper = pointInA + b
  B_lower = pointInB - B_radius
  B_upper = pointInB + B_radius

  xc = normxcorr2(sub(A, [l:u for (l,u) in zip(A_lower, A_upper)]...),
          sub(B, [l:u for (l,u) in zip(B_lower, B_upper)]...))

  r_max, ind = findmax(xc); 
  i_max, j_max = ind2sub(xc, ind)

  return [i_max-1-search_r; j_max-1-search_r], r_max, xc
  
end

function recompute_affine(meshset::MeshSet)
  params = meshset.params
  moving_pointslist = meshset.meshes[1].nodes
  fixed_pointslist = meshset.matches[1].dst_points
  moving_points = points_to_3xN_matrix(moving_pointslist)
  fixed_points = points_to_3xN_matrix(fixed_pointslist)
  tform = find_affine(moving_points, fixed_points)
  tform = adjust_affine_for_scaling(tform, params.scaling_factor)
  # convert column vector convention to row vector convention
  return tform
end

function affine_align_sections(moving_img_filename::String, 
                                                fixed_img_filename::String, 
                                                params=PARAMS_PREALIGNMENT)
  meshset = makeNewMeshSet(params)
  meshset.N = 2
  moving_mesh = Mesh(moving_img_filename)
  fixed_mesh = Mesh(fixed_img_filename)
  meshset.meshes = [fixed_mesh, moving_mesh]

  if meshset.N != 2
    error("Invalid Arguments")
  end

  moving_img = []
  fixed_img = []
  try
    fixed_img = getUfixed8Image(meshset.meshes[1])
    moving_img = getUfixed8Image(meshset.meshes[2])
  catch
    # mostly for testing purpose where name is the file path
    fixed_img = getUfixed8Image(meshset.meshes[1].name)
    moving_img = getUfixed8Image(meshset.meshes[2].name)
  end

  # Align
  tform, moving_points, fixed_points, residualIn1, residualIn2, rmsIn1, rmsIn2, rmsTotal = 
    affine_align_images(moving_img, fixed_img, params; return_points=true)
  println(tform)

  # Update matches object in meshset
  meshset.meshes[1].nodes = moving_points
  src_points_indices = 1:length(moving_points)
  dst_points = fixed_points
  matches = Matches(meshset.meshes[2].index, meshset.meshes[1].index, 
    length(moving_points), src_points_indices, dst_points, [], [], [])
  meshset.matches = [matches]
  println("Saving JLD for ", moving_img_filename[1:end-4], 
                                  " to ", fixed_img_filename[1:end-4])
  save(meshset)  

  warped_index = meshset.meshes[2].index
  # Write overlay thumbnail for review
  scale = 0.125
  s = [scale 0 0; 0 scale 0; 0 0 1]
  @time imgA, A_offset = imwarp(fixed_img, s)
  @time imgB, B_offset = imwarp(moving_img, tform*s)
  O, O_bb = imfuse(imgA, A_offset, imgB, B_offset)
  thumbnail_fn = string(join(warped_index[1:2], ","), "_prealigned_thumbnail.jpg")
  println("Writing ", thumbnail_fn)
  imwrite(restrict(O), joinpath(PREALIGNED_DIR, "review", thumbnail_fn))

  # Render image
  log_path = joinpath(PREALIGNED_DIR, "prealigned_offsets.txt")
  if !isfile(log_path)
    f = open(log_path, "w")
    close(f)
  end
  log_file = open(log_path, "a")
  warped_fn = string(join(warped_index[1:2], ","), "_prealigned.tif")
  println("Rendering ", warped_fn)
  @time warped_img, warped_offset = imwarp(moving_img, tform)
  println("Writing ", warped_fn)
  @time imwrite(warped_img, joinpath(PREALIGNED_DIR, warped_fn))
  log_line = join((warped_fn, warped_offset[1], warped_offset[2], 
                      size(warped_img,1), size(warped_img,2)), " ")
  write(log_file, log_line, "\n")
  close(log_file)
end

function compute_propogated_transform(index::Index)
  index = (index[1:2]..., 0, 0)
  filenames = filter(x -> x[end-2:end] == "jld", readdir(PREALIGNED_DIR))
  indices = [(parseName(x), x) for x in filenames]
  sort!(indices)
  if indices[1][1] == ones(Int, 4)
    error("Could not parse JLD filename to index: ", indices[1][2])
  end
  for k in 2:length(indices)
    if indices[k][1] > index
      break
    end
    if !isAdjacent(indices[k-1][1], indices[k][1])
      error("Missing section between ", indices[k-1][1], " and ", indices[k][1])
    end
  end
  T = diagm(ones(3))
  for k in 2:length(indices)
    ind, fn = indices[k]
    if ind > index
      break
    end
    meshset = load(joinpath(PREALIGNED_DIR, fn))["MeshSet"]
    A = recompute_affine(meshset)
    # row vector homogeneous point convention
    T *= A
  end
  return T
end

function prealignment()
    @time prealign_directory()
    @time render_prealignment_for_directory()
end
