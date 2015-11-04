module ImageRegistration

using Images
using Cairo

include("mesh.jl")
include("convolve.jl")
include("blockmatch.jl")
include("incidence_to_triangles.jl")
include("registration.jl")
include("transforms.jl")
include("boundingbox.jl")
include("imwarp.jl")
include("meshwarp.jl")
include("draw.jl")

export
  #types
  BoundingBox,
  Mesh,
  Matches,
  # core function
  default_params,
  blockmatch,
  matches_to_mesh,
  create_mesh,
  calculate_affine,
  calculate_rigid,
  calculate_translation,
  imwarp,
  meshwarp,
  normxcorr2,
  find_zero_indices,
  incidence_to_triangles,
  incidence_to_dict,
  dict_to_triangles,
  warp_pts,
  # auxilliary functions
  find_mesh_bb,
  tform_bb,
  snap_bb,
  sz_to_bb,
  bb_to_pts,
  +,
  -,
  load_uint8_img,
  load_test_images,
  # create_drawing,
  # create_contex,
  # draw_vectors,
  # draw_line,
  # draw_points,
  # draw_point,
  # draw_text,
  # draw_indices,
  # get_drawing,
  # draw_reference

end
