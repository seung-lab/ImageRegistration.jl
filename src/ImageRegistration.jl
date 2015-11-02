module ImageRegistration

using Images
using FixedPointNumbers
using ImageView
using Colors

include("mesh.jl")
include("convolve.jl")
include("blockmatch.jl")
include("incidence_to_triangles.jl")
include("registration.jl")
include("transforms.jl")
include("boundingbox.jl")
include("imwarp.jl")
include("meshwarp.jl")
include("visualize.jl")

export
  #types
  BoundingBox,
  Mesh,
  Matches,
  # core function
  default_params,
  blockmatch,
  matches2mesh,
  create_mesh,
  calculate_affine,
  calculate_rigid,
  calculate_translation,
  imwarp,
  meshwarp,
  normxcorr2,
  incidence_to_triangles,
  incidence_to_dict,
  dict_to_triangles,
  warp_pts,
  # auxilliary functions
  make_isotropic,
  draw_mesh,
  draw_vectors,
  find_mesh_bb,
  tform_bb,
  snap_bb,
  sz2bb,
  bb2pts,
  +,
  -,
  load_ufixed8_img,
  load_test_images

end