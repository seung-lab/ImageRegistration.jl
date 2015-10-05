module ImageRegistration

using Images
using FixedPointNumbers
using ImageView

include("mesh.jl")
include("convolve.jl")
include("blockmatch.jl")
include("incidence2triangles.jl")
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
  incidence2triangles,
  incidence2dict,
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
  -

end