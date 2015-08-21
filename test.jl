using Base.Test

include("render.jl")
include("meshwarp.jl")
include("pt2triangle.jl")
include("BoundingBox.jl")

function test()
  test_verts2triangle()
  test_xy2barycentric()
  test_pt2triangle()
  test_imfuse()
  test_padimage()
  test_mesh_warp()
  test_bb_operations()
  test_find_mesh_bb()
  test_tform_bb()
  test_snap_bb()
end