module ImageRegistration

# using Images
# using Cairo

include("mesh.jl")
include("convolve.jl")
include("blockmatch.jl")
include("incidence_to_triangles.jl")
include("registration.jl")
include("transforms.jl")
include("boundingbox.jl")
include("imwarp.jl")
include("meshwarp.jl")
include("geometry.jl")
# include("draw.jl")

export
  #types
  BoundingBox,
  Mesh,
  Matches,
  # core function
  default_params,
  get_max_xc_vector,
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
  fillpoly!,
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
  bb_to_sz,
  +,
  -,
  scale_bb,
  translate_bb,
  slice_to_bb,
  bb_to_slice,
  get_area,
  get_rect,
  get_bounds,
  get_size,
  get_offset,
  intersects,
  pt_in_poly,
  pt_on_line_segment,
  poly_contains_poly,
  poly_intersects,
  cross2,
  clip_polygon,
  pt_on_left,
  line_to_vector,
  find_line_segment_intersection,
  find_line_intersection,
  poly_area,
  load_uint8_img,
  load_test_images,
  make_rotation_matrix,
  make_translation_matrix,
  make_scale_matrix
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
