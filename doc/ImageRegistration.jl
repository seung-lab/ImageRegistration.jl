export Point, Points

global const NO_MATCH = [0; 0; -1];
global const NO_TRIANGLE = (0, 0, 0);
global const NO_RANGE = (0:0, 0:0);

typealias Point Array{Float64, 1};        # [i, j]
typealias Points Array{Point, 1};       # array of points
typealias Edges SparseMatrixCSC{Float64, Int64}     # sparse array for edges - columns represent edges and the rows represent the nodes

type Mesh
  offset::Point
  src_nodes::Points  # 2-by-n dense matrix of nodes, each column stores [i-coordinate, j-coordinate] in global coordinates
  dst_nodes::Points
  edges::Edges       # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
end

Mesh() = Mesh([0,0], [], [], spzeros(0,0))

type Matches
  mesh_indices
  src_points::Points
  dst_points::Points
end

function default_params()
  return Dict("scaling_factor" => 1.0, 
                "mesh_dist" => 750, 
                "gaussian_sigma" => 1.0,
                "block_size" => 200, 
                "search_r" => 150, 
                "min_r" => 0.25)
end

function blockmatch(src_img, dst_img, src_offset=[0,0], dst_offset=[0,0], 
                                                      params=default_params())
  if params["gaussian_sigma"] != 0
    sigma = params["gaussian_sigma"]
    src_img = imfilter_gaussian(src_img, [sigma, sigma])
    dst_img = imfilter_gaussian(dst_img, [sigma, sigma])
  end
  mesh = create_mesh(src_img, src_offset, params["mesh_dist"])
  return get_blockmatches(mesh, src_img, dst_img, 
                                  src_offset, dst_offset, params)
end

function get_triangles(mesh::Mesh)
  node_dict = incidence2dict(mesh.edges)
  return dict2triangles(node_dict)
end

function matches2mesh(matches, mesh)
  new_mesh = Mesh()
  mask = matches.mesh_indices
  new_mesh.offset = mesh.offset
  new_mesh.src_nodes = mesh.src_nodes[mask]
  new_mesh.dst_nodes = matches.dst_points
  new_mesh.edges = mesh.edges[mask, :]
  return new_mesh
end

function calculate_translation(matches)
  return calculate_translation(matches.src_points, matches.dst_points)
end

function calculate_rigid(matches)
  return calculate_rigid(matches.src_points, matches.dst_points)
end

function calculate_affine(matches)
  return calculate_affine(matches.src_points, matches.dst_points)
end

function count_nodes(mesh::Mesh)
  return length(src_nodes)
end