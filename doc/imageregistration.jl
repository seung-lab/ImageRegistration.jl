export Point, Points

global const NO_MATCH = [0; 0; -1];
global const NO_TRIANGLE = (0, 0, 0);
global const NO_RANGE = (0:0, 0:0);

typealias Point Array{Float64, 1};        # [i, j]
typealias Points Array{Point, 1};       # array of points

type Mesh
  offset::Point
  src_nodes::Points  # 2-by-n dense matrix of nodes, each column stores [i-coordinate, j-coordinate] in global coordinates
  dst_nodes::Points
  edges::Edges       # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
end

type Matches
  src_points::Points
  dst_points::Points
end

function default_params()
  return Dict("scaling_factor" => 1.0, 
                "mesh_length" => 750, 
                "gaussian_sigma" => 1.0,
                "block_size" => 200, 
                "search_r" => 150, 
                "min_r" => 0.25)
end

function blockmatch(imgA, imgB, offsetA=[0,0], offsetB=[0,0], params=default_params())
  if params["gaussian_sigma"] != 0
    sigma = params["gaussian_sigma"]
    imgA = imfilter_gaussian(imgA, [sigma, sigma])
    imgB = imfilter_gaussian(imgB, [sigma, sigma])
  end
  meshA = create_mesh(imgA, offsetA, params)
end

function create_mesh(img, offset, dist)
  nodes = Points(0)
  n, m = floor(Int64, collect(size(img))/dist)
  for i in 1:n
    for j in 1:m
      push!(nodes, [i,j]*dist+offset)
    end
  end
  edges = DelaunayTessellation()
  push!(edges, nodes)
  return Mesh(offset, nodes, nodes, edges)
end

function get_triangles(mesh::Mesh)
end

function matches2mesh(matches, offset=[0,0])
  edges = DelaunayTessellation()
  push!(edges, matches.src_points)
  return Mesh(offset, matches.src_points, matches.dst_points, edges)
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