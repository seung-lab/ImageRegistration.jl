"""
`DEFAULT_PARAMS` - Create dictionary of default parameters for blockmatching
"""
function default_params()
  return Dict(  "mesh_dist" => 750, 
                "block_size" => 200, 
                "search_r" => 150, 
                "min_r" => 0.25)
end

"""
`BLOCKMATCH` - Build mesh to produce nodes, and blockmatch between two images.
"""
function blockmatch(src_img, dst_img; src_offset=[0,0], dst_offset=[0,0], 
                                                      params=default_params())
  mesh = create_mesh(src_img, src_offset, params["mesh_dist"])
  matches = blockmatch(mesh.src_nodes, src_img, dst_img, 
                                        src_offset, dst_offset, params)
  return mesh, matches
end

"""
`MATCHES_TO_MESH` - Using matches and old mesh, create new mesh with just nodes
  related to the matches.

```
new_mesh = matches_to_mesh(matches::Matches, mesh::Mesh)
```

* matches: set of Matches
* mesh: Mesh from which the Matches were generated
* new_mesh: adjusted Mesh using matches.mesh_indices to filter only the nodes
    where matches were found. The edge-node incidence matrix will also be
    adjusted. 
"""
function matches_to_mesh(matches::Matches, mesh::Mesh)
  assert(length(matches.mesh_indices) == size(matches.dst_points, 1))
  new_mesh = Mesh()
  mask = matches.mesh_indices
  new_mesh.src_nodes = mesh.src_nodes[mask, :]
  new_mesh.dst_nodes = matches.dst_points
  edges = mesh.edges[mask, :]
  non_orphans = find_zero_indices(sum(edges, 1))
  new_mesh.edges = edges[:, non_orphans]
  return new_mesh
end

function meshwarp(img, mesh::Mesh, offset=[0,0])
  triangles = incidence_to_triangles(mesh.edges)
  return @time meshwarp(img, mesh.src_nodes, mesh.dst_nodes, triangles, offset)
end

"""
Transform Nx2 pts by 3x3 tform matrix
"""
function warp_pts(tform, pts)
    pts = hcat(pts, ones(size(pts,1)))
    tpts = pts * tform
    return tpts[:,1:2]
end

"""
Create array of UInt8 type from from image path
"""
function load_uint8_image(path::AbstractString)
  img = Images.load(path)
  return reinterpret(UInt8, Images.data(img)[:,:,1])'
end

"""
Load test images
"""
function load_test_images()
  dir = dirname(dirname(@__FILE__))
  pathA = joinpath(dir, "test", "test_images", "imgA.tif")
  pathB = joinpath(dir, "test", "test_images", "imgB.tif")
  return load_uint8_image(pathA), load_uint8_image(pathB)
end
