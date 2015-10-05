"""
`DEFAULT_PARAMS` - Create dictionary of default parameters for blockmatching
"""
function default_params()
  return Dict("scaling_factor" => 1.0, 
                "mesh_dist" => 750, 
                "block_size" => 200, 
                "search_r" => 150, 
                "min_r" => 0.25)
end

"""
`BLOCKMATCH` - Build mesh to produce nodes, and blockmatch between two images.
"""
function blockmatch(src_img, dst_img, src_offset=[0,0], dst_offset=[0,0], 
                                                      params=default_params())
  mesh = create_mesh(src_img, src_offset, params["mesh_dist"])
  matches = blockmatch(mesh.src_nodes, src_img, dst_img, 
                                        src_offset, dst_offset, params)
  return mesh, matches
end

"""
`MATCHES2MESH` - Using matches and old mesh, create new mesh with just nodes
  related to the matches.

```
new_mesh = matches2mesh(matches::Matches, mesh::Mesh)
```

* matches: set of Matches
* mesh: Mesh from which the Matches were generated
* new_mesh: adjusted Mesh using matches.mesh_indices to filter only the nodes
    where matches were found. The edge-node incidence matrix will also be
    adjusted. 
"""
function matches2mesh(matches::Matches, mesh::Mesh)
  assert(length(matches.mesh_indices) == size(matches.dst_points, 1))
  new_mesh = Mesh()
  mask = matches.mesh_indices
  new_mesh.src_nodes = mesh.src_nodes[mask, :]
  new_mesh.dst_nodes = matches.dst_points
  edges = mesh.edges[:, mask]
  non_orphans = find_zero_indices(sum(edges, 2))
  new_mesh.edges = edges[non_orphans,:]
  return new_mesh
end

function meshwarp(img, mesh::Mesh, offset=[0,0])
  triangles = incidence2triangles(mesh.edges)
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