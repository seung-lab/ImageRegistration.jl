src_nodes = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
dst_nodes = [5.0 1.0;
			6.0 3.0;
			7.0 3.0]
# incidence = [1 -1 0;
# 			1 0 -1;
# 			0 1 -1]
edges = spzeros(Int64, 3, 3)
edges[1:2, 1] = [1, -1]
edges[1, 2] = 1
edges[3, 2] = -1
edges[2:3, 3] = [1, -1]


mesh = Mesh(src_nodes, dst_nodes, edges)
mesh_indices = [1, 2]
matches = Matches(mesh_indices, src_nodes[1:2,:], dst_nodes[1:2,:])
new_mesh = matches_to_mesh(matches, mesh)
true_mesh = Mesh(src_nodes[1:2,:], dst_nodes[1:2,:], edges[1:2, 1])
@test new_mesh.src_nodes == true_mesh.src_nodes
@test new_mesh.dst_nodes == true_mesh.dst_nodes
@test new_mesh.edges == true_mesh.edges