# """
# `DRAW_MESH` - Display mesh on image

# ```
# an_lines = draw_mesh(imgc, img2, nodes, node_dict, color=RGB(1,1,1))  
# ```

# * imgc: ImageCanvas object
# * img2: ImageSliced2d object 
# * nodes: Nx2 array of mesh node positions
# * node_dict: dictionary of edges, indexed by node, containing connected nodes
# * color: optional color argument for mesh lines
# * an_lines: annotation object for the lines describing the mesh edges
# """ 
# function draw_mesh(imgc, img2, nodes, node_dict, color=RGB(1,1,1))  
#   lines = Array(Float64, 4, 0)
#   for k in sort(collect(keys(node_dict)))
#       for v in node_dict[k]
#           a = reverse(vec(nodes[k,:]))
#           b = reverse(vec(nodes[v,:]))
#           lines = hcat(lines, vcat(a, b))
#       end
#   end
#   an_lines = annotate!(imgc, img2, AnnotationLines(lines, color=color, 
#                                                           coord_order="yyxx"))
#   return an_lines
# end

# function draw_mesh(imgc, img2, mesh::Mesh)
#   node_dict = incidence_to_dict(mesh.edges)
#   an_src = draw_mesh(imgc, img2, mesh.src_nodes, node_dict, RGB(0,1,0))
#   an_dst = draw_mesh(imgc, img2, mesh.dst_nodes, node_dict, RGB(1,0,0))
#   return an_src, an_dst
# end

# """
# `DRAW_VECTORS` - Annotate ImageView with vectors (point and line segment)

# ```
# an_points, an_vectors = draw_vectors(imgc, img2, vectors, pt_color=RGB(0,0,1), 
#                                                         vec_color=RGB(1,0,1), 
#                                                         k=10)
# ```
# * imgc: ImageCanvas object
# * img2: ImageSliced2d object 
# * vectors: Nx4 array of vector start and end points
# * pt_color: optional color argument for points
# * vec_color: optional color argument for vectors
# * k: factor to scale the vectors by for easier viewing
# * an_points: annotation object for the points
# * an_vectors: annotation object for the vectors
# """
# function draw_vectors(imgc, img2, vectors, pt_color=RGB(0,0,1), 
#                                               vec_color=RGB(1,0,1), k=1)
#   vectors = [vectors[:,2]'; 
#               vectors[:,1]'; 
#               (vectors[:,4]'-vectors[:,2]')*k + vectors[:,2]'; 
#               (vectors[:,3]'-vectors[:,1]')*k + vectors[:,1]']
#   an_vectors = annotate!(imgc, img2, AnnotationLines(vectors, color=vec_color, 
#                                           coord_order="xxyy", linewidth=3))
#   an_points = annotate!(imgc, img2, AnnotationPoints(vectors[1:2,:], 
#                                                   color=pt_color, shape='*'))
#   return an_points, an_vectors
# end

# function draw_vectors(imgc, img2, matches::Matches)
#   vectors = hcat(matches.src_points, matches.dst_points)
#   return draw_vectors(imgc, img2, vectors)
# end

# function draw_vectors(imgc, img2, mesh::Mesh)
#   vectors = hcat(mesh.src_nodes, mesh.dst_nodes)
#   return draw_vectors(imgc, img2, vectors)
# end


