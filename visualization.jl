# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150810
#
# Visualize the meshes and displacement vectors with image annotations and 
# color plots.

using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "") # include current directory
@pyimport create_color_plots as ccp

function write_vector_color_plot(mesh_path, write_name, offset=[0,0])
# Write color plot to describe displacement vectors
# Args:
#   mesh_path: string to the JLD file
#   write_name: string to save color plot image
#   offset: 2-element array noting global offset of the images to adjust nodes
# Returns:
#   None
    v, vt, e = load_mesh(mesh_path) # original verts, final verts, incidence
    v[1,:] = v[1,:] - offset[1]
    v[2,:] = v[2,:] - offset[2]
    vt[1,:] = vt[1,:] - offset[1]
    vt[2,:] = vt[2,:] - offset[2]
    ccp.write_image_from_points(v', vt', write_name)
end

function draw_mesh(imgc, img2, nodes, node_dict, color=RGB(1,1,1))
# Display mesh on image
# Args:
#   imgc: ImageCanvas object
#   img2: ImageSliced2d object  
#   nodes: Nx2 array of mesh node positions
#   node_dict: dictionary of edges, indexed by node, containing connected nodes
# Returns:
#   imgc: ImageCanvas object
#   img2: ImageSliced2d object    
    lines = Array(Float64, 4, 0)
    for k in sort(collect(keys(node_dict)))
        for v in node_dict[k]
            a = reverse(vec(nodes[k,:]))
            b = reverse(vec(nodes[v,:]))
            lines = hcat(lines, vcat(a, b))
        end
    end
    annotate!(imgc, img2, AnnotationLines(lines, color=color, coord_order="xxyy"))
    return imgc, img2
end

function draw_mesh(img, nodes, node_dict)
    imgc, img2 = view(img)
    return draw_mesh(imgc, img2, nodes, node_dict)
end

function draw_vectors(img2, imgc, vectors, color=RGB(0,0,1))
# Display match displacement vectors on images
# Args:
#   imgc: ImageCanvas object
#   img2: ImageSliced2d object 
#   nodes: 4xN array of vector start and end points
# Returns:
#   imgc: ImageCanvas object
#   img2: ImageSliced2d object  
    lines = Array(Float64, 4, 0)
    points = Array(Float64, 2, 0)
    for j in 1:size(vectors, 2)
        points = hcat(points, vectors[1:2, j])
        lines = hcat(lines, vectors[:, j])
    end
    annotate!(imgc, img2, AnnotationPoints(points, color=color))
    annotate!(imgc, img2, AnnotationLines(lines, color=color, coord_order="xxyy"))
    return imgc, img2
end    

function draw_vectors(img, vectors)
    imgc, img2 = view(img)
    return draw_vectors(img2, imgc, vectors)
end

function draw_vectors(img, start_pts, end_pts)
    imgc, img2 = view(img)
    vectors = vcat(start_pts, end_pts)
    return draw_vectors(img2, imgc, vectors)
end

function draw_points(img2, imgc, pts, color=RGB(0,0,1))
# Display match displacement vectors on images
# Args:
#   imgc: ImageCanvas object
#   img2: ImageSliced2d object 
#   nodes: 2xN array of points
# Returns:
#   imgc: ImageCanvas object
#   img2: ImageSliced2d object  
    points = Array(Float64, 2, 0)
    for j in 1:size(pts, 2)
        points = hcat(points, pts[1:2, j])
    end
    annotate!(imgc, img2, AnnotationPoints(points, color=color))
    return imgc, img2
end 

function draw_points(img, pts)
    imgc, img2 = view(img)
    return draw_points(img2, imgc, pts)
end  

function draw_imfuse_meshes(Oc, O2, dst_nodes_A, SR_A, dst_nodes_B, SR_B)
# Incomplete
    SR_A = [0, 0]
    # SR_B = [7184.9, -178.7780] # NEED INTERPOLATION!
    SR_B = [7185, -179]
    SR_C = SR_B - SR_A
    if SR_C[1] > 0
        dst_nodes_B[:, 2] += SR_C[1]
    elseif SR_C[1] < 0
        dst_nodes_A[:, 2] -= SR_C[1]
    end 
    if SR_C[2] > 0
        dst_nodes_B[:, 1] += SR_C[2]
    elseif SR_C[2] < 0
        dst_nodes_A[:, 1] -= SR_C[2]
    end
    imgc, img2 = draw_mesh(imgc, img2, dst_nodes_B, node_dict_B, RGB(0,0,1))
    imgc, img2 = draw_mesh(imgc, img2, dst_nodes_A, node_dict_A, RGB(1,0,1))
end

function demo_color_plot()
# Write demo vector color plot to file
    mesh_path = joinpath(BUCKET, "EM_Images", "r4c3_solved.jld")
    write_name = "Tile_r4-c3_S2-W001_sec20_hexplot2.png"
    offset = [29090, 36251]
    write_vector_color_plot(mesh_path, write_name, offset)
end

function demo_mesh()
# Load a mesh and display it on a black background
    mesh_path = joinpath(BUCKET, "test_images", "solvedMesh.jld")
    offset = [21906, 36429]
    v, vt, e = load_mesh(mesh_path, offset)
    incidence = e
    initial_nodes = v
    nodes = xy2yx(initial_nodes')
    node_dict = incidence2dict(incidence)
    sz = round(Int,maximum(nodes[:,1]))+10, round(Int,maximum(nodes[:,2]))+10
    img = zeros(Bool, sz...)

    println("Original images and shapes")
    draw_mesh(make_isotropic(img), nodes, node_dict)
end

function demo_vectors()
# Load matches and display them on a black background
    vector_path = joinpath(BUCKET, "test_images", "vectors.jld")
end
