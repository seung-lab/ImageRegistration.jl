using HDF5
using JLD
using Images
using ImageView
using PiecewiseAffineTransforms
using Color
using FixedPointNumbers

include("incidence2triangles.jl")
include("piecewiseaffine_warp.jl")

const BUCKET = "/usr/people/tmacrina/seungmount/research/tommy/Julimaps"

xy2ij(shape, height) = [height .- shape[:, 2] shape[:, 1]]
rawdata(img) = convert(Array{Float64, 2}, data(separate(img)))

function demo_mesh()
    # file = "8000x8000_100px.jld"
    file = "r4c2-r4c3.jld"
    path = joinpath(BUCKET, file)
    d = load(path)
    v = convert(Array{Float64,2}, d["v"]) # nodes
    e = convert(Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, d["e"]) # edges
    k = d["k"] # spring constants
    l = d["l"] # resting lengths
    n = d["n"] # number of edges
    m = d["m"] # number of nodes
    incidence = e
    nodes = v
    # Display white mesh on black background
    node_dict = incidence2dict(incidence)
    triangles = dict2triangles(node_dict)
    sz = round(Int,maximum(nodes[1,:]))+10, round(Int,maximum(nodes[2,:]))+10
    img = Image(zeros(Bool, sz...)) # Needs to be Float
    shape = xy2ij(nodes', size(img, 1))

    println("Original images and shapes")
    # triplot(img, shape, triangles)
    # plot_edges2(img, shape, node_dict)
    plot_mesh(img, nodes', node_dict)
end

function plot_mesh(img, shape, node_dict)
    imgc, img2 = view(img)
    lines = Array(Float64, 4, 0)
    for k in sort(collect(keys(node_dict)))
        for v in node_dict[k]
            a = reverse(vec(shape[k,:]))
            b = reverse(vec(shape[v,:]))
            lines = hcat(lines, vcat(a, b))
        end
    end
    annotate!(imgc, img2, AnnotationLines(lines, coord_order="xxyy"))
end

function load_mesh()
    file = "solvedMesh.jld"
    path = joinpath(BUCKET, file)
    d = load(path)

    v = convert(Array{Float64,2}, d["nodes"]) # nodes
    vt = convert(Array{Float64,2}, d["nodes_t"]) # transformed nodes
    e = convert(Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, d["edges"]) # edges
    # k = d["k"] # spring constants
    # l = d["l"] # resting lengths

    # XWorldLimits: [2.1907e+04 2.9881e+04]
    # YWorldLimits: [3.6430e+04 4.4393e+04]
    initial_nodes = v;
    final_nodes = vt;
    initial_nodes[1,:] = v[1,:] - 21906
    initial_nodes[2,:] = v[2,:] - 36429
    final_nodes[1,:] = vt[1,:] - 21906
    final_nodes[2,:] = vt[2,:] - 36429    

    tileA = rawdata(imread(joinpath(BUCKET, "test_images", "Tile_r4-c2_S2-W001_sec20.tif")))
    tileB = rawdata(imread(joinpath(BUCKET, "test_images", "Tile_r4-c3_S2-W001_sec20.tif")))

    high = Int64(ceil(max(maximum(initial_nodes), maximum(final_nodes))))
    szA = size(tileA)
    tileA_padded = vcat(tileA, zeros(high-szA[1], szA[2]))
    szA = size(tileA_padded)
    tileA_padded = hcat(tileA_padded, zeros(szA[1], high-szA[2]))

    src_shape = xy2ij(initial_nodes', size(tileA_padded, 1))
    dst_shape = xy2ij(final_nodes', size(tileA_padded, 1))

    tileA_warped = warp_piecewiseaffine(tileA_padded, e, src_shape, dst_shape)     
    return tileA, tileA_warped
end

function warp_piecewiseaffine(img, incidence, initial_nodes, final_nodes)
    node_dict = incidence2dict(incidence)
    triangles = dict2triangles(node_dict)
    warped = pa_warp2(img, initial_nodes, final_nodes, triangles)
    plot_mesh(img, initial_nodes, node_dict)
    plot_mesh(warped, final_nodes, node_dict)
    return warped
end

function demo_warp2()
    img = imread(joinpath(BUCKET, "test_images", "turtle.jpg"))
    img = convert(Array{Float64, 3}, data(separate(img)))
    initial_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    320.0 290.0]
    final_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    400.0 460.0]
    incidence = [1 1 1 0 0 0 0 0;
                -1 0 0 1 1 0 0 0;
                0 0 0 -1 0 1 1 0;
                0 -1 0 0 0 0 -1 1;
                0 0 -1 0 -1 -1 0 -1]
    triangles = [1 2 5;
                1 4 5;
                2 3 5;
                3 4 5];
    src_shape = xy2ij(initial_nodes, size(img, 1))
    dst_shape = xy2ij(final_nodes, size(img, 1))
    node_dict = incidence2dict(incidence)
    plot_mesh(img, src_shape, node_dict)

    warp = pa_warp2(img, src_shape, dst_shape, triangles)
    plot_mesh(warp, dst_shape, node_dict)
end

function demo_warp3()
    img = imread(joinpath(BUCKET, "test_images", "turtle.jpg"))
    img = convert(Array{Float64, 3}, data(separate(img)))
    initial_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    320.0 290.0]
    final_nodes = [20.0 20.0;
                    620.0 20.0;
                    620.0 560.0;
                    20.0 560.0;
                    400.0 460.0]
    incidence = [1 1 1 0 0 0 0 0;
                -1 0 0 1 1 0 0 0;
                0 0 0 -1 0 1 1 0;
                0 -1 0 0 0 0 -1 1;
                0 0 -1 0 -1 -1 0 -1]
    triangles = [1 2 5;
                1 4 5;
                2 3 5;
                3 4 5];
    src_shape = xy2ij(initial_nodes, size(img, 1))
    dst_shape = xy2ij(final_nodes, size(img, 1))
    node_dict = incidence2dict(incidence)
    plot_mesh(img, src_shape, node_dict)

    warp = pa_warp3(img, src_shape, dst_shape, triangles)
    plot_mesh(warp, dst_shape, node_dict)
end