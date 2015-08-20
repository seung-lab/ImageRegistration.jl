# Adapted from PiecewiseAffineTransforms.jl
# https://github.com/dfdx/PiecewiseAffineTransforms.jl

include("BoundingBox.jl")

function source_point_mat(X, Y)
    x1, x2, x3 = X
    y1, y2, y3 = Y    
    F = [x1 x2 x3; y1 y2 y3; 1 1 1]
    F
end

function target_point_mat(U, V)
    u1, u2, u3 = U
    v1, v2, v3 = V
    T = [u1 u2 u3; v1 v2 v3]
    T
end

function affine_params(X, Y, U, V)
    F = source_point_mat(X, Y)
    T = target_point_mat(U, V)
    M = T * inv(F)
    M
end

function meshwarp{N}(img::Array{Float64, N},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                    trigs::Matrix{Int64}, bb = BoundingBox(), interp = true)
    println(dst)
    wbb = snap_bb(find_mesh_bb(dst))
    println(wbb)
    low, high = bounds2padding(size(img), minsandmax(wbb)...)
    println(low, high)
    img = padimage(img, low..., high...)
    src = src .+ low'
    dst = dst .+ low'
    wbb = BoundingBox(bb.x-low[1], bb.y-low[2], size(img)...)
    warped = zeros(eltype(img), size(img))

    # wbb = snap_bb(find_mesh_bb(dst))
    # warped = similar(img, Int64(bb.h), Int64(bb.w))

    if interp
        println("w interpolation")
    else
        println("w/o interpolation")
    end
    for t=1:size(trigs, 1)    
        tr = squeeze(trigs[t, :], 1)
        Y = src[tr, 1]
        X = src[tr, 2]
        
        V = dst[tr, 1]
        U = dst[tr, 2]
               
        # warp parameters from target (U, V) to source (X, Y)
        M = affine_params(U, V, X, Y)
        vs, us = poly2source(U, V)
        
        # for every pixel in target triangle we find corresponding pixel in source
        # and copy its value
        if interp
            for i=1:length(vs)
                u, v = us[i], vs[i]
                # x, y = M * [u, v, 1]
                x, y = M[1,1]*u + M[1,2]*v + M[1,3], M[2,1]*u + M[2,2]*v + M[2,3]

                fx, fy = floor(Int64, x), floor(Int64, y)
                wx, wy = x-fx, y-fy

                if 1 <= fy && fy+1 <= size(img, 1) && 1 <= fx && fx+1 <= size(img, 2)
                    # Expansion of p = [1-wy wy] * img[fy:fy+1, fx:fx+1] * [1-wx; wx]
                    p = ((1-wy)*img[fy,fx] + wy*img[fy+1,fx]) * (1-wx) + ((1-wy)*img[fy,fx+1] + wy*img[fy+1,fx+1]) * wx
                    warped[v, u] = p
                end
            end
        else
            for i=1:length(vs)
                u, v = us[i], vs[i]
                # x, y = M * [u, v, 1]
                x = round(Int64, M[1,1]*u + M[1,2]*v + M[1,3])
                y = round(Int64, M[2,1]*u + M[2,2]*v + M[2,3])
                if 1 <= y && y <= size(img, 1) && 1 <= x && x <= size(img, 2)
                    warped[v, u] = img[y, x]
                end
            end
        end
        
    end
    return warped, wbb
end

function poly2source(px, py)
# Run fillpoly2 and findn over the basic bounding box of a given triangle
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    top, bottom = floor(Int64,minimum(py)), ceil(Int64,maximum(py))
    mask = zeros(Bool, bottom-top+1, right-left+1)
    fill = fillpoly2!(mask, px-left+1, py-top+1, true)
    vs, us = findn(fill)
    vs += top-1
    us += left-1
    return vs, us
end

function fillpoly2!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T)
    @assert length(px) == length(py)    
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    for x=left:right     
        ys = Set{Int64}()
        j = length(px)
        for i=1:length(px)            
            if (px[i] <= x && x <= px[j]) || (px[j] <= x && x <= px[i])
                # special case: adding the whole cut to ys                            
                if px[i] == px[j]
                    push!(ys, ceil(Int64, py[i]))
                    push!(ys, ceil(Int64, py[j]))
                else
                    y = py[i] + (x - px[i]) / (px[j] - px[i]) * (py[j] - py[i])
                    push!(ys, ceil(Int64, y))
                end            
            end
            j = i
        end
        ys = sort([y for y in ys])
        # if there's an odd number of intersection points, add one imeginary point
        if length(ys) % 2 == 1
            push!(ys, ys[end])
        end
        for i=1:2:length(ys)           
            M[ys[i]:ys[i+1], x] = value  # <-- bounds error here!
        end
    end
    return M
end

function warp_pts(affine, pts)
    pts = hcat(pts, ones(size(pts,1)))
    tpts = pts * affine
    return tpts[:,1:2]
end

function test_mesh_warp()
    img = reshape(float(collect(1:121).%2), 11, 11) # 11x11 checkerboard
    # incidence = [1 1 0;
    #             -1 0 1;
    #             0 -1 -1]
    triangles = [1 2 3]

    src = [2.0 2.0;
            10.0 2.0;
            6.0 10.0]
    dst = [4.0 2.0;
            8.0 2.0;
            6.0 10.0]
    offset = [0, 0]
    img_warped, bb = meshwarp(img, src, dst, triangles)
    @test bb == BoundingBox(0,0,11,12)

    offset = [10, 10]
    iimg_warped, bb = meshwarp(img, src, dst, triangles, BoundingBox(10, 10, 11, 11))
    @test bb == BoundingBox(10,10,11,12)   

    # tform = [cos(pi/6) -sin(pi/6) 0;
    #         sin(pi/6) cos(pi/6) 0;
    #         0 0 1]
    # src = [-16.0 0.0;
    #         27.0 0.0;
    #         6.0 40.0]    
    # dst = warp_pts(tform, src)
    # offset = [0, 0]
    # iimg_warped, bb = meshwarp(img, src, dst, triangles)
    # @test bb == BoundingBox(0,0,11,11)

    tform = [1 0 0;
            0 1 0;
            10 10 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]   
    dst = warp_pts(tform, src)
    iimg_warped, bb = meshwarp(img, src, dst, triangles)
    @test bb == BoundingBox(0,0,20,20) 

    tform = [1 0 0;
            0 1 0;
            -4 -4 1]
    src = [1.0 1.0;
            1.0 4.0;
            4.0 0.0]   
    dst = warp_pts(tform, src)
    iimg_warped, bb = meshwarp(img, src, dst, triangles)
    @test bb == BoundingBox(-4,-5,16,15) 

    src = [0.0 0.0;
            0.0 10.0;
            20.0 20.0]  
    dst = [0.0 0.0;
            0.0 10.0;
            10.0 10.0]
    iimg_warped, bb = meshwarp(img, src, dst, triangles)
    @test bb == BoundingBox(-1,-1,12,12) 

    src = [-10.0 0.0;
            0.0 10.0;
            10.0 10.0]  
    dst = [0.0 0.0;
            0.0 10.0;
            10.0 10.0]
    img_warped, bb = meshwarp(img, src, dst, triangles)
    @test bb == BoundingBox(-1,-1,12,12) 
end