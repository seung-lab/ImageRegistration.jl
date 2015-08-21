# Adapted from PiecewiseAffineTransforms.jl
# https://github.com/dfdx/PiecewiseAffineTransforms.jl

include("BoundingBox.jl")

function create_affine(X, Y, U, V)
    F = [X[1] X[2] X[3]; 
        Y[1] Y[2] Y[3]; 
        1 1 1]
    T = [U[1] U[2] U[3]; 
        V[1] V[2] V[3]]
    return T * inv(F) # F*A=T => T*F^-1=A
end

"""
Find the extrema of a mesh, and generate a bounding box

Args:

* nodes: 2xN array of coordinates for a mesh

Returns:

* BoundingBox containing all of the mesh nodes

    BoundingBox(xlow, ylow, height, width) = find_bounds(nodes)
"""
function find_mesh_bb(nodes)
    xlow = floor(Int64,minimum(nodes[:,1]))
    ylow = floor(Int64,minimum(nodes[:,2]))
    xhigh = ceil(Int64,maximum(nodes[:,1]))
    yhigh = ceil(Int64,maximum(nodes[:,2]))
    return BoundingBox(xlow, ylow, xhigh-xlow, yhigh-ylow)
end

function meshwarp{N}(img::Array{Float64, N},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                     trigs::Matrix{Int64}, offset=[0,0], padded=false, interp=true)
    padded=false
    if padded
        wbb = snap_bb(find_mesh_bb(dst))
        low, high = bounds2padding(size(img), minsandmax(wbb)...)
        img = padimage(img, low..., high...)
        src = src .+ low'
        dst = dst .+ low'
        bb = BoundingBox(offset.-low..., size(img)...)
        warped = zeros(eltype(img), size(img))
    else
        bb = snap_bb(find_mesh_bb(dst))
        warped_img = similar(img, bb.h+1, bb.w+1)
        warped_offset = [bb.i, bb.j]
#        dst = dst .- [bb.i, bb.j]'
    end

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
        M = create_affine(U, V, X, Y)
        us, vs = poly2source(U, V)
        
        # for every pixel in target triangle we find corresponding pixel in source
        # and copy its value
        if interp
            for i=1:length(vs)
                u, v = us[i]-1+warped_offset[1], vs[i]-1+warped_offset[2]
                # x, y = M * [u, v, 1]
                x, y = M[1,1]*u + M[1,2]*v + M[1,3], M[2,1]*u + M[2,2]*v + M[2,3]
                x, y = x-offset[1]+1, y-offset[2]+1
                fx, fy = floor(Int64, x), floor(Int64, y)
                wx, wy = x-fx, y-fy

                if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
                    # Expansion of p = [1-wy wy] * img[fy:fy+1, fx:fx+1] * [1-wx; wx]
                    p = ((1-wy)*img[fx,fy] + wy*img[fx,fy+1]) * (1-wx) + ((1-wy)*img[fx+1,fy] + wy*img[fx+1,fy+1]) * wx
                    warped_img[u, v] = p
                end
            end
        else
            for i=1:length(vs)
                u, v = us[i], vs[i]
                # x, y = M * [u, v, 1]
                x = round(Int64, M[1,1]*u + M[1,2]*v + M[1,3])
                y = round(Int64, M[2,1]*u + M[2,2]*v + M[2,3])
                if 1 <= x && x <= size(img, 1) && 1 <= y && y <= size(img, 2)
                    warped_img[u, v] = img[x, y]
                end
            end
        end
        
    end
    return warped_img, [bb.i, bb.j]
end

function poly2source(px, py)
# Run fillpoly2 and findn over the basic bounding box of a given triangle
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    top, bottom = floor(Int64,minimum(py)), ceil(Int64,maximum(py))
    mask = zeros(Bool, bottom-top+1, right-left+1)
    fill = fillpoly2!(mask, px-left+1, py-top+1, true)
    us, vs = findn(fill)
    vs += left-1
    us += top-1
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

function test_meshwarp()
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
    pimg_warped, warped_offset = meshwarp(img, src, dst, triangles, offset, true)
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset, false)
    view(img_warped)
    view(pimg_warped)
    @test warped_offset == [0,0]

    offset = [10, 10]
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    @test warped_offset == [9,9]

    # tform = [cos(pi/6) -sin(pi/6) 0;
    #         sin(pi/6) cos(pi/6) 0;
    #         0 0 1]
    # src = [-16.0 0.0;
    #         27.0 0.0;
    #         6.0 40.0]    
    # dst = warp_pts(tform, src)
    # offset = [0, 0]
    # img_warped, bb = meshwarp(img, src, dst, triangles)
    # @test bb == BoundingBox(0,0,11,11)

    tform = [1 0 0;
            0 1 0;
            10 10 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]   
    dst = warp_pts(tform, src)
    pimg_warped, warped_offset = meshwarp(img, src, dst, triangles, offset, true)
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset, false)
    view(img_warped)
    view(pimg_warped)
    @test bb == BoundingBox(0,0,20,20) 

    tform = [1 0 0;
            0 1 0;
            -4 -4 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]    
    dst = warp_pts(tform, src)
    pimg_warped, warped_offset = meshwarp(img, src, dst, triangles, offset, true)
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset, false)
    view(img_warped)
    view(pimg_warped)
    @test bb == BoundingBox(-4,-5,16,15) 

    src = [0.0 0.0;
            0.0 10.0;
            20.0 20.0]  
    dst = [0.0 0.0;
            0.0 10.0;
            10.0 10.0]
    img_warped, bb = meshwarp(img, src, dst, triangles)
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
