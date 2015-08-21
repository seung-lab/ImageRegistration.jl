# Adapted from PiecewiseAffineTransforms.jl
# https://github.com/dfdx/PiecewiseAffineTransforms.jl

include("BoundingBox.jl")

function create_affine(U, V, X, Y)
    pts = [X[1] X[2] X[3]; 
            Y[1] Y[2] Y[3]; 
            1 1 1]
    tpts = [U[1] U[2] U[3]; 
            V[1] V[2] V[3];
            1 1 1]
    # tpts*A = pts => tpts^-1 * pts
    return inv(tpts) * pts
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
                     trigs::Matrix{Int64}, offset=[0,0], interp=true)
    # wbb = snap_bb(find_mesh_bb(dst))
    # low, high = bounds2padding(size(img), minsandmax(wbb)...)
    # img = padimage(img, low..., high...)
    # src = src .+ low'
    # dst = dst .+ low'
    # bb = BoundingBox(offset.-low..., size(img)...)
    # warped = zeros(eltype(img), size(img))

    bb = snap_bb(find_mesh_bb(dst))
    warped_img = similar(img, bb.h+1, bb.w+1)
    warped_offset = [bb.i, bb.j]

    if interp
        println("w interpolation")
    else
        println("w/o interpolation")
    end
    for t=1:size(trigs, 1)    
        tr = squeeze(trigs[t, :], 1)
        X = src[tr, 1]
        Y = src[tr, 2]
        
        U = dst[tr, 1]
        V = dst[tr, 2]
               
        # Create matrix to transform from warped image to original image
        M = create_affine(U, V, X, Y) # tform U,V to X,Y
        # Create list of coordinates in warped image that represent this
        # triangle (coordinates in the global space).
        us, vs = poly2source(U, V)
        
        # For every pixel in target triangle we find corresponding pixel in 
        # source and copy its value
        if interp
            for n=1:length(vs)
                # Convert warped coordinate to pixel space
                i, j = us[n]-warped_offset[1]+1, vs[n]-warped_offset[2]+1
                # Use warped coordinate in global space for transform
                u, v = us[n], vs[n]
                # x, y = M * [u, v, 1]
                x, y = M[1,1]*u + M[2,1]*v + M[3,1], M[1,2]*u + M[2,2]*v + M[3,2]
                # Convert original image coordinate to pixel space
                x, y = x-offset[1]+1, y-offset[2]+1
                fx, fy = floor(Int64, x), floor(Int64, y)
                wx, wy = x-fx, y-fy
                if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
                    # Expansion of p = [1-wy wy] * img[fy:fy+1, fx:fx+1] * [1-wx; wx]
                    p = ((1-wy)*img[fx,fy] + wy*img[fx,fy+1]) * (1-wx) + ((1-wy)*img[fx+1,fy] + wy*img[fx+1,fy+1]) * wx
                    warped_img[i, j] = p
                end
            end
        else
            for n=1:length(vs)
                i, j = us[n]-warped_offset[1]+1, vs[n]-warped_offset[2]+1
                u, v = us[n], vs[n]
                # x, y = M * [u, v, 1]
                x = M[1,1]*u + M[2,1]*v + M[3,1]
                y = M[1,2]*u + M[2,2]*v + M[3,2]
                x, y = round(Int64, x-offset[1]+1), round(Int64, y-offset[2]+1)
                if 1 <= x && x <= size(img, 1) && 1 <= y && y <= size(img, 2)
                    warped_img[i, j] = img[x, y]
                end
            end
        end
        
    end
    return warped_img, [bb.i, bb.j]
end

"""
Run fillpoly2 and findn over the basic bounding box of a given triangle
"""
function poly2source(pts_i, pts_j)
    # Find bb of vertices (vertices in global space)
    top, bottom = floor(Int64,minimum(pts_i)), ceil(Int64,maximum(pts_i))
    left, right = floor(Int64,minimum(pts_j)), ceil(Int64,maximum(pts_j))
    # Create image based on number of pixels in bb that will identify triangle
    mask = zeros(Bool, bottom-top+1, right-left+1)
    # Convert vertices into pixel space and fill the mask to identify triangle
    fillpoly2!(mask, pts_i-top+1, pts_j-left+1, true)
    # Create list of pixel coordinates that are contained by the triangle
    us, vs = findn(mask)
    # Convert that list of pixel coordinates back into global space
    us += top-1
    vs += left-1
    return us, vs
end

function fillpoly2!{T,P<:Number}(M::Matrix{T}, pts_i::Vector{P}, pts_j::Vector{P}, value::T)
    @assert length(pts_i) == length(pts_j)    
    left, right = floor(Int64,minimum(pts_j)), ceil(Int64,maximum(pts_j))
    for x=left:right     
        ys = Set{Int64}()
        j = length(pts_j)
        for i=1:length(pts_j)            
            if (pts_j[i] <= x && x <= pts_j[j]) || (pts_j[j] <= x && x <= pts_j[i])
                # special case: adding the whole cut to ys                            
                if pts_j[i] == pts_j[j]
                    push!(ys, ceil(Int64, pts_i[i]))
                    push!(ys, ceil(Int64, pts_i[j]))
                else
                    y = pts_i[i] + (x - pts_j[i]) / (pts_j[j] - pts_j[i]) * (pts_i[j] - pts_i[i])
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
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    view(img_warped)
    @test warped_offset == [0,0]

    offset = [10, 10]
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    @test warped_offset == [9,9]

    tform = [cos(pi/4) -sin(pi/4) 0;
            sin(pi/4) cos(pi/4) 0;
            0 0 1]
    src = [1.0 1.0;
            11.0 1.0;
            1.0 11.0]    
    dst = warp_pts(tform, src)
    offset = [0, 0]
    img_warped, bb = meshwarp(img, src, dst, triangles)
    view(img_warped)
    # @test bb == BoundingBox(0,0,11,11)

    tform = [1 0 0;
            0 1 0;
            0 0 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]   
    dst = warp_pts(tform, src)
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    view(img_warped)
    @test bb == BoundingBox(0,0,20,20) 

    tform = [1 0 0;
            0 1 0;
            -4 -4 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]    
    dst = warp_pts(tform, src)
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    view(img_warped)
    @test bb == BoundingBox(-4,-5,16,15) 

    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]  
    src = [0.0 0.0;
            0.0 9.5;
            10.0 0.0] 
    img_warped, bb = meshwarp(img, src, dst, triangles)
    view(img_warped)
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
