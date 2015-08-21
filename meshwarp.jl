# Adapted from PiecewiseAffineTransforms.jl
# https://github.com/dfdx/PiecewiseAffineTransforms.jl

include("BoundingBox.jl")

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

"""
MESHWARP Apply piecewise affine transform to image using bilinear interpolation

See definitions in IMWARP documentation for further help.

Args:

* img: 2D array, image (todo: extend to Image type)
* src: Nx2 array of mesh nodes that will be deformed _defined in global space_
* dst: Nx2 array of mesh nodes that have been deformed _defined in global space_
* trigs: Nx3 array defining list of triangles - each row contains indices 
    defining which nodes compose a triangle
* offset: 2-element array, position of img[1,1] in 2D space, so also the offset 
    of the src nodes (optional - default is [0,0])
* interp: bool determining whether to use bilinear interpolation or not
    (optional - default is true)

Returns:

* warped_img: with pixel values the same type as original image
    (for Int type, pixel values are rounded), and contained by the bounding
    box of the dst mesh
* warped_offset: 2-element array, position of warped_img[1,1] in 2D space 
"""
function meshwarp{N}(img::Array{Float64, N},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                    trigs::Matrix{Int64}, offset=[0,0], interp=true)
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
        # coordinates of the source triangle
        X = src[tr, 1]
        Y = src[tr, 2]
        # coordinates of the destination triangle
        U = dst[tr, 1]
        V = dst[tr, 2]
        # Create matrix to transform from destination triangle to source image
        src_tri = [X Y ones(3,1)]
        dst_tri = [U V ones(3,1)]
        M = dst_tri \ src_tri # dst_tri * M = src_tri
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
    fillpoly2!(mask, pts_j-left+1, pts_i-top+1, true)
    # Create list of pixel coordinates that are contained by the triangle
    us, vs = findn(mask)
    # Convert that list of pixel coordinates back into global space
    us += top-1
    vs += left-1
    return us, vs
end

"""
Update matrix to value at coordinates contained by the polygon defined by px, py

Args:

* M: 2D array
* px: 1D array of x-components of polygon vertices
* py: 1D array of y-components of polygon vertices
* value: the value to set an array element to if its contained in the polygon

Returns:

* (updated M)
"""
function fillpoly2!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T)
    @assert length(py) == length(px)    
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    # Scan poly from left to right
    for x=left:right     
        ys = Set{Int64}()
        m = length(px)
        for n=1:length(px)
            # Check if x falls between two successive vertices
            # Then add minimum y and maximum y of the polygon at that x            
            if (px[n] <= x && x <= px[m]) || (px[m] <= x && x <= px[n])
                # Vertices are directly on top of each other - add both y's                           
                if px[n] == px[m]
                    push!(ys, ceil(Int64, py[n]))
                    push!(ys, ceil(Int64, py[m]))
                else
                    # Equation of a line between two points, evaluated @ x
                    y = py[n] + (x-px[n]) * (py[m]-py[n])/(px[m]-px[n])
                    push!(ys, ceil(Int64, y))
                end            
            end
            m = n
        end
        # ys is now an array defining ordered point pairs defining slices inside
        # the polygon
        ys = sort([y for y in ys])
        # if odd number of intersection points, add duplicate end point
        if length(ys) % 2 == 1
            push!(ys, ys[end])
        end
        # Place value in matrix at all the y's between min and max for given x
        for n=1:2:length(ys)           
            M[ys[n]:ys[n+1], x] = value
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
    triangles = [1 2 3]

    tform = [1 0 0;
            0 1 0;
            0 0 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]   
    dst = warp_pts(tform, src)
    offset = [0,0]
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    @test warped_offset == [0,0] 

    tform = [1 0 0;
            0 1 0;
            0 0 1]
    src = [0.0 0.0;
            0.0 10.0;
            10.0 0.0]   
    dst = warp_pts(tform, src)
    offset = [5,10]
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    @test_approx_eq_eps img_warped zeros(11,11) 1e-10
    @test warped_offset == [0,0] 

    tform = [1 0 0;
            0 1 0;
            0 0 1]
    src = [0.0 0.0;
            0.0 20.0;
            20.0 0.0]   
    dst = warp_pts(tform, src)
    offset = [5,10]
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    @test img_warped[6,11] == 1.0 
    @test warped_offset == [0,0] 

    src = [2.0 2.0;
            10.0 2.0;
            6.0 10.0]
    dst = [4.0 2.0;
            8.0 2.0;
            6.0 10.0]
    offset = [0, 0]
    img_warped, warped_offset = meshwarp(img, src, dst, triangles, offset)
    @test warped_offset == [4,2]
end
