# Initially adapted from PiecewiseAffineTransforms.jl
# https://github.com/dfdx/PiecewiseAffineTransforms.jl

global MESHWARP_POLY2SOURCE_MASK = zeros(Bool, 0, 0);

"""
`MESHWARP` - Apply piecewise affine transform to image using bilinear interpolation

    warped_img, [bb.i, bb.j] = meshwarp(img, src, dst, trigs, offset, interp)

* `img`: 2D array, image (todo: extend to Image type)
* `src`: 2xN array, positions of mesh nodes in `img` (global space)
* `dst`: 2xN array, positions of mesh nodes in `warped_img` (global space)
* `trigs`: 3xN array, each column contains indices of nodes of a triangle
* `offset`: 2-element array, position of img[1,1] in global space (optional - default is [0,0])
* `interp`: bool determining whether to use bilinear interpolation or not
    (optional - default is true)
* `warped_img`: with pixel values the same type as original image (for Int type, pixel values are rounded)
* `warped_offset`: 2-element array, position of warped_img[1,1] in global space 

The bounding box of the warped image is defined as the smallest
integer-valued rectangle that contains the `dst` mesh.

This means that `warped_offset` is constrained to be integer-valued,
though `offset` is allowed to have floating point values.  The integer
constraint removes the need for further interpolation in any
subsequent fusing of multiple warped tiles.

See definitions in IMWARP documentation for further help.

""" 
function meshwarp{T}(img::SharedArray{T},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                    trigs::Matrix{Int64}, offset=[0,0], interp=true)

  bb = snap_bb(find_mesh_bb(dst))
  warped_img = SharedArray{T}(bb.h+1, bb.w+1)
  warped_offset = [bb.i, bb.j];

  Us = Array{Array{Float64, 1}}(size(trigs, 2));
  Vs = Array{Array{Float64, 1}}(size(trigs, 2));
  Ms = Array{Array{Float64, 2}}(size(trigs, 2));

  @everywhere gc();

  src_tri = ones(Float64, 3, 3)
  dst_tri = ones(Float64, 3, 3)

  @fastmath @inbounds for t in 1:size(trigs, 2);
    #println("$t / $(size(trigs, 1))")
    #tr = squeeze(trigs[t, :], 1) = necessary for 0.4.6
    tr = view(trigs, :, t)
    #= slow version
    # coordinates of the source triangle
    X = src[tr, 1]
    Y = src[tr, 2]
    # coordinates of the destination triangle
    U = dst[tr, 1]
    V = dst[tr, 2]
    # Create matrix to transform from destination triangle to source image
    src_tri = [X Y ones(3,1)]
    dst_tri = [U V ones(3,1)]
    =#
    U = dst[1, tr]
    V = dst[2, tr]
    # Create matrix to transform from destination triangle to source image
    src_tri[:, 1] = view(src, 1, tr)
    src_tri[:, 2] = view(src, 2, tr)

    dst_tri[:, 1] = U
    dst_tri[:, 2] = V

    M = dst_tri \ src_tri # dst_tri * M = src_tri
    # Create list of coordinates in warped image that represent this
    # triangle (coordinates in the global space).
    Us[t] = U;
    Vs[t] = V;
    Ms[t] = M;
  end

  function proc_range(idx, arr::Array)
    worker_procs = setdiff(procs(), myid());
    nchunks = length(worker_procs);
    if nchunks == 0 return 1:length(arr); end
    if idx == myid() return 1:0; end
    splits = [round(Int64, s) for s in linspace(0, length(arr), nchunks + 1)];
    return splits[findfirst(worker_procs, idx)]+1:splits[findfirst(worker_procs, idx) + 1]
  end

  #=@sync @fastmath @inbounds for t in 1:size(trigs, 1);
      @async remotecall_wait(procs()[rem(t, nprocs()-1)+2], calculate_pixels_in_trig!, Us[t], Vs[t], Ms[t], img, offset, warped_img, warped_offset)
  end=#

  @sync for p in procs()
          t = proc_range(p, Us);
          @async @fastmath @inbounds remotecall_wait(calculate_pixels_in_trig_chunk!, p, Us[t], Vs[t], Ms[t], img, offset, warped_img, warped_offset, interp ? computepixel_interp : computepixel)
  end

  return warped_img, [bb.i, bb.j]
end

function calculate_pixels_in_trig_chunk!(Us, Vs, Ms, img, offset, warped_img, warped_offset, interp)
 @simd for i in 1:length(Us)
	calculate_pixels_in_trig!(Us[i], Vs[i], Ms[i], img, offset, warped_img, warped_offset, interp)
  end
  global MESHWARP_POLY2SOURCE_MASK = zeros(Bool, 0, 0);
end

function calculate_pixels_in_trig!(U, V, M, img, offset, warped_img, warped_offset, computepixel_func)
  poly::Array{Bool,2}, i_oset::Int64, j_oset::Int64 = poly2source_new!(U, V)
    for j_poly in 1:size(poly,2)
    	for i_poly in 1:size(poly,1)
    @fastmath @inbounds begin
	#check if poly is true
      	if poly[i_poly, j_poly] == false continue end
        # Use warped coordinate in global space for transform - this changes back from the test polygon to global space
	u, v = i_poly + i_oset, j_poly + j_oset
        # Convert warped coordinate to pixel space
        i, j = u-warped_offset[1]+1, v-warped_offset[2]+1
        # x, y = M * [u, v, 1]
        @fastmath @inbounds x, y = M[1,1]*u + M[2,1]*v + M[3,1], M[1,2]*u + M[2,2]*v + M[3,2]
        # Convert original image coordinate to pixel space
        x, y = x-offset[1]+1, y-offset[2]+1
        fx, fy = floor(Int64, x), floor(Int64, y)
        wx, wy = x-fx, y-fy
        if 1 <= fx <= size(img, 1)-1 && 1 <= fy <= size(img, 2)-1
	computepixel_func(warped_img, img, i, j, wx, wy, fx, fy)
        end
      end
      end
    end
#=
   us, vs = poly2source!(U, V)
#    computepixel_func = interp ? computepixel_interp : computepixel
    @simd for ind in 1:length(us)
        # Convert warped coordinate to pixel space
    @fastmath @inbounds begin
        i, j = us[ind]-warped_offset[1]+1, vs[ind]-warped_offset[2]+1
        # Use warped coordinate in global space for transform
        u, v = us[ind], vs[ind]
        # x, y = M * [u, v, 1]
        @fastmath @inbounds x, y = M[1,1]*u + M[2,1]*v + M[3,1], M[1,2]*u + M[2,2]*v + M[3,2]
        # Convert original image coordinate to pixel space
        x, y = x-offset[1]+1, y-offset[2]+1
        fx, fy = floor(Int64, x), floor(Int64, y)
        wx, wy = x-fx, y-fy
        if 1 <= fx <= size(img, 1)-1 && 1 <= fy <= size(img, 2)-1
	computepixel_func(warped_img, img, i, j, wx, wy, fx, fy)
        end
      end
    end
    =#
end

    # Takes weights wx, wy that denotes how close the value is to [fx+1, fy+1] from [fx, fy] and writes the weighted average to [i, j]
    @inline function computepixel_interp(warped_img, img, i, j, wx, wy, fx, fy)
          # Expansion of p = [1-wy wy] * img[fy:fy+1, fx:fx+1] * [1-wx; wx]
          @fastmath @inbounds p1 = ((1-wy)*img[fx,fy] + wy*img[fx,fy+1])
	  @fastmath @inbounds p2 = ((1-wy)*img[fx+1,fy] + wy*img[fx+1,fy+1])
	  @fastmath p1 = p1 * (1-wx);
	  @fastmath p2 = p2 * (wx);
          writepixel(warped_img,i,j,p1+p2)
    end

    function computepixel(warped_img, img, i, j, wx, wy, fx, fy)
      	  @fastmath y = wy < 0.5 ? fy : fy + 1
      	  @fastmath x = wx < 0.5 ? fx : fx + 1
	  @inbounds v = img[x,y]
          writepixel(warped_img,i,j,v)
    end


"""
`POLY2SOURCE` - Run fillpoly and findn over the basic bounding box of a given triangle
""" 
function poly2source(pts_i, pts_j)
  # Find bb of vertices (vertices in global space)
  top, bottom = floor(Int64,minimum(pts_i)), ceil(Int64,maximum(pts_i))
  left, right = floor(Int64,minimum(pts_j)), ceil(Int64,maximum(pts_j))
  # Create image based on number of pixels in bb that will identify triangle
  mask = zeros(Bool, bottom-top+1, right-left+1)
  # Convert vertices into pixel space and fill the mask to identify triangle
  fillpoly!(mask, pts_j-left+1, pts_i-top+1, true; convex = false)
  # Create list of pixel coordinates that are contained by the triangle
  us, vs = findn(mask)
  # Convert that list of pixel coordinates back into global space
  us += top-1
  vs += left-1
  return us, vs
end

@inline function poly2source_new!(pts_i, pts_j)
  # Find bb of vertices (vertices in global space)
  top, bottom = floor(Int64,minimum(pts_i)), ceil(Int64,maximum(pts_i))
  left, right = floor(Int64,minimum(pts_j)), ceil(Int64,maximum(pts_j))
  # Create image based on number of pixels in bb that will identify triangle
  if size(MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2})[1] < bottom - top + 1 || size(MESHWARP_POLY2SOURCE_MASK)[2] < right - left + 1
    global MESHWARP_POLY2SOURCE_MASK = zeros(Bool, bottom-top+1, right-left+1)
  end
  (MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2})[:] = false;

  #mask = zeros(Bool, bottom-top+1, right-left+1)
  # Convert vertices into pixel space and fill the mask to identify triangle
  for i in 1:length(pts_i)
    @fastmath @inbounds pts_i[i] = pts_i[i] + 1 - top;
  end
  for i in 1:length(pts_j)
    @fastmath @inbounds pts_j[i] = pts_j[i] + 1 - left;
  end
  #fillpoly!(mask, pts_j-left+1, pts_i-top+1, true)
   fillpoly_convex!(MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2}, pts_j, pts_i, true)
   # return the offsets that have to be added
   return MESHWARP_POLY2SOURCE_MASK, top - 1, left - 1
   #=
  # Create list of pixel coordinates that are contained by the triangle
  us, vs = findn(MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2})
  # Convert that list of pixel coordinates back into global space
  for i in 1:length(us)
    @fastmath @inbounds us[i] = us[i] - 1 + top;
  end
  for i in 1:length(vs)
    @fastmath @inbounds vs[i] = vs[i] - 1 + left;
  end
#  us += top-1
#  vs += left-1
  return us, vs
  =#
end


function poly2source!(pts_i, pts_j)
  # Find bb of vertices (vertices in global space)
  top, bottom = floor(Int64,minimum(pts_i)), ceil(Int64,maximum(pts_i))
  left, right = floor(Int64,minimum(pts_j)), ceil(Int64,maximum(pts_j))
  # Create image based on number of pixels in bb that will identify triangle
  if size(MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2})[1] < bottom - top + 1 || size(MESHWARP_POLY2SOURCE_MASK)[2] < right - left + 1
    global MESHWARP_POLY2SOURCE_MASK = zeros(Bool, bottom-top+1, right-left+1)
  end
  (MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2})[:] = false;
  #mask = zeros(Bool, bottom-top+1, right-left+1)
  # Convert vertices into pixel space and fill the mask to identify triangle
  for i in 1:length(pts_i)
    @fastmath @inbounds pts_i[i] = pts_i[i] + 1 - top;
  end
  for i in 1:length(pts_j)
    @fastmath @inbounds pts_j[i] = pts_j[i] + 1 - left;
  end
  #fillpoly!(mask, pts_j-left+1, pts_i-top+1, true)
   fillpoly_convex!(MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2}, pts_j, pts_i, true)
  # Create list of pixel coordinates that are contained by the triangle
  us, vs = findn(MESHWARP_POLY2SOURCE_MASK::Array{Bool, 2})
  # Convert that list of pixel coordinates back into global space
  for i in 1:length(us)
    @fastmath @inbounds us[i] = us[i] - 1 + top;
  end
  for i in 1:length(vs)
    @fastmath @inbounds vs[i] = vs[i] - 1 + left;
  end
#  us += top-1
#  vs += left-1
  return us, vs
end

"""
`FILLPOLY!` - Fill pixels contained inside a polygon

    M = fillpoly!(M, px, py, value; reverse)

* `M`: 2D array, image
* `px`: 1D array, x-components of polygon vertices
* `py`: 1D array, y-components of polygon vertices
* `value`: fill value
* `convex`: Boolean kwarg, false by default; if true then assumes the polygon is convex and well-formed without duplicates
* `reverse`: Boolean kwarg, false by default; if true then fills in everything outside the polygon

Features:

* Intended to work for nonconvex polygons as well as convex polygons
* The vertices should be listed sequentially in px and py. Clockwise/counterclockwise ordering doesn't matter.
* Duplicate points in sequence are handled by removing until only one remains.
* Polygons need not be closed (i.e. last point needs not be the same as the first point) - if closed, they will be handled as duplicates initially.
* The code treats the "corner case," where the vertex of the polygon intersects the grid line at one point.
* The code treats the "edge case," where an edge of the polygon lies exactly on a grid line.
* The boundary of the polygon is filled whether or not reverse kwarg is set - i.e. the results with `reverse` set to true and false are not perfect complements due to the boundary.

Bugs:

* The original code (https://github.com/dfdx/PiecewiseAffineTransforms.jl/blob/master/src/polyline.jl) used implicit conversion to Int64 (presumably rounding).  Tommy/Shang replaced this by ceil.  This might produce inconsistent results, with the same pixel belonging to more than one triangle.
""" 
function fillpoly!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T; convex::Bool=false, reverse::Bool=false)
  convex ? (return fillpoly_convex!(M, px, py, value; reverse = reverse)) : (return fillpoly_nonconvex!(M, px, py, value; reverse = reverse))
end

function fillpoly!{T,P<:Number}(M::SharedArray{T}, px::Vector{P}, py::Vector{P}, value::T; convex::Bool=false, reverse::Bool=false)
  M = convert(Array, M)
  fillpoly!(M, px, py, value; convex=convex, reverse=reverse)
  M = convert(SharedArray, M)
end

function fillpoly_convex!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T; reverse::Bool=false)
  left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
  if reverse xrange = 1:size(M, 2)
  else xrange = left:right
  end
  l = length(px)
  # Scan poly from left to right
  for x=xrange     # loop over grid lines
    top = 0
    bot = 0
    m = length(px)
    xdir_last = sign(px[m] - px[m-1]) # direction of the line segment
    for n=1:length(px)  # loop over edges (m,1),(1,2),(2,3),...,(m-1,m)
      xdir_cur = sign(px[n] - px[m])
      # grid line intersects edge in one of two ways
      if (px[n] <= x <= px[m]) || (px[m] <= x <= px[n])
	if px[n] == px[m]  # intersection is entire edge, do nothing
	else# intersection is point
	  # do not add duplicate points (endpoint of the last segment), unless the direction is reversing in x or the last segment is vertical
	  if px[m] != x || xdir_last + xdir_cur == 0 || xdir_last == 0
            y = py[n] + (x-px[n]) * (py[m]-py[n])/(px[m]-px[n])
	    # deal with rounding error
	    if ceil(Int64, y) > size(M, 1)
	      val = floor(Int64, y)
	    else
              val = ceil(Int64, y)
	    end
	    if top == 0
	      top = val 
	    elseif top < val
	      bot = val
	    else
	      bot = top
	      top = val
	      break
	    end
       	  end            
      	end
      end
      xdir_last = xdir_cur
      m = n
    end

    if top == 0 || bot == 0 continue end

    if reverse
      @simd for y in 1:top
      @inbounds M[y, x] = value
      end
      @simd for y in bot:size(M,1)
      @inbounds M[y, x] = value
      end
    else
    # Place value in matrix at all the y's between min and max for given x
      @simd for y in top:bot
      @inbounds M[y, x] = value
      end
    end
  end
  return M
end

function fillpoly_nonconvex!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T; reverse::Bool=false)
  @assert length(py) == length(px)    
  left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
  if reverse xrange = 1:size(M, 2)
  else xrange = left:right
  end
  l = length(px)
  #force no duplicates, including at the end
    for i in 1:length(px)
	if px[i] == px[l] && py[i] == py[l]
	  px[l] = 0
	  py[l] = 0
	end
	  l = i
    end
    px = px[px .!= 0]
    py = py[py .!= 0]
  # Scan poly from left to right
  for x=xrange     # loop over grid lines
    ys = Array{Int64, 1}()
    tops = Set{Int64}() # tops of vertical edges
    signs = Array{Int64, 1}()
    ys_clean = Array{Int64, 1}()
    m = length(px)
    xdir_last = sign(px[m] - px[m-1]) # direction of the line segment
    for n=1:length(px)  # loop over edges (m,1),(1,2),(2,3),...,(m-1,m)
      xdir_cur = sign(px[n] - px[m])
      # grid line intersects edge in one of two ways
      if (px[n] <= x <= px[m]) || (px[m] <= x <= px[n])
	if px[n] == px[m]  # intersection is entire edge
#	  if xdir_last + xdir_cur == 0 # if two vertical edges in a row remove the last one
#	    pop!(ys)
#	  end
	  push!(tops, ceil(Int64, minimum((py[n], py[m]))))
	  #push!(signs, xdir_cur)
          #push!(ys, ceil(Int64, py[n]))
	else# intersection is point
	  # do not add duplicate points (endpoint of the last segment), unless the direction is reversing in x or the last segment is vertical
	  if px[m] != x || xdir_last + xdir_cur == 0 || xdir_last == 0
            y = py[n] + (x-px[n]) * (py[m]-py[n])/(px[m]-px[n])
	    push!(signs, xdir_cur)
	    # deal with rounding error
	    if ceil(Int64, y) > size(M, 1)
	      push!(ys, floor(Int64, y))
	    else
              push!(ys, ceil(Int64, y))
	    end
       	  end            
      	end
      end
      xdir_last = xdir_cur
      m = n
    end

    if length(ys) != 2
    # generically, two intersections for a convex polygon
    perm = sortperm(ys)  # sort the intersection points
    ys = ys[perm]  # sort the intersection points
    signs = signs[perm]  # sort the intersection points' signs
    i = 1
    while (i < length(ys))
	y_entry = ys[i]
	y_exit = ys[i]
	entry_sign = signs[i]
	while (i < length(ys))
	  i += 1
	  y_exit = ys[i]
	  exit_sign = signs[i]
	  if in(y_exit, tops) continue;
	  elseif exit_sign != entry_sign break; end
	end
	push!(ys_clean, y_entry)
	push!(ys_clean, y_exit)
	i += 1
    end
  else
    #if ys length is two then we can simply sort instead of going through the whole algorithm
    ys_clean = sort(ys)
    end
    
    if reverse
      unshift!(ys_clean, 1)
      push!(ys_clean, size(M,1))
    end
    # Place value in matrix at all the y's between min and max for given x
    for n=1:2:length(ys_clean)           
      @simd for y in ys_clean[n]:ys_clean[n+1]
      @inbounds M[y, x] = value
      end
    end
  end
  return M
end


function meshwarp{T}(img::Array{T},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                    trigs::Matrix{Int64}, offset=[0,0], interp=true)
  img = convert(SharedArray, img)
  meshwarp(img, src, dst, trigs, offset, interp)
end
