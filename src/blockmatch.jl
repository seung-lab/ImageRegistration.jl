const NO_MATCH = ([0 0], -1)
const NO_RANGE = (0:0, 0:0)

type Matches
  mesh_indices::Array{Int64, 1}
  src_points::Array
  dst_points::Array
end

"""
`GET_MAX_XC_VECTOR` - Find the cross correlation peak between two images

`[i, j, r_value], correlogram = get_max_xc_vector(A, B)`

* A: 2D array representing the first image
* B: 2D array representing the second image
* i: i displacement from the correlogram center to its peak
* j: j displacement from the correlogram center to its peak
* r_value: the value of the correlogram at its peak
* correlogram: 2D array representing the cross corelation between A & B
"""
function get_max_xc_vector(A, B)
  if std(A) == 0 || std(B) == 0
    return NO_MATCH, 0
  end
  xc = normxcorr2(A, B)
  r_max = maximum(xc)
  if isnan(r_max) 
    return NO_MATCH, 0
  end
  rad = round(Int64, (size(xc, 1) - 1)/ 2)  
  ind = findfirst(r_max .== xc)
  x1 = size(xc, 1)
  x2 = size(xc, 2)
  if ind == 0 
    return NO_MATCH, 0
  end
  (i_max, j_max) = (rem(ind, size(xc, 1)), cld(ind, size(xc, 1)))
  if i_max == 0 
    i_max = size(xc, 1)
  end
  return [i_max-1-rad j_max-1-rad], r_max, xc
end

"""
`BLOCKMATCH` - Find block matches between two images at node locations. 

`matches = get_blockmatches(nodes, src_img, dst_img, src_offset, dst_offset, params)`

* nodes: Nx2 array of nodes for locations in an image where blockmatches will be 
    conducted. Nodes are expected to be in global coordinates (not local to the 
    image's pixel space).
* src_img: 2D array representing the image that will be moving when warped
* dst_img: 2D array representing the image that will stay fixed (aka template)
* src_offset: 2-element array representing the location of src_img's [0,0] pixel
    in the global space.
* dst_offset: 2-element array representing the location of dst_img's [0,0] pixel
    in the global space.
* params: Dict object containing elements defining:
  * block_size: radius from a src node that will be sliced from src_img and used 
      in the cross correlation
  * search_r: additional radius beyond block_size that will be sliced from
      dst_img and used in the cross correlation
  * min_r: the minimum threshold of a cross correlation peak for a blockmatch to
      be accepted.

This method is parallelized. Make sure to start Julia session with additional 
processors:

`julia -p n`

"""
function blockmatch(nodes, src_img, dst_img, src_offset, dst_offset, params)
  n = size(nodes, 1)
  displacements = Array{Tuple{Array{Float64, 2}, Float64}, 1}(n)
  src_ranges = Array{Tuple{UnitRange{Int64}, UnitRange{Int64}}, 1}(n)
  dst_ranges = Array{Tuple{UnitRange{Int64}, UnitRange{Int64}}, 1}(n)

  block_size = params["block_size"]
  search_r = params["search_r"]
  min_r = params["min_r"]
  b_rad = block_size + search_r

  src_sz = size(src_img)
  dst_sz = size(dst_img)

  for i in 1:n
    pt = nodes[i]
    src_ranges[i] = get_range(src_sz, pt, src_offset, block_size)
    dst_ranges[i] = get_range(dst_sz, pt, dst_offset, b_rad)
  end

  num_procs = nprocs()
  k = 1
  nexti() = (i=k; k+=1; i)
  @sync begin
  for p in 1:num_procs
    if p != myid() || num_procs == 1
      @async begin
      while true
        i = nexti()
        if i > n
          break
        end
        if src_ranges[i] == NO_RANGE || dst_ranges[i] == NO_RANGE
          displacements[i] = NO_MATCH
          continue
        end

        src = src_img[src_ranges[i]...]
        dst = dst_img[dst_ranges[i]...]

        max_vect_xc = remotecall_fetch(get_max_xc_vector, p, src,  dst)
        displacements[i] = max_vect_xc[1:2]
      end
      end
    end
  end
  end

  mesh_indices = []
  src_points = []
  dst_points = []
  for i in 1:n
    v = displacements[i]
    if v == NO_MATCH 
      continue
    end
    if v[2] < min_r 
      continue
    end
    push!(mesh_indices, i)
    push!(src_points, nodes[i, :])
    push!(dst_points, nodes[i, :] + v[1])
  end
  if length(dst_points) == 0
    return Void
  end
  return Matches(mesh_indices, vcat(src_points...), vcat(dst_points...))
end

"""
Is a point contained within an image perimeter, less a certain radius?
"""
function is_internal(sz, pt, radius)
  above_min_i = pt[1] > radius
  above_min_j = pt[2] > radius
  below_max_i = pt[1] <= (sz[1]-radius)
  below_max_j = pt[2] <= (sz[2]-radius)
  return above_min_i && above_min_j && below_max_i && below_max_j
end

"""
From global point and offset, calculate range of pixels within radius to slice
"""
function get_range(sz, pt, offset, radius)
  if !is_internal(sz, pt-offset, radius)
    return NO_RANGE
  end
  pt_local = pt-offset
  i = round(Int64, ceil(pt_local[1]))
  j = round(Int64, ceil(pt_local[2]))
  i_range = i-radius:i+radius
  j_range = j-radius:j+radius
  return i_range, j_range
end