import Base: +, -, ==

immutable BoundingBox{T}
  i::T
  j::T
  h::T    # height
  w::T    # width
end

eps = 1e-6

BoundingBox() = BoundingBox(0,0,0,0)
BoundingBox(a, b, c, d) = BoundingBox(promote(a,b,c,d)...)

"""
Add bounding boxes to find BB of their union
"""
function +(bbA::BoundingBox, bbB::BoundingBox)
  i = min(bbA.i, bbB.i)
  j = min(bbA.j, bbB.j)
  h = max(bbA.h+bbA.i, bbB.h+bbB.i) - i
  w = max(bbA.w+bbA.j, bbB.w+bbB.j) - j
  return BoundingBox(i,j,h,w)
end

"""
Subtract bounding boxes to find BB of their intersection
"""
function -(bbA::BoundingBox, bbB::BoundingBox)
  i = max(bbA.i, bbB.i)
  j = max(bbA.j, bbB.j)
  h = min(bbA.h+bbA.i, bbB.h+bbB.i)-i
  w = min(bbA.w+bbA.j, bbB.w+bbB.j)-j
  if h < 0 || w < 0
    bb = BoundingBox(NaN, NaN, NaN, NaN)
  else
    bb = BoundingBox(i,j,h,w)
  end
  return bb
end

"""
Boolean if bounding boxes intersect
"""
function intersects(bbA::BoundingBox, bbB::BoundingBox)
  bb = bbA - bbB
  return !isnan(bb.i)
end

"""
Test if BoundingBoxes have the same origin and dimensions
"""
function ==(bbA::BoundingBox, bbB::BoundingBox)
  return bbA.i == bbB.i && bbA.j == bbB.j && bbA.w == bbB.w && bbA.h == bbB.h
end

"""
Convert bounding box object to a polygon point list (counter-clockwise)
"""
function bb_to_pts(r)
  return [r.i r.j;
          r.i+r.h r.j;
          r.i+r.h r.j+r.w;
          r.i r.j+r.w;
          r.i r.j];
end

"""
Return xmin, ymin, xmax, ymax of a bounding box (opposing corner definition)
"""
function minsandmax(bb)
  return bb.i, bb.j, bb.i+bb.h, bb.j+bb.w
end

"""
`FIND_MESH_BB` - Find bounding box around mesh

    BoundingBox(ilow, jlow, height, width) = find_mesh_bb(nodes)

* `nodes`: 2xN matrix of mesh nodes
* `BoundingBox`: smallest integer-valued rectangle containing all mesh nodes

""" 
function find_mesh_bb(nodes)
    ilow = floor(Int64,minimum(nodes[:,1]))
    jlow = floor(Int64,minimum(nodes[:,2]))
    ihigh = ceil(Int64,maximum(nodes[:,1]))
    jhigh = ceil(Int64,maximum(nodes[:,2]))
    return BoundingBox(ilow, jlow, ihigh-ilow, jhigh-jlow)
end

"""
Snap bounding box to integer values (smallest rectangle containing original)
Returns BoundingBox of integers
"""
function snap_bb(bb)
  r = bb_to_pts(bb)
  i = floor(Int,bb.i)
  j = floor(Int,bb.j)
  h = ceil(Int,r[2,1]) - i
  w = ceil(Int,r[3,2]) - j
  return BoundingBox{Int}(i, j, h, w)
end

"""
Apply affine transform to points in a bounding box & find new bounding box
"""
function tform_bb(bb, tform)
  tform_pts = [bb_to_pts(bb) ones(size(bb_to_pts(bb),1),1)] * tform
  i = minimum(tform_pts[:,1])
  j = minimum(tform_pts[:,2])
  h = maximum(tform_pts[:,1])-i
  w = maximum(tform_pts[:,2])-j
  return BoundingBox(i, j, h, w)
end

"""
Convert Tuple for image size into a BoundingBox at (1,1)
"""
function sz_to_bb(sz)
  return BoundingBox(0, 0, sz[1], sz[2])
end

"""
Convert a BoundingBox to its snapped sizes
"""
function bb_to_sz(bb)
  bb = snap_bb(bb);
  return bb.h, bb.w
end

"""
Get bounding box offset
"""
function get_offset(bb::BoundingBox)
  return [bb.i, bb.j]
end

"""
Get height & width of bounding box
"""
function get_size(bb::BoundingBox)
  return [bb.h, bb.w]
end

"""
Get four-tuple of bounding box upper-left & lower-right in i,j coordinates
"""
function get_bounds(bb::BoundingBox)
  return (bb.i, bb.j, bb.i+bb.h, bb.j+bb.w)
end

"""
Get four-tuple of bounding box upper-left & dimensions in x,y coordinates

Convention used by Cairo package
"""
function get_rect(bb::BoundingBox)
  return (bb.j, bb.i, bb.w, bb.h)
end

"""
Return the product of the width & height of bounding box
"""
function get_area(bb::BoundingBox)
  return bb.w*bb.h
end

"""
Convert bounding box to tuple of ranges for easy array slicing
"""
function bb_to_slice(bb::BoundingBox{Int64})
  return (bb.i+1):(bb.i+bb.h), (bb.j+1):(bb.j+bb.w)
end

function bb_to_slice(bb::BoundingBox{Float64})
  return round(Int64, bb.i+1):round(Int64, bb.i+bb.h), 
              round(Int64, bb.j+1):round(Int64, bb.j+bb.w)
end

"""
Convert tuple of ranges to bounding box
"""
function slice_to_bb(slice)
  return BoundingBox(slice[1][1], slice[2][1], 
                      slice[1][end]-slice[1][1], slice[2][end]-slice[2][1])
end

"""
Shift bounding box by 2-element array
"""
function translate_bb(bb::BoundingBox, offset)
  return BoundingBox(bb.i + offset[1], bb.j + offset[2], bb.h, bb.w)
end

"""
Transform bounding box by scaling matrix and snap for nearest integer
"""
function scale_bb(bb::BoundingBox{Int64}, scale)
  tform = make_scale_matrix(scale)
  tbb = tform_bb(bb, tform)
  return snap_bb(tbb)
end

"""
Transform bounding box by scaling matrix
"""
function scale_bb(bb::BoundingBox{Float64}, scale)
  tform = make_scale_matrix(scale)
  return tform_bb(bb, tform)
end

"""
Returns true if point, pt, is contained by the polygon defined by points in poly

Pulled from:
https://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

Based on Jordan curve theorem.

Draw a semi-infinite ray from the point (x infinite, y fixed) & count how many
edges that ray crosses (or actually, how many edges it is not parallel to). If
the ray crosses an odd number of edges when its slope is larger than the slope
of the edge, its contained. Only consider edges for which the point falls within
its y-range.

The semi-infinite ray method is inconsistent on boundaries, so include a test
for the border at the very end.
"""
function pt_in_poly(pt, poly)
  is_contained = false
  N = size(poly,1)
  # only count the distinct vertices
  if poly[1,:] == poly[N,:]
    N = N-1
  end
  j = N
  for i = 1:N
    i_vert = poly[i,:][:]
    j_vert = poly[j,:][:]
    # is point is on same side of edges endpoints in y?
    if (pt[2] > i_vert[2]) != (pt[2] > j_vert[2])
      # is the line from pt->i_vert parallel to j_vert->i_vert?
      if pt[1] < (j_vert[1]-i_vert[1])*(pt[2]-i_vert[2])/(j_vert[2]-i_vert[2]) + i_vert[1]
        is_contained = !is_contained
      end
    end
    j = i
  end

  if !is_contained
    j = N
    for i = 1:N
      i_vert = poly[i,:][:]
      j_vert = poly[j,:][:]
      on_line = pt_on_line_segment(pt, (i_vert, j_vert))
      if on_line
        return true
      end
    end
  end

  return is_contained
end

"""
If point exists along the line segment, the triangle it defines should be flat.
"""
function pt_on_line_segment(pt, line)
  epsilon = 1e-4
  a, b = line
  return -epsilon < norm(a-b) - (norm(a-pt) + norm(b-pt)) < epsilon
end

"""
Returns true if polyA contains polyB
"""
function poly_contains_poly(polyA, polyB)
  for i in 1:size(polyB, 1)
    pt = polyB[i,:][:]
    if pt_in_poly(pt, polyA)
      return true
    end
  end
  return false
end

"""
Returns true if polyA & polyB intersect
"""
function poly_intersects(polyA, polyB)
  return poly_contains_poly(polyA, polyB) || poly_contains_poly(polyB, polyA)
end