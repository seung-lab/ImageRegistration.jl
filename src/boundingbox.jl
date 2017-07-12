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
  if h <= 0 || w <= 0
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
    ilow = floor(Int64,minimum(nodes[1,:]))
    jlow = floor(Int64,minimum(nodes[2,:]))
    ihigh = ceil(Int64,maximum(nodes[1,:]))
    jhigh = ceil(Int64,maximum(nodes[2,:]))
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
  return (bb.i):(bb.i+bb.h-1), (bb.j):(bb.j+bb.w-1)
end

function bb_to_slice(bb::BoundingBox{Float64})
  return round(Int64, bb.i):round(Int64, bb.i+bb.h-1), 
              round(Int64, bb.j):round(Int64, bb.j+bb.w-1)
end

"""
Convert tuple of ranges to bounding box
"""
function slice_to_bb(slice)
  return BoundingBox(slice[1][1], slice[2][1], 
                      slice[1][end]-slice[1][1]+1, slice[2][end]-slice[2][1]+1)
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
