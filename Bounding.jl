module Bounding

using Base.Test
import Base: +, ==

type BoundingBox{T<:Real}
  i::T
  j::T
  h::T    # height
  w::T    # width
end

BoundingBox() = BoundingBox(0,0,0,0)

"""
Add bounding boxes by finding BoundingBox that encompasses both
"""
function +(bbA::BoundingBox, bbB::BoundingBox)
  i = min(bbA.i, bbB.i)
  j = min(bbA.j, bbB.j)
  h = max(bbA.h+bbA.i, bbB.h+bbB.i) - i
  w = max(bbA.w+bbA.j, bbB.w+bbB.j) - j
  return BoundingBox(i,j,h,w)
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
function bb2pts(r)
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

    BoundingBox(xlow, ylow, height, width) = find_mesh_bb(nodes)

* `nodes`: 2xN matrix of mesh nodes
* `BoundingBox`: smallest integer-valued rectangle containing all mesh nodes

""" 
function find_mesh_bb(nodes)
    xlow = floor(Int64,minimum(nodes[:,1]))
    ylow = floor(Int64,minimum(nodes[:,2]))
    xhigh = ceil(Int64,maximum(nodes[:,1]))
    yhigh = ceil(Int64,maximum(nodes[:,2]))
    return BoundingBox(xlow, ylow, xhigh-xlow, yhigh-ylow)
end

"""
Snap bounding box to integer values (smallest rectangle containing original)
Returns BoundingBox of integers
"""
function snap_bb(bb)
  r = bb2pts(bb)
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
  tform_pts = [bb2pts(bb) ones(size(bb2pts(bb),1),1)] * tform
  i = minimum(tform_pts[:,1])
  j = minimum(tform_pts[:,2])
  h = maximum(tform_pts[:,1]) - i
  w = maximum(tform_pts[:,2]) - j
  return BoundingBox(i, j, h, w)
end

"""
Convert Tuple for image size into a BoundingBox at (0,0)
"""
function sz2bb(sz)
  return BoundingBox(0, 0, sz[1]-1, sz[2]-1);  # should origin be at 1,1, 0,0, or 0.5,0.5?
end

function test_bb_operations()
  A = BoundingBox(10,10,30,20)
  B = BoundingBox(5,20,20,10)
  C = A+B
  D = BoundingBox(5,10,35,20)
  @test C == D

  @test_throws InexactError BoundingBox(0,0,-1,-1)
end

function test_find_mesh_bb()
    dst = [4.0 2.0;
            8.0 2.0;
            6.0 10.0]
    bounds = find_bounds(dst')
    @test bounds == (3, 1, 9, 11)

    dst = [0.0 0.0;
            0.0 0.0;
            0.0 0.0]
    bounds = find_bounds(dst')
    @test bounds == (-1, -1, 1, 1)
end

function test_tform_bb()
  sz = (100, 100)
  tform = [1 0 0;
          0 1 0;
          100 100 1];
  bb = BoundingBox(100, 100, 100, 100)
  tbb = tform_bb(sz2bb(sz), tform)
  @test bb2pts(bb) == bb2pts(tbb)

  sz = (100, 100)
  tform = [cos(pi/6) -sin(pi/6) 0;
          sin(pi/6) cos(pi/6) 0;
          0 0 1];
  bb = BoundingBox(0, -50, 136.60254037844388, 136.60254037844388)
  tbb = tform_bb(sz2bb(sz), tform)
  @test_approx_eq bb2pts(bb) bb2pts(tbb)
end

function test_snap_bb()
  bb = snap_bb(sz2bb((100.5, 200.75)))
  tbb = sz2bb((101, 201))
  @test bb2pts(bb) == bb2pts(tbb)

  bb = snap_bb(sz2bb((100.9, 200.1)))
  tbb = sz2bb((101, 201))
  @test bb2pts(bb) == bb2pts(tbb)

  bb = snap_bb(BoundingBox(100.9, 200.1, 100.5, 100.2))
  tbb = sz2bb((101, 201))
  tbb = BoundingBox(100, 200, 102, 101)
  @test bb2pts(bb) == bb2pts(tbb)  
end

end