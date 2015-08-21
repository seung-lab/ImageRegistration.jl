using Base.Test

type BoundingBox
  x::Float64
  y::Float64
  h::Float64
  w::Float64
end

BoundingBox() = BoundingBox(0,0,0,0)

"""
Add bounding boxes by finding BoundingBox that encompasses both
"""
function +(bbA::BoundingBox, bbB::BoundingBox)
  x = min(bbA.x, bbB.x)
  y = min(bbA.y, bbB.y)
  w = max(bbA.w+bbA.x, bbB.w+bbB.x) - x
  h = max(bbA.h+bbA.y, bbB.h+bbB.y) - y
  return BoundingBox(x,y,w,h)
end

"""
Test if BoundingBoxes have the same origin and dimensions
"""
function ==(bbA::BoundingBox, bbB::BoundingBox)
  return bbA.x == bbB.x && bbA.y == bbB.y && bbA.w == bbB.w && bbA.h == bbB.h
end

"""
Convert bounding box object to a polygon point list (counter-clockwise)
"""
function bb2pts(r)
  return [r.x r.y;
          r.x+r.h r.y;
          r.x+r.h r.y+r.w;
          r.x r.y+r.w;
          r.x r.y];
end

"""
Return xmin, ymin, xmax, ymax of a bounding box (opposing corner definition)
"""
function minsandmax(bb)
  return bb.x, bb.y, bb.x+bb.w, bb.y+bb.h
end

"""
Snap bounding box to the nearest exterior pixels
"""
function snap_bb(bb)
  r = bb2pts(bb)
  bb.x = floor(bb.x)
  bb.y = floor(bb.y)
  h = ceil(r[2,1]) - bb.x
  w = ceil(r[3,2]) - bb.y
  return BoundingBox(bb.x, bb.y, h, w)
end

"""
Apply affine transform to points in a bounding box & find new bounding box
"""
function tform_bb(bb, tform)
  tform_pts = tform * [bb2pts(bb) ones(size(bb2pts(bb),1),1)]'
  i = minimum(tform_pts[1,:])
  j = minimum(tform_pts[2,:])
  w = maximum(tform_pts[2,:]) - j
  h = maximum(tform_pts[1,:]) - i
  return BoundingBox(i, j, w, h)
end

"""
Convert Tuple for image size into a BoundingBox at (0,0)
"""
function sz2bb(sz)
  return BoundingBox(0, 0, sz...);  # should origin be at 1,1, 0,0, or 0.5,0.5?
end

"""
Find the extrema of a mesh, and generate a bounding box

Args:

* nodes: 2xN array of coordinates for a mesh

Returns:

* BoundingBox containing all of the mesh nodes

    BoundingBox(xlow, ylow, xlow+xhigh, ylow+yhigh) = find_bounds(nodes)
"""
function find_mesh_bb(nodes)
    xlow = Int64(floor(minimum(nodes[:,1])))-1
    ylow = Int64(floor(minimum(nodes[:,2])))-1
    xhigh = Int64(ceil(maximum(nodes[:,1])))
    yhigh = Int64(ceil(maximum(nodes[:,2])))
    return BoundingBox(xlow, ylow, xhigh-xlow, yhigh-ylow)
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

function test_img()
    tfm=tformrotate(0.2)
    img=testimage("mandrill")
    view(imwarp(img,tfm),xy=["y","x"])
end

function test()
  test_bb_operations()
  test_find_mesh_bb()
  test_tform_bb()
  test_snap_bb()
end