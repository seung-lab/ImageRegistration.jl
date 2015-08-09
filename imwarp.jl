using Images
using ImageView
using AffineTransforms
using Grid
using Base.Test

type Rect
  i::Float64
  j::Float64
  h::Float64
  w::Float64
end

type SpatialRef
  i::Float64
  j::Float64
end

function rect2pts(r)
  # Create polygon point list, counter-clockwise
  return [r.i r.j;
          r.i+r.h r.j;
          r.i+r.h r.j+r.w;
          r.i r.j+r.w;
          r.i r.j];
end

# Apply geometric transform to image
function imwarp(img, tform)
  bb = sz2bb(size(img))
  tbb = snap_bb(tform_bb(bb, tform))
  out_img = Array(Int64, Int64(tbb.h), Int64(tbb.w))

  i = 1:size(img,1)
  j = 1:size(img,2)
  zi = CoordInterpGrid((i,j), img.data, 0.0, InterpLinear)

  i_off = Int64(tbb.i)
  j_off = Int64(tbb.j)

  for i = 1:Int64(bb.h)
    for j = 1:Int64(bb.w)
      pt = [i+i_off j+j_off 1]
      ipt = (pt*tform^-1)[:,1:2]
      out_img[i, j] = zi[ipt...]
    end
  end  
  return out_img
end

function snap_bb(bb)
  r = rect2pts(bb)
  bb.i = floor(bb.i)
  bb.j = floor(bb.j)
  h = ceil(r[2,1]) - bb.i
  w = ceil(r[3,2]) - bb.j
  return Rect(bb.i, bb.j, h, w)
end

function tform_bb(bb, tform)
  bb_pts = hcat(rect2pts(bb), ones(Int64, 5,1))
  tform_pts = bb_pts * tform
  i = minimum(tform_pts[:,1])
  j = minimum(tform_pts[:,2])
  w = maximum(tform_pts[:,2]) - j
  h = maximum(tform_pts[:,1]) - i
  return Rect(i, j, w, h)
end

function sz2bb(sz)
  return Rect(0, 0, sz...);
end

function test_tform_bb()
  sz = (100, 100)
  tform = [1 0 0;
          0 1 0;
          100 100 1];
  bb = Rect(100, 100, 100, 100)
  tbb = tform_bb(sz2bb(sz), tform)
  @test rect2pts(bb) == rect2pts(tbb)

  sz = (100, 100)
  tform = [cos(pi/6) -sin(pi/6) 0;
          sin(pi/6) cos(pi/6) 0;
          0 0 1];
  bb = Rect(0, -50, 136.60254037844388, 136.60254037844388)
  tbb = tform_bb(sz2bb(sz), tform)
  @test_approx_eq rect2pts(bb) rect2pts(tbb)
end

function test_snap_bb()
  bb = snap_bb(sz2bb((100.5, 200.75)))
  tbb = sz2bb((101, 201))
  @test rect2pts(bb) == rect2pts(tbb)

  bb = snap_bb(sz2bb((100.9, 200.1)))
  tbb = sz2bb((101, 201))
  @test rect2pts(bb) == rect2pts(tbb)

  bb = snap_bb(Rect(100.9, 200.1, 100.5, 100.2))
  tbb = sz2bb((101, 201))
  tbb = Rect(100, 200, 102, 101)
  @test rect2pts(bb) == rect2pts(tbb)  
end

function test_img()
  imread("test_images/turtle.png")
end

function test()
  test_tform_bb()
  test_snap_bb()
end
