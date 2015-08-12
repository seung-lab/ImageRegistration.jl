using Images
using ImageView
using TestImages
using Interpolations
using AffineTransforms
using Grid
using Base.Test
using Color
using FixedPointNumbers

type Rect
  i::Float64
  j::Float64
  h::Float64
  w::Float64
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
    out_img = similar(img, Int64(tbb.h), Int64(tbb.w))

    zi = interpolate(img, BSpline(Constant),OnGrid)
#    zi = CoordInterpGrid((1:size(img,1),1:size(img,2)), img, 0.0, InterpLinear)

    i_off = Int64(tbb.i)
    j_off = Int64(tbb.j)
    II=[i for i = 1:Int64(tbb.h), j = 1:Int64(tbb.w)][:]'
    JJ=[j for i = 1:Int64(tbb.h), j = 1:Int64(tbb.w)][:]'
    pt = [II+i_off; JJ+j_off]
    ipt = tforminv(tform,pt)
    for t=1:length(II)
        out_img[II[t], JJ[t]] = zi[ipt[:,t]...]
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
  tform_pts = tformfwd(tform,rect2pts(bb)');
  i = minimum(tform_pts[1,:])
  j = minimum(tform_pts[2,:])
  w = maximum(tform_pts[2,:]) - j
  h = maximum(tform_pts[1,:]) - i
  return Rect(i, j, w, h)
end

function sz2bb(sz)
  return Rect(0, 0, sz...);  # should origin be at 1,1 or 0,0?
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
    tfm=tformrotate(0.2)
    img=testimage("mandrill")
    view(imwarp(img,tfm),xy=["y","x"])
end

function test()
  test_tform_bb()
  test_snap_bb()
end
