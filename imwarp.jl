using Images
using ImageView
using AffineTransforms
using Grid
using Base.Test

type Rect
  x::Float64
  y::Float64
  h::Float64
  w::Float64
end

type SpatialRef
  x::Float64
  y::Float64
end

function rect2pts(rect::Rect)
  pts = [rect.x rect.y;
        rect.x + rect.h rect.y;
        rect.x + rect.h rect.y + rect.w;
        rect.x rect.y + rect.w];
  return pts
end

function transform_images(imgA, imgB, tformA, tformB)

end

# Apply geometric transform to image
function imwarp(img, tform)
  img_size = size(img)
  data = data(img)


end

function warp_bb(bb, tform)
  rect_pts = hcat(rect2pts(bb), ones(Int64, 4,1))
  tform_pts = rect_pts * tform
  x = minimum(tform_pts[:, 1])
  y = minimum(tform_pts[:, 2])
  h = maximum(tform_pts[:, 1]) - x
  w = maximum(tform_pts[:, 2]) - y
  return Rect(x, y, h, w)
end

function sz2bb(sz)
  return Rect(0, 0, sz...)
end

function test_warp_bb()
  sz = (100, 100)
  tform = [1 0 0;
          0 1 0;
          100 100 1];
  bb = Rect(100, 100, 100, 100)
  tbb = warp_bb(sz2bb(sz), tform)
  @test rect2pts(bb) == rect2pts(tbb)

  sz = (100, 100)
  tform = [cos(pi/6) -sin(pi/6) 0;
          sin(pi/6) cos(pi/6) 0;
          0 0 1];
  bb = Rect(0, -50, 136.6, 136.6)
  tbb = warp_bb(sz2bb(sz), tform)
  @test rect2pts(bb) == rect2pts(tbb)
end

function test_img()
  imread("test_images/turtle.png")
end

function test()
  test_warp_bb()
end
