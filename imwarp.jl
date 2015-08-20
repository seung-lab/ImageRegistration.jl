using Images
using ImageView
using TestImages
using AffineTransforms
using Base.Test
using Color
include("BoundingBox.jl")

# Not used - slow.
function bilinear(img, x, y)
    fx, fy = floor(Int64, x), floor(Int64, y)
    wx, wy = x-fx, y-fy
    if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
      # Expansion of p = [1-wx wx] * img[fx:fx+1, fy:fy+1] * [1-wy; wy]
      p = ((1-wx)*img[fx,fy] + wx*img[fx+1,fy]) * (1-wy) + ((1-wx)*img[fx,fy+1] + wx*img[fx+1,fy+1]) * wy
      return p
    else
      0
    end
end

"""
Apply geometric transform to image

Args:

* img: 2D array of an image
* tform: 3x3 array for an affine transform of the form

  [a b 0;
   c d 0;
   x y 1]

Returns:

* out_img: the new image transformed from the input image (for Integer type images, 
  the output is the same type with pixel values rounded)
* tbb: BoundingBox of out_img (origin relative to img)

  out_img, tbb = imwarp(img, tform)
"""
function imwarp{T}(img::Array{T}, tform)
    bb = sz2bb(size(img))
    tbb = snap_bb(tform_bb(bb, tform))
    out_img = similar(img, Int64(tbb.h), Int64(tbb.w))

    i_off = Int64(tbb.i)
    j_off = Int64(tbb.j)
    M = inv(tform);
    for j = 1:Int64(tbb.w)
      for i = 1:Int64(tbb.h)
          u, v = i+i_off, j+j_off
          # x, y = M * [u, v, 1]
          x, y = M[1,1]*u + M[1,2]*v + M[1,3], M[2,1]*u + M[2,2]*v + M[2,3]

          # Slow...
          #out_img[i,j] = round(Uint8, bilinear(img, x, y))
          #out_img[i,j] = bilinear(img, x, y)

          # Bilinear interpolation
          fx, fy = floor(Int64, x), floor(Int64, y)
          wx, wy = x-fx, y-fy
          if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
              # Expansion of p = [1-wx wx] * img[fx:fx+1, fy:fy+1] * [1-wy; wy]
              p = ((1-wx)*img[fx,fy] + wx*img[fx+1,fy]) * (1-wy) + ((1-wx)*img[fx,fy+1] + wx*img[fx+1,fy+1]) * wy
              if isa(out_img[1], Integer)
                out_img[i,j] = round(T, p);
              else
                out_img[i,j] = p
              end
          #else
          #    out_img[i,j] = 0
          end
      end
    end
    return out_img, tbb
end
