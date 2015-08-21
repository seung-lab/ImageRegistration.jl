using Images
using ImageView
using TestImages
# using AffineTransforms    # incompatible with MATLAB convention
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

## Definitions
Bounding boxes of images contain the image by running through the upper 
left corner pixel and lower right corner pixel.

    (1,1)  ___________
          |           |
    height|           |
    (i) m |           |
          |___________|
                n   (m,n)
               (j)
              width
Offset between an internal point, A (i,j), and the upper leftcorner of the 
image will be (i-1, j-1).

Coordinates are in i,j format, to match Julia's column-first storage.

Affine transform matrix form

    A = [a b 0;
         c d 0;
         x y 1] --> homogoneous coordinates are in 3rd column (A[:,3])
                    so, right-hand matrix multiplication:
                      [x, y, 1] * M

Args:

* img: 2D array of an image
* tform: 3x3 array for an affine transform for right-handed multiplication
* offset: 2-element array for displacement of the upper left-hand pixel of img
  from the global origin (optional)

Returns:

* warped_img: the new image transformed from the input image 
    (for Int type images, the output is the same type with pixel values rounded)
* offset: updated 2-element array for displacement of image from global origin

    warped_img, warped_offset = imwarp(img, tform, offset)
"""
function imwarp{T}(img::Array{T}, tform, offset=[0.0,0.0])
  # img bb set at offset, with width and height matching the image
  bb = BoundingBox(offset..., size(img, 1)-1, size(img, 2)-1)
  # transform initial bb, and generate new bb (may contain continous values)
  wbb = tform_bb(bb, tform)
  # snap transformed bb to the nearest exterior int values for pixel handling
  tbb = snap_bb(wbb)
  # set warped_img to be zeros (same Type as img), with dimensions of tbb
  warped_img = similar(img, Int64(tbb.h)+1, Int64(tbb.w)+1)

  # offset of warped_img from the global origin
  warped_offset = [Int64(tbb.i), Int64(tbb.j)]
  M = inv(tform)
  # cycle through all the pixels in warped_img
  for j = 1:size(warped_img,2)
    for i = 1:size(warped_img,1) # cycle through column-first for speed
        # shift from pixel space to global space before we reverse tform
        # (we index to zero, then add on the offset)
        u, v = i-1+warped_offset[1], j-1+warped_offset[2]
        # transform from warped image to initial image
        # x, y = [u, v, 1] * M --> but writing it out moves faster
        x, y = M[1,1]*u + M[2,1]*v + M[3,1], M[1,2]*u + M[2,2]*v + M[3,2]
        # shift back from global space to pixel space for the original img
        # (subtract off its offset, then index up to one)
        x, y = x-offset[1]+1, y-offset[2]+1

        # Slow...
        #warped_img[i,j] = round(Uint8, bilinear(img, x, y))
        #warped_img[i,j] = bilinear(img, x, y)

        # Bilinear interpolation
        fx, fy = floor(Int64, x), floor(Int64, y)
        wx, wy = x-fx, y-fy
        if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
            # Expansion of p = [1-wx wx] * img[fx:fx+1, fy:fy+1] * [1-wy; wy]
            p = ((1-wx)*img[fx,fy] + wx*img[fx+1,fy]) * (1-wy) + ((1-wx)*img[fx,fy+1] + wx*img[fx+1,fy+1]) * wy
            if isa(warped_img[1], Integer)
              warped_img[i,j] = round(T, p);
            else
              warped_img[i,j] = p
            end
        #else
        #    warped_img[i,j] = 0 # Fill value set to zero based on similar function above
        end
    end
  end
  return warped_img, [tbb.i, tbb.j], bb, wbb, tbb
end

function test_imwarp()
  img = reshape(float(collect(1:121).%2), 11, 11) # 11x11 checkerboard
  warp = img
  warp[11,:] = zeros(11)
  warp[:,11] = zeros(11)

  tform = [1 0 0;
          0 1 0;
          0 0 1]
  offset = [0, 0]
  img_warped, warped_offset = imwarp(img, tform, offset)
  @test_approx_eq warp img_warped
  @test warped_offset == [0, 0]

  tform = [1 0 0;
          0 1 0;
          10 10 1]
  offset = [0, 0]
  img_warped, warped_offset = imwarp(img, tform, offset)
  @test_approx_eq warp img_warped
  @test warped_offset == [10, 10]

  tform = [cos(0.5) -sin(0.5) 0;
          sin(0.5) cos(0.5) 0;
          0 0 1]
  offset = [0, 0]
  img_warped, warped_offset = imwarp(img, tform, offset)
  @test warped_offset == [0, -6]
end
