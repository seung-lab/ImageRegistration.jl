using Images
using ImageView
using TestImages
# using AffineTransforms    # incompatible with MATLAB convention
using Base.Test
using Color
import Bounding: BoundingBox, snap_bb, tform_bb

"""
`IMWARP` - Apply affine transform to image using bilinear interpolation

    warped_img, warped_offset = imwarp(img, tform, offset)

* `img`: 2D array, image (todo: extend to Image type)
* `tform`: 3x3 matrix, affine transform (defined as row vector x matrix)
* `offset`: 2-element array, position of img[1,1] in global space
* `warped_img`: with pixel values the same type as original image (for Int type, 
  pixel values are rounded)
* `warped_offset`: 2-element array, position of warped_img[1,1] in global space

The bounding box of the warped image is defined as the smallest
integer-valued rectangle that contains the affine transform of the
bounding box of the original image.

This means that `warped_offset` is constrained to be integer-valued,
though `offset` is allowed to have floating point values.  The integer
constraint removes the need for further interpolation in any
subsequent fusing of multiple warped tiles.

Bilinear interpolation, with extrapolation using zero fill value.

### Definitions

Global position of `img` pixels (analogous definitions for `warped_img`): 
  
* [1,1] pixel has position (offset[1], offset[2]) in global space
* [i,j] pixel has position (offset[1]+i-1, offset[2]+j-1)  

Affine transform of a global position:

* homogeneous coordinates [x, y, 1]  [ax + by + c, dx + ey + f, 1]
* or equivalently [x, y, 1]  [x, y, 1] * tform

        where `tform` = [a d 0;  
                         b e 0;  
                         c f 1]

Note that:

1. transform is *row* vector x matrix (right hand side matrix multiplication)
2. definition compatible with [MATLAB](http://www.mathworks.com/help/images/ref/affine2d-class.html), but not AffineTransforms.jl
3. meaning of transform depends on whether the image is in ij or xy format

Affine transform of an image (two equivalent definitions):

1. The value of `img` at a position is equal to the value of
the `warped_img` at the transformed position. 
2. The value of `warped_img` at a position is equal to the value of
`img` at the inverse transformed position.  

We apply definition 2 as it's compatible with gridding and interpolation.

Bounding box of an image of size (m,n):

* smallest rectangle in global space that contains the positions of the [1,1] 
  and [m,n] pixels
* represented by the 4-tuple (offset[1],offset[2],m-1,n-1)

    (offset[1],offset[2])  ___________
                          |           |
                    height|           |
                     m-1  |           |
                          |___________|
                               n-1   (offset[1]+m-1,offset[2]+n-1)
                              width

""" 
function imwarp{T}(img::Array{T}, tform, offset=[0,0])
  # img bb rooted at offset, with height and width calculated from image
  bb = BoundingBox(offset..., size(img, 1)-1.0, size(img, 2)-1.0)
  # transform original bb to generate new bb (may contain continuous values)
  wbb = tform_bb(bb, tform)
  # snap transformed bb to the nearest exterior integer values
  tbb = snap_bb(wbb)
  # construct warped_img, pixels same Type as img, size calculated from tbb
  # WARNING: should have zero values, but unclear whether guaranteed by similar
  warped_img = similar(img, tbb.h+1, tbb.w+1)
  # Check once if the image type is integer, because we'll need to round
  is_int_image = isa(T, Integer)
  # offset of warped_img from the global origin
  warped_offset = [tbb.i, tbb.j]
  M = inv(tform)
  # cycle through all the pixels in warped_img
  for j = 1:size(warped_img,2)
    for i = 1:size(warped_img,1) # cycle through column-first for speed
        # shift from pixel space to global space before we reverse tform
        # (we index to zero, then add on the offset)
        u, v = i-1+warped_offset[1], j-1+warped_offset[2]
        # transform from warped image to initial image
        # x, y = [u, v, 1] * M - but writing it out moves faster
        x, y = M[1,1]*u + M[2,1]*v + M[3,1], M[1,2]*u + M[2,2]*v + M[3,2]
        # shift back from global space to pixel space for the original img
        # (subtract off its offset, then index up to one)
        x, y = x-offset[1]+1.0, y-offset[2]+1.0

        # Slow...
        #warped_img[i,j] = round(Uint8, bilinear(img, x, y))
        #warped_img[i,j] = bilinear(img, x, y)

        # Bilinear interpolation
        fx, fy = floor(Int64, x), floor(Int64, y)
        wx, wy = x-fx, y-fy
        # if 1 <= fx && fx+1 <= size(img, 1) && 1 <= fy && fy+1 <= size(img, 2)
        inside = true
        if 1 <= fx && fx+1 <= size(img, 1)
            if 1 <= fy && fy+1 <= size(img, 2)   # normal case
            # Expansion of p = [1-wx wx] * img[fx:fx+1, fy:fy+1] * [1-wy; wy]
                p = ((1-wx)*img[fx,fy] + wx*img[fx+1,fy]) * (1-wy) + ((1-wx)*img[fx,fy+1] + wx*img[fx+1,fy+1]) * wy
            elseif fy == size(img, 2) && wy==0   # edge case
                p = (1-wx)*img[fx,fy] + wx*img[fx+1,fy]
            else
                inside=false
            end
        elseif fx == size(img, 1) && wx==0
            if 1 <= fy && fy+1 <= size(img, 2)   # edge case
                p = img[fx,fy] * (1-wy) + img[fx,fy+1] * wy
            elseif fy == size(img, 2) && wy==0  # corner case
                p = img[fx,fy]
            else
                inside=false
            end
        else
            inside=false
        end
        if inside==true
            if is_int_image
              warped_img[i,j] = round(T, p);
            else
              warped_img[i,j] = p
            end
        #else
        #    warped_img[i,j] = 0 # Fill value set to zero based on similar function above
        end
    end
  end
  return warped_img, warped_offset
end

function test_imwarp()
  img = reshape(float(collect(1:121).%2), 11, 11) # 11x11 checkerboard
  warp = img

  tform = [1 0 0;
          0 1 0;
          0 0 1]
  offset = [0, 0]
  img_warped, warped_offset = imwarp(img, tform, offset)
  @test_approx_eq warp img_warped
  @test warped_offset == [0, 0]

  tform = [1 0 0;
          0 1 0;
          0 0 1]
  offset = [10, 20]
  img_warped, warped_offset = imwarp(img, tform, offset)
  @test_approx_eq warp img_warped
  @test warped_offset == [10, 20]

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
  offset = [0.0, 0.0]
  img_warped, warped_offset = imwarp(img, tform, offset)
  @test warped_offset == [0, -5]
end

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
