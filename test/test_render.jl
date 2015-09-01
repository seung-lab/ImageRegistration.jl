using Base.Test

# test_padimage()
o = ones(5,2)
po = padimage(o, 0, 0, 0, 0)
@test size(po, 1) == 5
@test size(po, 2) == 2

o = ones(5,2)
po = padimage(o, 1, 2, 3, 4)
@test size(po, 1) == 11
@test size(po, 2) == 6
@test po[1,1] == 0

o = ones(5,2)
po = padimage(o, 0, 1, 1, 0)
@test size(po, 1) == 6
@test size(po, 2) == 3
@test po[6,3] == 0    

o = convert(Array{Int64, 2}, ones(5,2))
po = padimage(o, 1, 2, 3, 4)
@test size(po, 1) == 11
@test size(po, 2) == 6

# test_padimage()
img = ones(1,1)
offset = [0,0]
bb = BoundingBox(0,0,1,1)
@test padimage(img, offset, bb) == img

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(1,1,1,1)
@test padimage(img, offset, bb) == img

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(0,0,2,2)
@test padimage(img, offset, bb) == [0 0; 0 1]

img = ones(1,1)
offset = [1,1]
bb = BoundingBox(0,0,3,3)
@test padimage(img, offset, bb) == [0 0 0; 0 1 0; 0 0 0]

img = ones(1,1)
offset = [0,0]
bb = BoundingBox(-1,-1,2,2)
@test padimage(img, offset, bb) == [0 0; 0 1]

# test_imfuse()
A = rand(5,5)
B = rand(5,5)
BB_A = [0, 0]
BB_B = [0, 0]
O, BB_O = imfuse(A, BB_A, B, BB_B)
@test O.channels[1] == A
@test O.channels[2] == B
@test BB_O == BB_A
@test BB_O == BB_B

A = rand(5,5)
B = rand(5,5)
BB_A = [0, 0]
BB_B = [2, 3]
O, BB_O = imfuse(A, BB_A, B, BB_B)
@test O.channels[1][1:5, 1:5] == A
@test O.channels[2][3:7, 4:8] == B
@test BB_O == BB_A
