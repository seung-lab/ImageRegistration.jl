# test_bb_operations()
A = BoundingBox(10,10,30,20)
B = BoundingBox(5,20,20,10)
C = A+B
D = BoundingBox(5,10,35,20)
@test C == D

# @test_throws InexactError BoundingBox(0,0,-1,-1)

# test_find_mesh_bb()
  dst = [4.0 2.0;
          8.0 2.0;
          6.0 10.0]
  bounds = find_mesh_bb(dst)
  @test bounds == BoundingBox(4,2,4,8)

  dst = [0.0 0.0;
          0.0 0.0;
          0.0 0.0]
  bounds = find_mesh_bb(dst)
  @test bounds == BoundingBox(0,0,0,0)

# test_tform_bb()
tform = [1 0 0;
        0 1 0;
        0 0 1];
bb = BoundingBox(100,100,100,100)
tbb = tform_bb(bb, tform)
@test bb == tbb

tform = [1 0 0;
        0 1 0;
        100 100 1];
bb = BoundingBox(100,100,100,100)
tbb = tform_bb(bb, tform)
@test BoundingBox(200,200,100,100) == tbb

sz = (100, 100)
tform = [cos(pi/6) -sin(pi/6) 0;
        sin(pi/6) cos(pi/6) 0;
        0 0 1];
bb = BoundingBox(0, -50, 136.60254037844388, 136.60254037844388)
tbb = tform_bb(sz_to_bb(sz), tform)
@test_approx_eq bb_to_pts(bb) bb_to_pts(tbb)

# test_snap_bb()
bb = snap_bb(sz_to_bb((100.5, 200.75)))
tbb = sz_to_bb((101, 201))
@test bb_to_pts(bb) == bb_to_pts(tbb)

bb = snap_bb(sz_to_bb((100.9, 200.1)))
tbb = sz_to_bb((101, 201))
@test bb_to_pts(bb) == bb_to_pts(tbb)

bb = snap_bb(BoundingBox(100.9, 200.1, 100.5, 100.2))
tbb = sz_to_bb((101, 201))
tbb = BoundingBox(100, 200, 102, 101)
@test bb_to_pts(bb) == bb_to_pts(tbb)  

a = BoundingBox(0,0,10,10)
b = BoundingBox(0,0,5,5)
@test a+b == BoundingBox(0,0,10,10)
@test a-b == BoundingBox(0,0,5,5)

a = BoundingBox(2,2,10,10)
b = BoundingBox(0,0,5,5)
@test a+b == BoundingBox(0,0,12,12)
@test a-b == BoundingBox(2,2,3,3)

a = BoundingBox(20,20,10,10)
b = BoundingBox(0,0,5,5)
@test a+b == BoundingBox(0,0,30,30)
@test isnan((a-b).i)

a = [0.0 0.0;
     0.0 10.0;
     10.0 10.0;
     10.0 0.0;
     0.0 0.0];
@test in_poly([5,5], a)
@test in_poly([0,0], a)
@test !in_poly([-1,0], a)
@test in_poly([1,1], a)
@test !in_poly([11,1], a)
@test in_poly([10,1], a)
@test in_poly([9.99,1], a)

b = [5.0 0.0;
     9.0 10.0;
     3.0 11.0;
     4.0 5.0];
@test in_poly([5,1], b)
@test !in_poly([0,0], b)
@test !in_poly([-1,0], b)
@test !in_poly([1,1], b)
@test !in_poly([11,1], b)
@test !in_poly([10.01,1], b)
@test !in_poly([9.99,1], b)
@test in_poly([4,5], b)
@test in_poly([4,6], b)

@test poly_intersects(a, b)
c = [0.0 0.0;
     3.0 0.0;
     3.0 4.0]
d = [3.0 0.0;
     3.0 4.0;
     6.0 0.0]
@test poly_intersects(c, d)
@test poly_intersects(d, c)
@test !poly_intersects(b, c)
@test !poly_intersects(c, b)
@test poly_intersects(c, a)
@test !poly_intersects(b, d)
@test poly_intersects(a, d)

@test on_line_segment([0,0], ([-1,0], [1,0]))
@test on_line_segment([0,0], ([-1,0], [1,0.01]))
@test on_line_segment([0,0], ([-1,0], [1.0001,0]))
@test on_line_segment([0,0], ([-1,0], [1.001,0]))
@test on_line_segment([0,0], ([-1,0], [1.01,0]))
@test on_line_segment([0,0], ([-1,0], [1.1,0]))
@test !on_line_segment([0,0], ([-1,0], [1.1,0.1]))

@test on_line_segment([0,0], ([0,-1], [0,1]))
@test on_line_segment([0,0], ([0,-1], [0.01,1]))
@test on_line_segment([0,0], ([0,-1], [0,1.0001]))
@test on_line_segment([0,0], ([0,-1], [0,1.001]))
@test on_line_segment([0,0], ([0,-1], [0,1.01]))
@test on_line_segment([0,0], ([0,-1], [0,1.1]))
@test !on_line_segment([0,0], ([0,-1], [0.1,1.1]))

@test on_line_segment([1,1], ([0,0], [10,10]))

e = [3.0 3.0;
     6.0 3.0;
     6.0 6.0;
     3.0 6.0;
     3.0 3.0];
@test poly_intersects(a, e)
@test poly_intersects(e, a)
