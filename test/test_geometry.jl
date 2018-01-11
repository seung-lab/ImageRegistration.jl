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

a = [3,4]
b = [4,5]
@test cross2(a, b) == -1
@test cross2(b, a) == 1

a = [0,0,2,2]
b = [0,2,2,0]
@test find_line_segment_intersection(a,b) == [1,1]
@test find_line_intersection(a,b) == [1,1]
b = [0,4,2,2]
@test find_line_segment_intersection(a,b) == [2,2]
@test find_line_intersection(a,b) == [2,2]
b = [0,10,2,8]
@test find_line_segment_intersection(a,b) == nothing
@test find_line_intersection(a,b) == [5,5]
a = [1,0,1,2]
b = [0,2,2,0]
@test find_line_segment_intersection(a,b) == [1,1]
@test find_line_intersection(a,b) == [1,1]
a = [1,0,1,2]
b = [1,2,1,1]
@test find_line_segment_intersection(a,b) == [1,1]
@test find_line_intersection(a,b) == [1,1]

a = [0,0,3,0]
b = [4,5]
@test pt_on_left(b, a) == true
b = [1000,1000]
@test pt_on_left(b, a) == true
b = [-1,-1]
@test pt_on_left(b, a) == false
a = [3,0,0,0]
@test pt_on_left(b, a) == true
b = [4,5]
@test pt_on_left(b, a) == false

subject = [0.0 0.0
           10.0 0.0
           10.0 10.0
           0.0 10.0
           0.0 0.0]
clip = [-1.0 3.0
      11.0 3.0
      11.0 8.0
      -1.0 8.0
      -1.0 3.0]
output = [10.0 3.0
      10.0 8.0
      0.0 8.0
      0.0 3.0
      10.0 3.0]
@test_approx_eq clip_polygon(subject, clip) output

subject = [0.0 0.0
           10.0 0.0
           10.0 10.0
           0.0 10.0
           0.0 0.0]
clip = [3.0 3.0
      8.0 3.0
      8.0 8.0
      3.0 8.0
      3.0 3.0]
output = [3.0 3.0
      8.0 3.0
      8.0 8.0
      3.0 8.0
      3.0 3.0]
@test_approx_eq clip_polygon(subject, clip) output

a = [0.0 0.0
    4.0 0.0
    4.0 3.0
    0.0 0.0]
@test_approx_eq poly_area(a) 6
a = [0.0 0.0
    -4.0 0.0
    -4.0 -3.0
    0.0 0.0]
@test_approx_eq poly_area(a) 6