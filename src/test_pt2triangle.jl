# function test_pt2triangle()
pt = [0 0]
verts = [0 0; 3 0; 3 4]
triangles = [1 2 3]
t = pt2triangle(pt, verts, triangles)
@test_approx_eq t[1] 1
@test t[2] == tuple(1,0,0)

pt = [1.5/3+5 1/3]
verts = [0 0; 3 0; 3 4; 5 0; 6 0; 5.5 1]
triangles = [1 2 3; 4 5 6]
t = pt2triangle(pt, verts, triangles)
@test_approx_eq t[1] 2
@test_approx_eq t[2][1] 1/3
@test_approx_eq t[2][2] 1/3
@test_approx_eq t[2][3] 1/3

# function test_xy2barycentric()
pt = [0 0]
verts = [0 0; 3 0; 3 4]
tw = xy2barycentric(pt, verts)
w = [1 0 0]
@test_approx_eq w tw

pt = [1.5/3 1/3]
verts = [0 0; 1 0; 0.5 1]
tw = xy2barycentric(pt, verts)
w = [1/3 1/3 1/3]
@test_approx_eq w tw

#function test_verts2triangle()
verts = [0 0; 3 0; 3 4]
t = verts2triangle(verts)
@test getx(geta(t)) == 1.0
@test gety(geta(t)) == 1.0
@test getx(getb(t)) == 2.0
@test gety(getb(t)) == 1.0