# function test_incidence_to_triangles()
D = [1 -1 0 0 0;
   1 0 -1 0 0;
   0 1 -1 0 0;
   0 1 0 -1 0;
   0 1 0 0 -1;
   0 0 0 1 -1];
triangles = [1 2 3;
      2 4 5];
tri = incidence_to_triangles(D)
@test triangles == tri

D = [1 -1 0 0 0;
   1 0 -1 0 0;
   1 0 0 -1 0;
   1 0 0 0 -1;
   0 1 -1 0 0;
   0 1 0 -1 0;
   0 1 0 0 -1;
   0 0 1 -1 0;
   0 0 1 0 -1;
   0 0 0 1 -1];
triangles = [1 2 3;
      1 2 4;
      1 2 5;
      1 3 4;
      1 3 5;
      1 4 5;
      2 3 4;
      2 3 5;
      2 4 5;
      3 4 5];
tri = incidence_to_triangles(D)
@test triangles == tri

D = [1 -1 0 0 0;
   0 0 1 -1 0;
   0 0 1 0 -1;
   1 0 -1 0 0;
   1 0 0 -1 0;
   1 0 0 0 -1;
   0 1 -1 0 0;
   0 1 0 -1 0;
   0 1 0 0 -1;
   0 0 0 1 -1];
triangles = [1 2 3;
      1 2 4;
      1 2 5;
      1 3 4;
      1 3 5;
      1 4 5;
      2 3 4;
      2 3 5;
      2 4 5;
      3 4 5];
tri = incidence_to_triangles(D)
@test triangles == tri

# function test_incidence_to_dict()
D = [1 -1 0 1 0;
   0 0 1 -1 0];
node_dict = Dict(1 => Set(2), 3 => Set(4))
nd = incidence_to_dict(D)
@test node_dict == nd

D = [1 -1 0 1 0;
   1 -1 1 1 0;
   1 0 1 1 0];
node_dict = Dict(1 => Set([2, 3]))
nd = incidence_to_dict(D)
@test node_dict == nd 