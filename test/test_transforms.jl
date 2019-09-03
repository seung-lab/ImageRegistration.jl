using LinearAlgebra
println("testing transforms")

# function test_affine
moving_pts = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
fixed_pts = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
tform = calculate_affine(moving_pts, fixed_pts)
@test tform ≈ Matrix(1.0I,2,2)

tform = calculate_rigid(moving_pts, fixed_pts)
@test tform ≈ Matrix(1.0I,2,2)

tform = calculate_translation(moving_pts, fixed_pts)
@test tform ≈ Matrix(1.0I,3,3)

moving_pts = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
fixed_pts = [5.0 1.0;
			6.0 3.0;
			7.0 3.0]
shift = [1 0 0; 0 1 0; 4 0 1]

tform = calculate_translation(moving_pts, fixed_pts)
@test tform ≈ shift
