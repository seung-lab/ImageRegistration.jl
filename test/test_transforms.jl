# function test_affine
moving_pts = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
fixed_pts = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
tform = calculate_affine(moving_pts, fixed_pts)
@test_approx_eq tform eye(3)

tform = calculate_rigid(moving_pts, fixed_pts)
@test_approx_eq tform eye(3)

tform = calculate_translation(moving_pts, fixed_pts)
@test_approx_eq tform eye(3)

moving_pts = [1.0 1.0;
			2.0 3.0;
			3.0 3.0]
fixed_pts = [5.0 1.0;
			6.0 3.0;
			7.0 3.0]
shift = [1 0 0; 0 1 0; 4 0 1]
tform = calculate_affine(moving_pts, fixed_pts)
@test_approx_eq tform shift

tform = calculate_rigid(moving_pts, fixed_pts)
@test_approx_eq tform shift

tform = calculate_translation(moving_pts, fixed_pts)
@test_approx_eq tform shift