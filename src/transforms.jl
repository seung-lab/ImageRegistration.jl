"""
`CALCULATE_AFFINE` - Compute left-hand affine transform from two Nxd point sets

```
tform = calculate_affine(moving_pts, fixed_pts)
```

* moving_pts: Nxd array, each row being homogeneous coords of a point, of the 
    point set that is moving to fit the fixed point set.
* fixed_pts: Nxd array, each row being homogeneous coords of a point, of the 
    reference point set.
* tform: affine transformation that produces a mapping of the moving_pts to 
    fixed_pts.

```
fixed_pts = moving_pts * tform
```
"""
function calculate_affine(moving_pts, fixed_pts)
  return moving_pts \ fixed_pts
end

"""
`CALCULATE_RIGID` - Compute right-hand rigid transform from two Nxd point sets

```
tform = calculate_rigid(moving_pts, fixed_pts)
```

* moving_pts: Nxd array, each row being homogeneous coords of a point, of the 
    point set that is moving to fit the fixed point set.
* fixed_pts: Nxd array, each row being homogeneous coords of a point, of the 
    reference point set.
* tform: rigid transformation that produces a mapping of the moving_pts to 
    fixed_pts.

```
fixed_pts = moving_pts * tform
```
"""
function calculate_rigid(moving_pts, fixed_pts)
  moving_pts = moving_pts[:,1:end-1]
  fixed_pts = fixed_pts[:,1:end-1]
  n, dim = size(moving_pts)
  moving_bar = mean(moving_pts, 1)
  fixed_bar = mean(fixed_pts, 1)
  moving_centered = moving_pts .- moving_bar
  fixed_centered = fixed_pts .- fixed_bar
  C = fixed_centered.' * moving_centered / n
  U,s,V = svd(C)
  # Rotation R in least squares sense: 
  # moving_pts - moving_bar = (fixed_pts - fixed_bar)*R
  R = (U * diagm(vcat(ones(dim-1), det(U*V.'))) * V.' ).'
  t = fixed_bar - moving_bar*R
  return [R zeros(dim); t 1]
end

"""
`CALCULATE_TRANSLATION` - Compute left-hand translation from two Nxd point sets

```
tform = calculate_translation(moving_pts, fixed_pts)
```

* moving_pts: Nxd array, each row being nonhomogeneous coords of a point, of the 
    point set that is moving to fit the fixed point set.
* fixed_pts: Nxd array, each row being nonhomogeneous coords of a point, of the 
    reference point set.
* tform: translation transformation that produces a mapping of the moving_pts to 
    fixed_pts.

```
fixed_pts = moving_pts * tform
```
"""
function calculate_translation(moving_pts, fixed_pts)
  n, dim = size(moving_pts)
  moving_bar = mean(moving_pts, 1)
  fixed_bar = mean(fixed_pts, 1)
  t = fixed_bar - moving_bar
  return [eye(2) zeros(dim); t 1]
end

function calculate_affine(matches::Matches)
  return calculate_affine(matches.src_points, matches.dst_points)
end

function calculate_rigid(matches::Matches)
  return calculate_rigid(matches.src_points, matches.dst_points)
end

function calculate_translation(matches::Matches)
  return calculate_translation(matches.src_points, matches.dst_points)
end

function calculate_affine(mesh::Mesh)
  return calculate_affine(mesh.src_nodes, mesh.dst_nodes)
end

function calculate_rigid(mesh::Mesh)
  return calculate_rigid(mesh.src_nodes, mesh.dst_nodes)
end

function calculate_translation(mesh::Mesh)
  return calculate_translation(mesh.src_nodes, mesh.dst_nodes)
end

"""
Matrix for rotation in degrees (clockwise rotation with left-handed Cartesian)
"""
function make_rotation_matrix(theta, image_size = nothing)
  angle = deg2rad(theta)
  tform = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1]
  tform[abs(tform) .< eps] = 0
  if image_size != nothing
    	rotated_bb = snap_bb(tform_bb(sz_to_bb([image_size[1], image_size[2]]), tform))
    	rotation_offset = [rotated_bb.i, rotated_bb.j]
	translation_matrix = make_translation_matrix(-rotation_offset);
	tform = tform * translation_matrix
  end
  return tform
end

function make_translation_matrix(offset)
  return [1 0 0; 0 1 0; offset[1] offset[2] 1]
end

function make_scale_matrix(scale)
  return [scale 0 0; 0 scale 0; 0 0 1]
end
