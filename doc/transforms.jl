"""
`calculate_rigid` - Calculate left-hand rigid transform from two N x d point sets

Inputs:
  moving_pts, fixed_pts:
    2D arrays with each row being coords of a point of dimension d-1 (homogeneous coords) or d.
Returns:
  tform: 
    Rigid transformation T that produces  fixed_pts = moving_pts * T
    (if moving_pts and fixed_pts are in homogeneous coords)
"""
function calculate_rigid(moving_pts, fixed_pts)
  n, dim = size(moving_pts)
  moving_bar = mean(moving_pts, 1)
  fixed_bar = mean(points2, 1)
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
"""
function calculate_affine(moving_pts, fixed_pts)
  moving = hcat(moving_pts, ones(size(moving_pts,1)))
  fixed = hcat(fixed_pts, ones(size(fixed_pts,1)))
  return moving \ fixed
end

"""
"""
function calculate_translation(moving_pts, fixed_pts)
  moving_bar = mean(moving_pts, 1)
  fixed_bar = mean(points2, 1)
  t = fixed_bar - moving_bar
  return [eye(2) zeros(dim); t 1]
end