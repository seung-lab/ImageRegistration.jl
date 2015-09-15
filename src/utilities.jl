"""
Count number of displacement vector lengths k-sigma from the group mean

Args:

* vectors: 4xN array of points (start and end points of displacement vectors)
* k: number of sigmas to set the threshold for counting

Returns:

* integer for number of displacement vectors above the k-sigma threshold

  num = count_outliers(vectors, k)
"""
function count_outliers(vectors, k)
  d = sum((vectors[1:2,:] - vectors[3:4,:]).^2, 1).^(1/2)
  d_mean = mean(d)
  d_std = mean((d.-d_mean).^2).^(1/2)
  return sum((d.-d_mean)./d_std .> k)
end

"""
Transform Nx3 pts by 3x3 tform matrix
"""
function warp_pts(tform, pts)
    pts = hcat(pts, ones(size(pts,1)))
    tpts = pts * tform
    return tpts[:,1:2]
end

"""
Make N array of 1x2 points into Nx3 homogenous points
"""
function points_to_Nx3_matrix(points)
  points = hcat(points...)
  if size(points,1)==2
    points = [points; ones(eltype(points), 1, size(points,2))]'
  end
  return points
end

"""
Edit the offset_log text file associated with prealignment and alignment

log_path: file path of the log file, a .txt file
image_name: string not including the file extension
offset: 2-element collection for the i,j offset
sz: 2-element collection for the i,j height and width
"""
function update_offset_log!(log_path, image_name, offset, sz)
  if !isfile(log_path)
    f = open(log_path, "w")
    close(f)
    offset_log = [image_name, offset..., sz...]'
  else  
    offset_log = readdlm(log_path)
    idx = findfirst(offset_log[:,1], image_name)
    if idx != 0
      offset_log[idx, 2:3] = collect(offset)
      offset_log[idx, 4:5] = collect(sz)
    else
      log_line = [image_name, offset..., sz...]
      offset_log = vcat(offset_log, log_line')
    end
  end
  writedlm(log_path, offset_log)
end





