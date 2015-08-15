using Images
using ImageView

include("render.jl")

mesh_set = load(joinpath(BUCKET, "EM_images", "Test.jld"))["MeshSet"]
src_pts, dst_pts = load_matches(mesh_set.matches[1])
# draw_points(make_isotropic(rawdata(imread(tile_pathA))), src_pts)
imgc, img2, h = draw_points(make_isotropic(rawdata(imread(tile_pathB))), dst_pts)

function right_click(imgc, img2, h, x, y)
	xu,yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
	xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
	# println(xi, ", ", yi)
	limit = (img2.zoombb.xmax - img2.zoombb.xmin) * 0.0125 # 100/8000
	pts = h.ann.data.pts
	idx = find_idx_of_closest_pt(pts, [xi, yi], limit)
	if idx > 0
		println(idx, ": ", pts[:, idx])
		h.ann.data.pts = hcat(pts[:,1:idx-1], pts[:,idx+1:end])
		ImageView.redraw(imgc)
	end
end

function find_idx_of_closest_pt(pts, pt, limit)
	d = sum((pts.-pt).^2, 1).^(1/2)
	idx = eachindex(d)'[d .< limit]
	if length(idx) == 1
		return idx[1]
	else
		return 0
	end
end

c = canvas(imgc)
bind(c, "<Button-3>", (c, x, y)->right_click(imgc, img2, h, parse(Int, x), parse(Int, y)))


