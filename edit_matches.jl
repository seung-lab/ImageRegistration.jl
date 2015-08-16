using Images
using ImageView

include("render.jl")

function find_idx_of_closest_pt(pts, pt, limit)
	d = sum((pts.-pt).^2, 1).^(1/2)
	idx = eachindex(d)'[d .< limit]
	if length(idx) == 1
		return idx[1]
	else
		return 0
	end
end

function edit_matches(imgc, img2, annotation)
	e = Condition()

	pts_to_remove = []
	pts = copy(annotation.ann.data.pts)

	c = canvas(imgc)
	win = Tk.toplevel(c)
	bind(c, "<Button-3>", (c, x, y)->right_click(parse(Int, x), parse(Int, y)))
	bind(win, "<Return>", path->end_edit())
	bind(win, "<KP_Enter>", path->end_edit())
	bind(win, "<Destroy>", path->end_edit())

	function right_click(x, y)
		xu,yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
		xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
		# println(xi, ", ", yi)
		limit = (img2.zoombb.xmax - img2.zoombb.xmin) * 0.0125 # 100/8000
		annpts = annotation.ann.data.pts
		annidx = find_idx_of_closest_pt(annpts, [xi, yi], limit)
		if annidx > 0
			idx = find_idx_of_closest_pt(pts, [xi, yi], limit)
			println(idx, ": ", [xi, yi])
			annotation.ann.data.pts = hcat(annpts[:,1:annidx-1], annpts[:,annidx+1:end])
			ImageView.redraw(imgc)
			push!(pts_to_remove, idx)
		end
	end

	function end_edit()
		println("End edit")
		notify(e)
		bind(c, "<Button-3>", path->path)
		bind(win, "<Return>", path->path)
		bind(win, "<KP_Enter>", path->path)
		bind(win, "<Destroy>", path->path)
	end

	println("Right click to remove correspondences, then press enter.")
	wait(e)

	return pts_to_remove
end

function demo()
	mesh_set = load(joinpath(BUCKET, "EM_images", "Test.jld"))["MeshSet"]
	src_pts, dst_pts = load_matches(mesh_set.matches[1])
	tile_pathB = joinpath(BUCKET, "EM_images", "Tile_r4-c3_S2-W001_sec20.tif")
	# draw_points(make_isotropic(rawdata(imread(tile_pathA))), src_pts)
	imgc, img2, annotation = draw_points(make_isotropic(rawdata(imread(tile_pathB))), dst_pts)
	a = edit_matches(imgc, img2, annotation)
end
