include("convolve.jl")


using Images
include("imwarp.jl")
include("visualize.jl")


function test()
	getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)));
	sec1 = getimage("./sections/S2-W001_sec24_0.175.tif")
	#sec2 = getimage("./sections/S2-W001_sec24_0.175.tif")
	sec2 = getimage("./sections/S2-W001_sec25_0.175.tif")
	println(size(sec1))
	println(size(sec2))
	#trans, points1, points2 = AffineAlignSections(sec1, sec2, 0.175)
	points, half_block_size, search_radius = GenerateMatchPoints(sec1, sec2)
	points1, points2 = GetBlockMatches(sec1, sec2, points, half_block_size, search_radius, 0.3)
	trans = FindAffine(points1, points2)
	A = AdjustAffineForScaling(trans, 0.175)

	println(A)
	#out_img, offset = imwarp(sec2, trans)
	#imwrite(sec1, joinpath(".","test_outputs", string("sec1", ".tif")));
	#imwrite(sec2, joinpath(".","test_outputs", string("sec2_", offset[1], "_", offset[2], ".tif")));
	#imwrite(out_img, joinpath(".","test_outputs", string("warped_", offset[1], "_", offset[2], ".tif")));
	points1 = [points1[2,:]; points1[1,:]]
	points2 = [points2[2,:]; points2[1,:]]
	draw_points(sec1, points1)
	draw_points(sec2, points2)
	#draw_points(out_img, points1)
end


function AffineAlignSections(img1::Array{}, img2::Array{}, downsample_ratio)
	points, half_block_size, search_radius = GenerateMatchPoints(img1, img2)
	points1, points2 = GetBlockMatches(img1, img2, points, half_block_size, search_radius, 0.3)
	A = FindAffine(points1, points2)
	A = AdjustAffineForScaling(A, downsample_ratio)
end


function GenerateMatchPoints(img1::Array{}, img2::Array{})
	border_ratio = 0.1;
	radius_ratio = 0.05;
	grid_size = 3;
	half_block_size = 150;
	overlap = [min(size(img1), size(img2))...]
	#overlap = size(img1)

	border = round(Int, border_ratio * minimum(overlap));
	search_radius = round(Int, minimum(overlap) * radius_ratio);
	search_radius = search_radius > 200 ? search_radius : 200;
	println("search_radius: ", search_radius)

	border = border > search_radius + half_block_size ? border : search_radius + half_block_size;
	println(border)

	# four corners for now
	#points = [[x; y] for x = (1+border,size(img1,1)-border), y = (1+border,size(img1,2)-border)]
	#points = [[x; y] for x = (1+border,overlap[1]-border), y = (1+border,overlap[2]-border)]

	step = floor(Int, (overlap - 2*border - 1) ./ (grid_size-1));
	points = [[x; y] for x = 1+border:step[1]:overlap[1]-border, y = 1+border:step[2]:overlap[2]-border]
	points = points[:];
	println("n points: ", size(points))

	return points, half_block_size, search_radius
end

function GetBlockMatches(img1::Array{}, img2::Array{}, points, half_block_size, search_radius, accept_xcorr = 0.3)
	accept_xcorr = 0.;
	points1 = Array{Int,1}[]
	points2 = Array{Int,1}[]
	for p = points
		offset,r = BlockMatchAtPoint(img1, p, img2, p, half_block_size, search_radius)
		println(p, offset, r)
		if r > accept_xcorr
			push!(points1, [p; 1])
			push!(points2, [p + offset; 1])
		end
	end
	points1 = hcat(points1[:]...)
	points2 = hcat(points2[:]...)
	return points1, points2
end


function AdjustAffineForScaling(A, scale)
	# Column vector convention. I.e. p_transformed = A * p.
	ss = 1. / scale;
	B = copy(A)
	# equivalent of inv(S) * A * S, where S = [1 0 0; 0 1 0; 0 0 1/scale]
	B[1:end-1, end] *= ss
	return B
end


function FindAffine(points1, points2)
# Inputs:
# 		2D arrays with each column being a point in homogeneous coords.
# Returns:
#    	Affine transformation A, that produces points2 = A * points1

	A = points2 / points1
end

function BlockMatchAtPoint(A, pointInA, B, pointInB, half_block_size, search_r)
# A is the template
# Returns:
#	offset: offset from pointInB to pointInA's actual match in B.
#	r_max:  correlation value at the match point
#	xc:		raw cross-correlation map

	b = half_block_size
	B_radius = half_block_size + search_r;

	A_lower = pointInA - b
	A_upper = pointInA + b
	B_lower = pointInB - B_radius
	B_upper = pointInB + B_radius

	xc = normxcorr2(sub(A, [l:u for (l,u) in zip(A_lower, A_upper)]...),
					sub(B, [l:u for (l,u) in zip(B_lower, B_upper)]...));

	r_max, ind = findmax(xc); 
	i_max, j_max = ind2sub(xc, ind);

	return [i_max - 1 - search_r; j_max - 1 - search_r], r_max, xc;
	
end
