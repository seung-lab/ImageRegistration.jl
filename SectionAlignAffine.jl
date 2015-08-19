include("convolve.jl")

function test()
	getimage(path) = convert(Array{Float64, 2}, data(imread(path)));
	sec1 = getimage("./sections/S2-W001_sec24_0.175.tif")
	#sec2 = getimage("./sections/S2-W001_sec24_0.175.tif")
	sec2 = getimage("./sections/S2-W001_sec25_0.175.tif")
	AffineAlignSections(sec1, sec2)
end


function AffineAlignSections(img1::Array{}, img2::Array{})
	# Assuming same size
	border_ratio = 0.1;
	border = round(Int, border_ratio * minimum(size(img1)));
	search_radius = 200;
	half_block_size = 40;

	border = border > search_radius + half_block_size ? border : search_radius + half_block_size;
	println(border)

	accept_xcorr = 0.3;

	# four corners for now
	points = [[x; y] for x = (1+border,size(img1,1)-border), y = (1+border,size(img1,2)-border)]
	points = points[:];
	println(size(points))
	points1 = Array{Int,1}[]
	points2 = Array{Int,1}[]
	for p = points
		offset,r = BlockMatchAtPoint(img1, p, img2, p, half_block_size, search_radius)
		println(p, offset, r)
		if r > accept_xcorr
			push!(points1, p)
			push!(points2, p + offset)
		end
	end
	points1 = hcat(points1[:]...)
	points2 = hcat(points2[:]...)
	return FindAffine(points1, points2)
end


function FindAffine(points1, points2)
# Inputs:
# 		2D arrays with each column being a point.
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

	dim1 = size(xc, 1);
	(i_max, j_max) = (rem(ind, dim1), cld(ind, dim1));
	if i_max == 0 i_max = dim1; end

	return [i_max - 1 - search_r; j_max - 1 - search_r], r_max, xc;
	
end
