#=
include("convolve.jl")


using Images
include("imwarp.jl")
include("visualize.jl")
#include("render.jl")	# cyclic inclusion
=#

function test()
	getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)));
	sec1 = getimage("../sections/S2-W001_sec24_0.175.tif")
	#sec2 = getimage("./sections/S2-W001_sec24_0.175.tif")
	sec2 = getimage("../sections/S2-W001_sec25_0.175.tif")
	println(size(sec1))
	println(size(sec2))
	trans, points1, points2, res1, res2 = AffineAlignSections(sec1, sec2, PARAMS_PREALIGNMENT; return_points=true)
	points1 = points_to_3xN_matrix(points1)
	points2 = points_to_3xN_matrix(points2)

	println(trans)
	#out_img, offset = imwarp(sec2, trans)
	#imwrite(sec1, joinpath(".","test_outputs", string("sec1", ".tif")));
	#imwrite(sec2, joinpath(".","test_outputs", string("sec2_", offset[1], "_", offset[2], ".tif")));
	#imwrite(out_img, joinpath(".","test_outputs", string("warped_", offset[1], "_", offset[2], ".tif")));
	p22 = points2 + res2;
	p11 = points1 + res1;
	points1 = points1[1:2,:]
	points2 = points2[1:2,:]
	draw_vectors(sec2, vcat(points2[1:2,:], p22[1:2,:]))
	#ccp.write_image_from_points(points1[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
	draw_points(sec1, points1)
	draw_points(sec2, points2)
	#draw_points(out_img, points1)
end


function test2()
	getimage(path) = convert(Array{Float64, 2}, convert(Array, imread(path)));
	#getimage(path) = convert(Array{Ufixed8, 2}, convert(Array, imread(path)));
	#getimage(path) = reinterpret(UInt8, convert(Array, imread(path)));
	sec1 = getimage("../output_images/(1,1)_montage.tif")
	sec2 = getimage("../output_images/(1,2)_montage.tif")
	println(size(sec1))
	println(size(sec2))
	trans, points1, points2, res1, res2 = AffineAlignSections(sec1, sec2, PARAMS_PREALIGNMENT; return_points=true)
	points1 = points_to_3xN_matrix(points1)
	points2 = points_to_3xN_matrix(points2)

	println(trans)
	downsample = 4
	sec1 = sec1[1:downsample:end, 1:downsample:end];
	sec2 = sec2[1:downsample:end, 1:downsample:end];

	p11 = points1 + res1;
	p22 = points2 + res2;

	p11 = p11[1:2,:]
	p22 = p22[1:2,:]
	points1 = points1[1:2,:]
	points2 = points2[1:2,:]

	p11 = p11/downsample;
	p22 = p22/downsample;
	points1 = points1/downsample;
	points2 = points2/downsample;
	
	trans = AdjustAffineForScaling(trans.', downsample).'
	out_img, offset = imwarp(sec2, inv(trans))
	println(offset)
	fused, fused_offset = imfuse(sec1, [0,0], out_img, offset)
	println(fused_offset)

	block_radius = 150;
	scalebar = [1; 1; 2*block_radius/downsample; 2*block_radius/downsample]
	draw_vectors(make_isotropic(sec1), hcat(vcat(points1, p11), scalebar))
	draw_vectors(make_isotropic(sec2), hcat(vcat(points2, p22), scalebar+100))
	draw_vectors(make_isotropic(fused), hcat(vcat(points1.-fused_offset, p11.-fused_offset), scalebar)) # todo: check offset
    #ccp.write_image_from_points(points1[1:2,:].', p11[2:-1:1,:].', "test_outputs/write_name.jpg")
end

function test3()
	sec1 = "../output_images/(1,1)_montage.tif"
	sec2 = "../output_images/(1,2)_montage.tif"
	A, meshset = AffineAlignSections(sec1, sec2)
	return meshset
end

function points_to_3xN_matrix(points)
	points = hcat(points...)
	if size(points,1)==2
		points = [points; ones(eltype(points), 1, size(points,2))]
	end
	return points
end

"""
`AffineAlignSections`

Returns the affine transformation A, such that p_in_B = p_in_A * A, where p is point coordinates in row vector.
"""
function AffineAlignSections(img1::Array{}, img2::Array{}, params=PARAMS_PREALIGNMENT; return_points=false)
	# points here are in column vector convention
	points = GenerateMatchPoints(img1, img2, params)
	points1list, points2list = GetBlockMatches(img1, img2, points, params)

	n_matches = length(points1list)
	println("n_matches: ", n_matches)

	points1 = points_to_3xN_matrix(points1list)
	points2 = points_to_3xN_matrix(points2list)

	A = FindAffine(points1, points2)
	residualIn2 = A*points1 - points2;
	rmsIn2 = mean( sum(residualIn2.^2, 1) )^0.5
	points2in1 = inv(A)*points2;
	residualIn1 = points2in1 - points1;
	rmsIn1 = mean( sum(residualIn1.^2, 1) )^0.5

	tolerance_ratio = 0.001
	rmsThres = tolerance_ratio * mean([size(img1)..., size(img2)...])
	rmsTotal = mean([rmsIn1, rmsIn2].^2)^0.5

	if n_matches < 0.5 * length(points)
		println("WARNING [AffineAlignSections]: # of matches is small. n_matches: ", n_matches)
	end
	if  rmsTotal > rmsThres
		println("WARNING [AffineAlignSections]: high residual. RMS error: ", rmsTotal)
	end

	A = AdjustAffineForScaling(A, params.scaling_factor)

	# convert column vector convention to row vector convention
	A = A.'

	if return_points
		return A, points1list, points2list, residualIn1, residualIn2, rmsIn1, rmsIn2, rmsTotal
	else
		return A
	end
end


function GenerateMatchPoints(img1::Array{}, img2::Array{}, params=PARAMS_PREALIGNMENT)
	border_ratio = 0.1;
	grid_size = minimum(floor(Int, collect(size(img1))/params.mesh_length))
	block_radius = params.block_size
	search_radius = params.search_r
	overlap = [min(size(img1), size(img2))...]
	#overlap = size(img1)

	border = round(Int, border_ratio * minimum(overlap))
	border = border > search_radius + block_radius ? border : search_radius + block_radius
	println("border excluded: ", border)

	# four corners for now
	#points = [[x; y] for x = (1+border,size(img1,1)-border), y = (1+border,size(img1,2)-border)]
	#points = [[x; y] for x = (1+border,overlap[1]-border), y = (1+border,overlap[2]-border)]

	step = floor(Int, (overlap - 2*border - 1) ./ (grid_size-1));
	points = [[x; y] for x = 1+border:step[1]:overlap[1]-border, y = 1+border:step[2]:overlap[2]-border]
	points = points[:];
	println("n points: ", size(points))
	return points
end

function GetBlockMatches(img1::Array{}, img2::Array{}, points, params=PARAMS_PREALIGNMENT)
	points1 = []
	points2 = []
	for p = points
		offset,r = BlockMatchAtPoint(img1, p, img2, p, params.block_size, params.search_r)
		println(p, offset, r)
		if r >= params.min_r
			push!(points1, collect(p))
			push!(points2, collect(p + offset))
		end
	end
	return points1, points2
end


function AdjustAffineForScaling(A, scale)
	# Column vector convention. I.e. p_transformed = A * p.
	B = copy(A)
	# equivalent of inv(S) * A * S, where S = [1 0 0; 0 1 0; 0 0 1/scale]
	B[1:end-1, end] *= scale
	return B
end


function FindAffine(points1, points2)
# Inputs:
# 		2D arrays with each column being a point in homogeneous coords.
# Returns:
#    	Affine transformation A, that produces points2 = A * points1

	A = points2 / points1
end

function BlockMatchAtPoint(A, pointInA, B, pointInB, block_radius, search_r)
# A is the template
# Returns:
#	offset: offset from pointInB to pointInA's actual match in B.
#	r_max:  correlation value at the match point
#	xc:		raw cross-correlation map

	b = block_radius
	B_radius = block_radius + search_r;

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

function BlockMatchSections
end

function recompute_affine(meshset::MeshSet)
	params = meshset.params
	points1list = meshset.meshes[1].nodes
	points2list = meshset.matches[1].dst_points
	points1 = points_to_3xN_matrix(points1list)
	points2 = points_to_3xN_matrix(points2list)
	A = FindAffine(points1, points2)
	A = AdjustAffineForScaling(A, params.scaling_factor)
	# convert column vector convention to row vector convention
	return A.'
end

function Meshes2SectionMatches!(meshset::MeshSet, params; return_points=false)
	if meshset.N != 2
		error("Invalid Arguments")
	end

	img1 = []
	img2 = []
	try
		img1 = getFloatImage(meshset.meshes[1])
		img2 = getFloatImage(meshset.meshes[2])
	catch
		# mostly for testing purpose where name is the file path
		img1 = getFloatImage(meshset.meshes[1].name)
		img2 = getFloatImage(meshset.meshes[2].name)
	end

	A, points1, points2, residualIn1, residualIn2, rmsIn1, rmsIn2, rmsTotal = 
		AffineAlignSections(img1, img2, params; return_points=true)

	# in the future this might not be a good idea where nodes could be the fine alignment nodes
	meshset.meshes[1].nodes = points1
	src_points_indices = 1:length(points1)
	dst_points = points2

	matches = Matches(meshset.meshes[1].index, meshset.meshes[2].index, length(points1),
		src_points_indices, dst_points, [], [], []);
	meshset.matches = [matches]
	return A, meshset
end


"""
`AffineAlignSections`

Arguments:

 * sec1: index or name, that can be used to construct a Mesh object
 * sec2: same as above

Returns:

 * Affine transform and the MeshSet containing the matches.

"""
function AffineAlignSections(sec1, sec2, params=PARAMS_PREALIGNMENT)
	if isa(sec1, Index) && isa(sec2, Index) && length(sec1)==2 && length(sec2)==2
		sec1 = sec1..., MONTAGED_INDEX, MONTAGED_INDEX
		sec2 = sec2..., MONTAGED_INDEX, MONTAGED_INDEX
	end
	meshset = makeNewMeshSet(params)
	meshset.N = 2
	mesh1 = Mesh(sec1)
	mesh2 = Mesh(sec2)
	meshset.meshes = [mesh1, mesh2]
	return Meshes2SectionMatches!(meshset, params)
end

function prealign_directory()
    img_filenames = filter(x -> x[end-2:end] == "tif", readdir(MONTAGED_DIR))
    for filename_A in img_filenames
        img_preceding = filter(x -> parseName(x)[2]-1 == parseName(filename_A)[2], img_filenames)
        if length(img_preceding) > 0
            filename_B = img_preceding[1]
            println("Prealigning ", filename_B[1:end-4], " to ", filename_A[1:end-4])
            @time tform, meshset = AffineAlignSections(filename_A, filename_B)
   	  		# src_index = meshset.meshes[1].index
			# dst_index = meshset.meshes[2].index
   			# meshset.meshes[1].index = (src_index[1:2]..., PREALIGNED_INDEX, PREALIGNED_INDEX)
			# meshset.meshes[2].index = (dst_index[1:2]..., PREALIGNED_INDEX, PREALIGNED_INDEX)
            println("Saving JLD for ", filename_B[1:end-4], " to ", filename_A[1:end-4])
            save(meshset)         
        end
    end
end

function compute_propogated_transform(index::Index)
	index = (index[1:2]..., 0, 0)
	filenames = filter(x -> x[end-2:end] == "jld", readdir(PREALIGNED_DIR))
	indices = [(parseName(x), x) for x in filenames]
	sort!(indices)
	if indices[1][1] == ones(Int, 4)
		error("Could not parse JLD filename to index: ", indices[1][2])
	end
	for k in 2:length(indices)
		if indices[k][1] > index
			break;
		end
		if !isAdjacent(indices[k-1][1], indices[k][1], true)
			error("Missing section between ", indices[k-1][1], " and ", indices[k][1])
		end
	end
	T = diagm(ones(3));
	for k in 2:length(indices)
		ind, fn = indices[k]
		if ind > index
			break;
		end
		meshset = load(joinpath(PREALIGNED_DIR, fn))["MeshSet"]
		A = recompute_affine(meshset)
		# row vector homogeneous point convention
		T *= A
	end
	return T
end

function prealignment()
    @time prealign_directory()
    @time render_prealignment_for_directory()
end
