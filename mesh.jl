using Images
using HDF5
using JLD
include("convolve.jl")

type Tile
	index::Tuple{Integer, Integer}
	array::Array{FloatingPoint, 2}
	dx::FloatingPoint
	dy::FloatingPoint
end

# EVENTUALLY ALL MESHES SHOULD EXTEND SOME TYPE - SAY 'MESHABLE' - WITH FUNCTIONS THAT HAVE CORRESPONDING METHODS
type Mesh
	index::Tuple{Integer, Integer}	# sectionIndex, tileIndex as a tuple - if tileIndex is 0 then denotes entire section / collection of multiple meshes

	di::FloatingPoint # displacement of the tile in the i direction
	dj::FloatingPoint # displacement of the tile in the j direction

	dims::Tuple{Integer, Integer} # mesh dimensions in terms of nodes in the i direction, j direction. (0, 0) if the mesh is not a primitive mesh
	offsets::Tuple{FloatingPoint, FloatingPoint} # mesh offset in terms of the top left, from the image
	dists::Tuple{FloatingPoint, FloatingPoint} # mesh distance in each direction

	n::Integer # number of nodes in the mesh
	m::Integer # number of edges in the mesh

	nodes::Array{FloatingPoint, 2} # 2-by-n dense matrix of nodes, each column stores [i-coordinate, j-coordinate]
	nodes_t::Array{FloatingPoint, 2} # 2-by-n dense matrix of nodes after transformation, in the same order as nodes - by default same as nodes
	nodes_fixed::Array{Bool, 1} # 1-by-n dense matrix of nodes that denote whether points are fixed

	edges::SparseMatrixCSC{FloatingPoint, Integer} # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
	                            		       # with the k-th row denoting the node stored in the k-th column of nodes - there will not be any columns with all zeros
	edge_lengths::Array{FloatingPoint, 1}  # 1-by-m dense matrix of floats that stores the rest lengths.
	edge_coeffs::Array{FloatingPoint, 1}   # 1-by-m dense matrix of floats that stores the spring coefficients.
end

type Matches
	src_index::Tuple{Integer, Integer} 	# source index, sectionIndex, tileIndex
	dst_index::Tuple{Integer, Integer}	# destination index, sectionIndex, tileIndex

	n::Integer				# number of points in the mesh

	src_pointIndices::Array{Integer, 1}	# keeps track of which points have been matched
	src_points::Array{FloatingPoint, 2} 	# 2-by-n matrix of points at the source, each column stores [i-coordinate, j-coordinate]
	dst_points::Array{FloatingPoint, 2} 	# 2-by-n dense matrix of points at the destination, each column stores [i-coordinate, j-coordinate]

	dispVectors::Array{FloatingPoint, 2}
end

function Mesh2JLD(filename, mesh::Mesh)
	fid = jldopen(filename, "w");
	fid["index"] = mesh.index;
	fid["di"] = mesh.di;
	fid["dj"] = mesh.dj;
	fid["n"] = mesh.n;
	fid["m"] = mesh.m;
	fid["nodes"] = mesh.nodes;
	fid["nodes_t"] = mesh.nodes_t;
	fid["nodes_fixed"] = mesh.nodes_fixed;
	fid["edges"] = mesh.edges;
	fid["edge_lengths"] = mesh.edge_lengths;
	fid["edge_coeffs"] = mesh.edge_coeffs;
	close(fid);
end


# makeTileMesh(A, di, dj, mesh_length, mesh_coeff)
# returns a TileMesh for a matrix A positioned at (di, dj) with the specified mesh length and mesh spring coefficient.
function makeTileMesh(A, di, dj, mesh_length, mesh_coeff)

	(Ai, Aj) = size(A);
	(Mi, Mj) = (mesh_length * sin(pi / 3), mesh_length);

	mA = (convert(Integer, div(Ai, Mi) + 1), convert(Integer, div(Aj, Mj) + 1));
 	rA = (rem(Ai, Mi)/2, rem(Aj, Mj)/2);

	n = getMeshIndex(mA, mA[1], mA[2]);
	if n == 0
		n = getMeshIndex(mA, mA[1], mA[2]-1);
	end

	m_upperbound = 3 * n;

	nodes = Array(FloatingPoint, 2, n);
	nodes_fixed = fill(false, n);
	edges = spzeros(FloatingPoint, n, m_upperbound);
	edge_lengths = fill(convert(FloatingPoint, mesh_length), m_upperbound);
	edge_coeffs = fill(convert(FloatingPoint, mesh_coeff), m_upperbound);

	m = 0;	

	for i in 1:mA[1], j in 1:mA[2]
		k = getMeshIndex(mA, i, j);
		if k == 0
			continue;
		end

		(pi, pj) = getMeshCoord(mA, (rA[1] + di, rA[2] + dj), Mi, Mj, i, j);
		nodes[:, k] = [pi, pj]

		if (j != 1)
			m += 1;
			edges[k, m] = -1;
			edges[getMeshIndex(mA, i, j-1), m] = 1;
		end

		if (i != 1)
			if !((i % 2 == 1) && (j == mA[2]))
				m += 1;
				edges[k, m] = -1;
				edges[getMeshIndex(mA, i-1, j), m] = 1;
			end

			if (i % 2 == 0) && (j != mA[2]) 			
				m += 1;
				edges[k, m] = -1;
				edges[getMeshIndex(mA, i-1, j+1), m] = 1;
			elseif (i % 2 == 1) && (j != 1)
				m += 1;
				edges[k, m] = -1;
				edges[getMeshIndex(mA, i-1, j-1), m] = 1;
			end
			
			if (i % 2 == 1) && ((j == 1) || (j == mA[2]))
				m += 1;
				edges[k, m] = -1;
				edges[getMeshIndex(mA, i-2, j), m] = 1;
				edge_lengths[m] = 2 * Mi;
			end
		end
	end

	edges = edges[:, 1:m];
	edge_lengths = edge_lengths[1:m];
	edge_coeffs = edge_coeffs[1:m];

	mesh = Mesh((0, 0), di, dj, mA, rA, (Mi, Mj), n, m, nodes, nodes, nodes_fixed, edges, edge_lengths, edge_coeffs);
	return mesh;
end

function getMeshIndex(mA, i, j)
	if (i % 2 == 0) && (j == mA[2]) 
		return 0;
	end

	if (i < 1) || (j < 1) || (i > mA[1]) || (j > mA[2])
		return 0;
	end
	
	ind = 0

	ind += div(i-1, 2) * (mA[2] - 1); #even rows
	ind += div(i, 2) * mA[2]; #odd rows
	ind += j;

	ind = convert(Integer, ind);

	return ind;
end

function getMeshCoord(mA, rA, Mi, Mj, i, j)
	if (i % 2 == 0) && (j == mA[2])
	return (0, 0);
	end;

	pi = (i-1) * Mi + rA[1];
	if i % 2 == 0
		pj = (j-0.5) * Mj + rA[2];
	else
		pj = (j-1) * Mj + rA[2];
	end

	return (pi, pj);
end

function isInternal(A, i, j, d)
	if i > d && i <= (size(A, 1) - d) && j > d && j <= (size(A, 2) - d)
		return true
	end
	return false
end

# i, j are in world coordinates (as per the mesh coordinate specification)
function getMaxVector_normxcorr2(A, Ai_g, Aj_g, dAi, dAj, B, dBi, dBj, block_size, search_r)
	b_rad = block_size + search_r;
	noMatch = [0 0 -1];

	# convert to matrix coordinates
	Ai = round(Integer, Ai_g - dAi);
	Aj = round(Integer, Aj_g - dAj);

	Bi = round(Integer, Ai_g - dBi);
	Bj = round(Integer, Aj_g - dBj);

	if (!isInternal(A, Ai, Aj, block_size) || !isInternal(B, Bi, Bj, b_rad))
		return noMatch;
	end

	#	println("($Ai, $Aj)->($Bi, $Bj): ");

	Ai_range = ceil(Ai)-block_size:ceil(Ai)+block_size;
	Aj_range = ceil(Aj)-block_size:ceil(Aj)+block_size;
	Bi_range = ceil(Bi)-b_rad:ceil(Bi)+b_rad;
	Bj_range = ceil(Bj)-b_rad:ceil(Bj)+b_rad;

	xc = normxcorr2(A[Ai_range, Aj_range], B[Bi_range, Bj_range]);
	r_max = maximum(xc);
	ind = find(xc .== r_max);
	if (size(ind) == 0)
		return noMatch;
	end
	(i_max, j_max) = (rem(ind[1], size(xc, 2)), cld(ind[1], size(xc, 2)));

	return [i_max - 1 - search_r; j_max - 1 - search_r; r_max];
end

function extractMatchesAtMesh(A, A_Mesh, B, dBi, dBj, block_size, search_r, min_r)
	src_index = A_Mesh.index;
	dst_index = (A_Mesh.index[1], A_Mesh.index[2] + 1);

	n = 0;
	n_upperbound = convert(Integer, A_Mesh.n);

	src_points = Array(FloatingPoint, 2, n_upperbound);
	src_pointIndices = Array(Integer, n_upperbound);
	dst_points = Array(FloatingPoint, 2, n_upperbound);
	dispVectors = Array(FloatingPoint, 2, n_upperbound);

	for j in 1:n_upperbound
		Ai = A_Mesh.nodes[1, j];
		Aj = A_Mesh.nodes[2, j];
	#	print("Matching mesh point $j of $n_upperbound: ");
		v = getMaxVector_normxcorr2(A, Ai, Aj, A_Mesh.di, A_Mesh.dj, B, dBi, dBj, block_size, search_r);
		if v[3] < min_r
	#		println("No match found.");
			continue
		end
	#	println("$v");
		n += 1;
		src_pointIndices[n] = j;
		dispVectors[:, n] = v[1:2];
		src_points[:, n] = A_Mesh.nodes[:, j];
		dst_points[:, n] = src_points[:, n] + dispVectors[:, n];
	end

	src_pointIndices = src_pointIndices[1:n];
	src_points = src_points[:, 1:n];
	dst_points = dst_points[:, 1:n];
	dispVectors = dispVectors[:, 1:n];

	println("$n total matches found of $n_upperbound mesh points.");

	matches = Matches(src_index, dst_index, n, src_pointIndices, src_points, dst_points, dispVectors);
	return matches;
end

# add fixed nodes to the input mesh with the corresponding matches - nodes_t
function addFixedMatches!(A_Mesh, matches, spring_coeff)

	n_new = A_Mesh.n + matches.n;
	m_new = A_Mesh.m + matches.n;

	nodes_add = matches.dst_points;
	nodes_fixed_add = fill(true, matches.n); 

	edges_new = similar(A_Mesh.edges, n_new, m_new);
	edges_new[1:A_Mesh.n, 1:A_Mesh.m] = A_Mesh.edges[:, :]

	edge_lengths_add = fill(convert(FloatingPoint, 0), matches.n);
	edge_coeffs_add = fill(convert(FloatingPoint, spring_coeff), matches.n);

	for j in [1:matches.n]
		edges_new[matches.src_pointIndices[j], A_Mesh.m + j] = -1;
		edges_new[A_Mesh.n + j, A_Mesh.m + j] = 1;
	end

	A_Mesh.n = n_new;
	A_Mesh.m = m_new;
	A_Mesh.nodes = hcat(A_Mesh.nodes, nodes_add);
	A_Mesh.nodes_t = hcat(A_Mesh.nodes_t, nodes_add);
	A_Mesh.nodes_fixed = vcat(A_Mesh.nodes_fixed, nodes_fixed_add);
	A_Mesh.edges = edges_new;
	A_Mesh.edge_lengths = vcat(A_Mesh.edge_lengths, edge_lengths_add);
	A_Mesh.edge_coeffs = vcat(A_Mesh.edge_coeffs, edge_coeffs_add);
end

# find the triangular mesh indices for a given point in A
function findMeshTriangle(Am, i, j)

	noTriangle = (0, 0, 0);

	# convert to A's local coordinates, and displace by the mesh offset
	Ai = i - Am.di - Am.offsets[1];
	Aj = j - Am.dj - Am.offsets[2];

	# find which rows the point belongs to
	i0 = round(Integer, Ai / Am.dists[1] + 1);

	if isodd(i0)
		j0 = round(Integer, Aj / Am.dists[2] + 1);
	else
		j0 = round(Integer, Aj / Am.dists[2] + 0.5);
	end

	ind0 = getMeshIndex(Am.dims, i0, j0);
	
	if ind0 == 0
		return noTriangle;
	end

	node0 = Am.nodes[:, ind0];

	#println("$i0, $j0 at ($node0)");

	di = i - node0[1];
	dj = j - node0[2];

	#println("$di, $dj");

	theta = abs(atan(di / dj));
	if (theta < pi / 3)
		if (dj >= 0)
			ind1 = getMeshIndex(Am.dims, i0, j0 + 1);
			if (di >= 0)
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0+1, j0);
				else
					ind2 = getMeshIndex(Am.dims, i0+1, j0+1);
				end
			else
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0-1, j0);
				else
					ind2 = getMeshIndex(Am.dims, i0-1, j0+1);
				end
			end
		else
			ind1 = getMeshIndex(Am.dims, i0, j0 - 1);
			if (di >= 0)
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0+1, j0-1);
				else
					ind2 = getMeshIndex(Am.dims, i0+1, j0);
				end
			else
				if isodd(i0)
					ind2 = getMeshIndex(Am.dims, i0-1, j0-1);
				else
					ind2 = getMeshIndex(Am.dims, i0-1, j0);
				end
			end
		end
	else
		if (di >= 0)
			if isodd(i0)
				ind1 = getMeshIndex(Am.dims, i0+1, j0-1);
				ind2 = getMeshIndex(Am.dims, i0+1, j0);
			else
				ind1 = getMeshIndex(Am.dims, i0+1, j0);
				ind2 = getMeshIndex(Am.dims, i0+1, j0+1);
			end
		else
			if isodd(i0)
				ind1 = getMeshIndex(Am.dims, i0-1, j0-1);
				ind2 = getMeshIndex(Am.dims, i0-1, j0);
			else
				ind1 = getMeshIndex(Am.dims, i0-1, j0);
				ind2 = getMeshIndex(Am.dims, i0-1, j0+1);
			end
		end
	end

	if (ind1 == 0 || ind2 == 0)
		return noTriangle;
	end
	return (ind0, ind1, ind2);
end

# Convert Cartesian coordinate to triple of barycentric coefficients
function triangleWeights(Am, triangle, pi, pj)
	R = vcat(Am.nodes[:, triangle[1]]', Am.nodes[:, triangle[2]]', Am.nodes[:, triangle[3]]')
	R = hcat(R, ones(FloatingPoint, 3, 1));
	r = hcat(pi, pj, 1.0);
	return r * R^-1;
end

function makeMatchedMesh(A, Am, B, Bm, block_size, search_r, min_r, match_coeff)
	A2B = extractMatchesAtMesh(A, Am, B, Bm.di, Bm.dj, block_size, search_r, min_r);
	B2A = extractMatchesAtMesh(B, Bm, A, Am.di, Am.dj, block_size, search_r, min_r);
	n = Am.n + Bm.n;

	Ae_padded = vcat(Am.edges, spzeros(Bm.n, Am.m));
	Be_padded = vcat(spzeros(Am.n, Bm.m), Bm.edges);
	
	edges_new = spzeros(n, A2B.n + B2A.n); 

	j = 0;
	for i in 1:A2B.n
		pi = A2B.dst_points[1, i]; 
		pj = A2B.dst_points[2, i];
		triangle = findMeshTriangle(Bm, pi, pj);
		#print(triangle);
		if triangle == (0, 0, 0)
			continue;
		end;
		j += 1;
		weights = triangleWeights(Bm, triangle, pi, pj);
		edges_new[A2B.src_pointIndices[i], j] = -1;
		for k in 1:3
			edges_new[Am.n + triangle[k], j] = weights[k];
		end;
	end
	for i in 1:B2A.n
		pi = B2A.dst_points[1, i]; 
		pj = B2A.dst_points[2, i];
		triangle = findMeshTriangle(Am, pi, pj);
		if triangle == (0, 0, 0)
			continue;
		end
		j += 1;
		weights = triangleWeights(Am, triangle, pi, pj);
		edges_new[Am.n + B2A.src_pointIndices[i], j] = -1;
		for k in 1:3
			edges_new[triangle[k], j] = weights[k];
		end;
	end

	m = Am.m + Bm.m + j;

	nodes = hcat(Am.nodes, Bm.nodes);
	nodes_t = hcat(Am.nodes_t, Bm.nodes_t);
	nodes_fixed = vcat(Am.nodes_fixed, Bm.nodes_fixed);

	edges = hcat(Ae_padded, Be_padded, edges_new[:, 1:j]);
	edge_lengths = vcat(Am.edge_lengths, Bm.edge_lengths, zeros(j));
	edge_coeffs = vcat(Am.edge_coeffs, Bm.edge_coeffs, fill(convert(FloatingPoint, match_coeff), j));

	mesh = Mesh((0, 0), 0, 0, (0, 0), (0,0), (0, 0), n, m, nodes, nodes, nodes_fixed, edges, edge_lengths, edge_coeffs);
	return mesh;
end
################################# SCRIPT FOR TESTING ###################################
tic();
A = convert(Array{Float64, 2}, data(imread("./Test/Tile_r4-c2_S2-W001_sec20.tif")));
toc();
dAi = 21906;
dAj = 36429;

tic();
B = convert(Array{Float64, 2}, data(imread("./Test/Tile_r4-c3_S2-W001_sec20.tif")));
toc();
dBi = 29090; # 2908.6;
dBj = 36251; # 3624.3;

block_size = 20;
search_r = 80;
min_r = 0.5;
mesh_length = 100;
mesh_coeff = 0.5;
match_coeff = 1.0;
#=

A = rand(11, 11);
dAi = 0;
dAj = 0;
mesh_length = 2.0;
mesh_coeff = 1.0;
=#
@time Am = makeTileMesh(A, dAi, dAj, mesh_length, mesh_coeff);
@time Bm = makeTileMesh(B, dBi, dBj, mesh_length, mesh_coeff);
@time C = makeMatchedMesh(A, Am, B, Bm, block_size, search_r, min_r, match_coeff);

@time Mesh2JLD("r4c2-r4c3.jld", C);
#=
@time addFixedMatches!(Am, matches, match_coeff); 
@time Mesh2JLD("r4c2-r4c3.jld", Am);
=#

#=

A = rand(11, 11);
A_Mesh = makeTileMesh(A, 0, 0, 2, 0.5);


matches = extractMatchesAtMesh(A, A_Mesh, dAi, dAj, B, dBi, dBj, block_size, search_r, min_r);
=#
