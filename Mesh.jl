module MeshModule

export Mesh, Matches, MeshSet
export Tile2Mesh, getMeshImage, Meshes2Matches, makeNewMeshSet, addMesh2MeshSet!, addMatches2MeshSet!, solveMeshSet!, JLD2Mesh, Mesh2JLD

using Images
using HDF5
using JLD
include("convolve.jl")
include("MeshSolve.jl")

####### global variables ########
global noMatch = [0; 0; -1];
global noTriangle = (0, 0, 0);

####### typealiases for sanity #######
typealias Index Tuple{Int64, Int64, Int64};			# (wafer, section, tile)

typealias Triangle Tuple{Int64, Int64, Int64};			# index of three points of the triangle for some point
typealias Triangles Array{Triangle, 1};				# index of three points of the triangle for some point

typealias Weight Tuple{Float64, Float64, Float64};		# weights for respective triangle
typealias Weights Array{Weight, 1};				# weights for respective triangle

typealias Pair Tuple{Int64, Int64};				# useful for abstraction

typealias Point Array{Float64, 1};				# [i; j]
typealias Points Array{Point, 1};				# array of points
typealias BinaryProperty Array{Bool, 1};			# array of bools

typealias Edges SparseMatrixCSC{Float64, Int64}			# sparse array for edges - columns represent edges and the rows represent the nodes
typealias FloatProperty Array{Float64, 1}

type Mesh
	path;							# path to the image file

	index::Index						# wafer, section, tile index of the mesh				tileIndex = 0 if the tile is a whole section
	grid::Pair						# row, column of the tile within the section				(0, 0) if the tile is a whole section
	disp::Point						# displacement of the tile within the section				(0, 0) if the tile starts at the top left corner

	dims::Pair				 		# mesh dimensions in terms of nodes in the i direction, j direction. 	(0, 0) if the mesh is not a regular mesh
	offsets::Point				                # mesh offset in terms of the top left, from the image. 		[0; 0] if the mesh is not a regular mesh
	dists::Point				                # mesh distance in each direction. 					[0; 0] if the mesh is not a regular mesh

	n::Int64 						# number of nodes in the mesh
	m::Int64	 					# number of edges in the mesh

	nodes::Points						# 2-by-n dense matrix of nodes, each column stores [i-coordinate, j-coordinate] in global coordinates
	nodes_t::Points			 			# 2-by-n dense matrix of nodes after transformation, in the same order as nodes - by default same as nodes
	nodes_fixed::BinaryProperty				# 1-by-n dense vector of nodes that denote whether points are fixed

	edges::Edges				                # n-by-m sparse incidence matrix, each column is zero except at two elements with value -1 and 1
	edge_lengths::FloatProperty	  			# 1-by-m dense vector of floats that stores the rest lengths.
	edge_coeffs::FloatProperty	   			# 1-by-m dense vector of floats that stores the spring coefficients.
end

type Matches
	src_mesh::Mesh						# source mesh
	dst_mesh::Mesh						# destination mesh

	n::Int64						# number of matches

	src_pointIndices::Array{Int64, 1}			# index of points in src_mesh.points
	dst_points::Points					# location of points in the destination
	dst_triangles::Triangles				# index of the triangles
	dst_weights::Weights					# weights at each

	dispVectors::Points					# displacement vector src->dest
end


type MeshSet
	N::Int64						# number of meshes in the set
	M::Int64						# number of matches in the set - (a -> b) and (b -> a) are distinct

	indices::Array{Index, 1}				# wafer, section, tile index as a tuple - if tileIndex is 0 then denotes entire section

	n::Int64						# number of nodes in the set across the whole set
	m::Int64						# number of edges in the set across the whole set
	m_i::Int64						# number of internal edges in the set across the whole set
	m_e::Int64						# number of edges between meshes
	
	meshes::Array{Mesh, 1}					# vector of meshes in the set
	nodes_indices::Array{Int64, 1}				# vector of number of nodes before the n-th mesh. To get index at i-th node of n-th mesh, nodes_indices[n] + i.

	matches::Array{Matches, 1}				# vector of matches in the set
	matches_pairs::Array{Pair, 1}				# vector of index (in meshes) - (a, b) means the match is between (meshes[a] -> meshes[b])
end
#=
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
=#

function getMeshImage(mesh::Mesh)
	return convert(Array{Float64, 2}, data(imread(mesh.path)));
end

# Tile2Mesh(path, grid, di, dj, mesh_length, mesh_coeff)
function Tile2Mesh(path, index, grid, di, dj, tile_fixed, mesh_length, mesh_coeff)

	A = convert(Array{Float64, 2}, data(imread(path)));

	(Ai, Aj) = size(A);

	dists = [mesh_length * sin(pi / 3); mesh_length];
	dims = (convert(Int64, div(Ai, dists[1]) + 1), convert(Int64, div(Aj, dists[2]) + 1));
 	offsets = [rem(Ai, dists[1])/2; rem(Aj, dists[2])/2];
	disp = [di; dj];

	n = maximum([getMeshIndex(dims, dims[1], dims[2]); getMeshIndex(dims, dims[1], dims[2]-1)]);

	m_upperbound = 3 * n;

	nodes = Points(n);
	nodes_fixed = BinaryProperty(n); nodes_fixed[:] = tile_fixed;
	edges = spzeros(Float64, n, m_upperbound);
	edge_lengths = FloatProperty(m_upperbound); edge_lengths[:] = convert(Float64, mesh_length);
	edge_coeffs = FloatProperty(m_upperbound); edge_coeffs[:] = convert(Float64, mesh_coeff);


	m = 0;	

	for i in 1:dims[1], j in 1:dims[2]
		k = getMeshIndex(dims, i, j);
		if k == 0 continue; end

		nodes[k] = getMeshCoord(dims, disp+offsets, dists, i, j);
		
		if (j != 1)
			m += 1;	edges[k, m] = -1; edges[getMeshIndex(dims, i, j-1), m] = 1;
		end

		if (i != 1)
			if iseven(i) || j != dims[2]
				m += 1;	edges[k, m] = -1; edges[getMeshIndex(dims, i-1, j), m] = 1;
			end
			if iseven(i) && (j != dims[2]) 			
				m += 1; edges[k, m] = -1; edges[getMeshIndex(dims, i-1, j+1), m] = 1;
			end
			if isodd(i) && (j != 1)
				m += 1; edges[k, m] = -1; edges[getMeshIndex(dims, i-1, j-1), m] = 1;
			end
			if isodd(i) && ((j == 1) || (j == dims[2]))
				m += 1; edges[k, m] = -1; edges[getMeshIndex(dims, i-2, j), m] = 1;
				edge_lengths[m] = 2 * dists[1];
			end
		end
	end

	edges = edges[:, 1:m];
	edge_lengths = edge_lengths[1:m];
	edge_coeffs = edge_coeffs[1:m];

	return Mesh(path, index, grid, disp, dims, offsets, dists, n, m, nodes, nodes, nodes_fixed, edges, edge_lengths, edge_coeffs);
end

function getMeshIndex(dims, i, j)
	ind = 0;
	
	if iseven(i) && (j == dims[2]) return ind; end
	if ((i < 1) || (j < 1) || (i > dims[1]) || (j > dims[2])) return ind; end
	
	ind += div(i-1, 2) * (dims[2] - 1); #even rows
	ind += div(i, 2) * dims[2]; #odd rows
	ind += j;

	ind = convert(Int64, ind);
	return ind;
end

function getMeshCoord(dims, total_offset, dists, i, j)
	if iseven(i) && (j == dims[2]) return (0, 0); end
	
	pi = (i-1) * dists[1] + total_offset[1];

	if iseven(i)	pj = (j-0.5) * dists[2] + total_offset[2];
	else		pj = (j-1) * dists[2] + total_offset[2];
	end

	return [pi; pj];
end

function isInternal(A, i, j, d)
	if i > d && i <= (size(A, 1) - d) && j > d && j <= (size(A, 2) - d)
		return true
	end
	return false
end

# i, j are in world coordinates (as per the mesh coordinate specification)
function getBlockMatchAtPoint(A, Am, i, j, B, Bm, block_size, search_r)
	b_rad = block_size + search_r;

	# convert to matrix coordinates
	Ai = round(Int64, i - Am.disp[1]);
	Aj = round(Int64, j - Am.disp[2]);

	Bi = round(Int64, i - Bm.disp[1]);
	Bj = round(Int64, j - Bm.disp[2]);

	if (!isInternal(A, Ai, Aj, block_size) || !isInternal(B, Bi, Bj, b_rad))
		return noMatch;
	end

	Ai_range = ceil(Ai)-block_size:ceil(Ai)+block_size;
	Aj_range = ceil(Aj)-block_size:ceil(Aj)+block_size;
	Bi_range = ceil(Bi)-b_rad:ceil(Bi)+b_rad;
	Bj_range = ceil(Bj)-b_rad:ceil(Bj)+b_rad;

	xc = normxcorr2(sub(A, Ai_range, Aj_range), sub(B, Bi_range, Bj_range));
	r_max = maximum(xc); ind = find(xc .== r_max);
	if (size(ind) == 0) return noMatch; end
	(i_max, j_max) = (rem(ind[1], size(xc, 2)), cld(ind[1], size(xc, 2)));

	return [i_max - 1 - search_r; j_max - 1 - search_r; r_max];
end

function Meshes2Matches(A, Am, B, Bm, block_size, search_r, min_r)

	src_mesh = Am;
	dst_mesh = Bm;

	n = 0;
	n_upperbound = Am.n;

	src_pointIndices = Array(Int64, n_upperbound);
	dst_points = Points(n_upperbound);
	dst_triangles = Triangles(n_upperbound);
	dst_weights = Weights(n_upperbound);
	dispVectors = Points(n_upperbound);

	for j in 1:n_upperbound
		(Ai, Aj) = Am.nodes[j]
		v = getBlockMatchAtPoint(A, Am, Ai, Aj, B, Bm, block_size, search_r);
		if v[3] < min_r continue; end
		
		dispVector = v[1:2];
		dst_point = Am.nodes[j] + dispVector;
		dst_triangle = findMeshTriangle(Bm, dst_point[1], dst_point[2]); 
		if dst_triangle == (0, 0, 0) continue; end
		n += 1;
		src_pointIndices[n] = j;
		dispVectors[n] = dispVector;
		dst_points[n] = Am.nodes[j] + dispVectors[n];
		dst_triangles[n] = dst_triangle;
		dst_weights[n] = getTriangleWeights(Bm, dst_triangle, dst_point[1], dst_point[2]);
	end

	src_pointIndices = src_pointIndices[1:n];
	dst_points = dst_points[1:n];
	dispVectors = dispVectors[1:n];
	dst_triangles = dst_triangles[1:n];
	dst_weights = dst_weights[1:n]; 

	println("$n total matches found of $n_upperbound mesh points.");

	matches = Matches(src_mesh, dst_mesh, n, src_pointIndices, dst_points, dst_triangles, dst_weights, dispVectors);
	return matches;
end

# find the triangular mesh indices for a given point in A
function findMeshTriangle(Am, i, j)

	# convert to A's local coordinates, and displace by the mesh offset
	Ai = i - Am.disp[1] - Am.offsets[1];
	Aj = j - Am.disp[2] - Am.offsets[2];

	# find which rows the point belongs to
	i0 = round(Int64, Ai / Am.dists[1] + 1);

	if isodd(i0)
		j0 = round(Int64, Aj / Am.dists[2] + 1);
	else
		j0 = round(Int64, Aj / Am.dists[2] + 0.5);
	end

	ind0 = getMeshIndex(Am.dims, i0, j0);
	
	if ind0 == 0
		return noTriangle;
	end

	node0 = Am.nodes[ind0];

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
function getTriangleWeights(Am, triangle, pi, pj)
	R = vcat(Am.nodes[triangle[1]]', Am.nodes[triangle[2]]', Am.nodes[triangle[3]]')
	R = hcat(R, ones(Float64, 3, 1));
	r = hcat(pi, pj, 1.0);

	V = r * R^-1;

	return (V[1], V[2], V[3]);
end

function findNode(Ms, mesh_ind, node_ind);
	return Ms.nodes_indices[mesh_ind] + node_ind;
end

function findIndex(Ms, mesh_index_tuple)
	return findfirst(this -> mesh_index_tuple == this.index, Ms.meshes)
end

function makeNewMeshSet()
	N = 0;
	M = 0;

	indices = Array{Index, 1}(0);
	
	n = 0;
	m = 0;
	m_i = 0;
	m_e = 0;

	meshes = Array{Mesh, 1}(0);
	nodes_indices = Array{Int64, 1}(0);

	matches = Array{Matches, 1}(0);
	matches_pairs = Array{Pair, 1}(0);

	return MeshSet(N, M, indices, n, m, m_i, m_e, meshes, nodes_indices, matches, matches_pairs);
end

function addMesh2MeshSet!(Am, Ms)
	Ms.N += 1;

	push!(Ms.indices, Am.index);

	Ms.m_i += Am.m;
	Ms.m += Am.m;

	push!(Ms.meshes, Am);
	if length(Ms.nodes_indices) == 0 push!(Ms.nodes_indices, 0);
	else push!(Ms.nodes_indices, Ms.n);
	end

	Ms.n += Am.n;
	return;
end

function addMatches2MeshSet!(M, Ms)
	Ms.M += 1;
	
	Ms.n;
	Ms.m += M.n;
	Ms.m_e += M.n;

	push!(Ms.matches, M);
	push!(Ms.matches_pairs, (findIndex(Ms, M.src_mesh.index), findIndex(Ms, M.dst_mesh.index)));
	return;
end

function solveMeshSet!(Ms, eta, n_steps, n_grad, show_plot)
	nodes = Points(0);
	nodes_fixed = BinaryProperty(0);
	edges = spzeros(Float64, Ms.n, 0);
	edges_padded = Array{Edges}(0); 
	edge_lengths = FloatProperty(0);
	edge_coeffs = FloatProperty(0);
	for i in 1:Ms.N
		cur_mesh = Ms.meshes[i];
		append!(nodes, cur_mesh.nodes);
		append!(nodes_fixed, cur_mesh.nodes_fixed);
		append!(edge_lengths, cur_mesh.edge_lengths);
		append!(edge_coeffs, cur_mesh.edge_coeffs);
		if (i == Ms.N)	
			push!(edges_padded, vcat(spzeros(Float64, Ms.nodes_indices[i], cur_mesh.m), 
						 cur_mesh.edges));
		else
			push!(edges_padded, vcat(spzeros(Float64, Ms.nodes_indices[i], cur_mesh.m), 
						 cur_mesh.edges, 
						 spzeros(Float64, Ms.n - Ms.nodes_indices[i] - cur_mesh.n, cur_mesh.m)));
		end
	end

	for i in 1:Ms.M
		cur_matches = Ms.matches[i];
		append!(edge_lengths, fill(0.0, cur_matches.n));
		append!(edge_coeffs, fill(convert(Float64, match_coeff), cur_matches.n));
		edges_toadd = spzeros(Float64, Ms.n, cur_matches.n);

		for j in 1:Ms.matches[i].n
			edges_toadd[findNode(Ms, Ms.matches_pairs[i][1], cur_matches.src_pointIndices[j]), j] = -1;
			edges_toadd[findNode(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][1]), j] = cur_matches.dst_weights[j][1];
			edges_toadd[findNode(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][2]), j] = cur_matches.dst_weights[j][2];
			edges_toadd[findNode(Ms, Ms.matches_pairs[i][2], cur_matches.dst_triangles[j][3]), j] = cur_matches.dst_weights[j][3];
		end
		push!(edges_padded, edges_toadd);
	end

	nodes = hcat(nodes...);
	edges = hcat(edges_padded...);

	SolveMesh!(nodes, nodes_fixed, edges, edge_coeffs, edge_lengths, eta, n_steps, n_grad, show_plot);
	nodes_t = Points(0);
	for i in 1:size(nodes, 2)
       		push!(nodes_t, vec(nodes[:, i]))
       	end
	for i in 1:Ms.N
		cur_mesh = Ms.meshes[i];
		cur_mesh.nodes_t = nodes_t[Ms.nodes_indices[i] + 1:cur_mesh.n];
	end

end

function MeshSet2JLD(filename, Ms)
	jldopen(filename, "w") do file
		write(file, "MeshSet", Ms);
	end
end

function JLD2MeshSet(filename)
	jldopen(filename, "r") do file
		Ms = file["MeshSet"];
	end

	return Ms
end

end

function test()
	Ap = "./EM_images/Tile_r4-c2_S2-W001_sec20.tif";
	dAi = 21906;
	dAj = 36429;

	Bp = "./EM_images/Tile_r4-c3_S2-W001_sec20.tif";
	dBi = 29090; # 2908.6;
	dBj = 36251; # 3624.3;

	block_size = 20;
	search_r = 80;
	min_r = 0.65;
	mesh_length = 100;
	mesh_coeff = 0.5;
	match_coeff = 1.0;
	eta = 0.15
	n_steps = 80;
	n_grad = 40;

	@time Am = Tile2Mesh(Ap, (1, 2, 42), (4, 2), dAi, dAj, false, mesh_length, mesh_coeff);
	@time Bm = Tile2Mesh(Bp, (1, 2, 43), (4, 3), dBi, dBj, false, mesh_length, mesh_coeff);
	@time A = getMeshImage(Am);
	@time B = getMeshImage(Bm);
	@time Mab = Meshes2Matches(A, Am, B, Bm, block_size, search_r, min_r);
	@time Mba = Meshes2Matches(B, Bm, A, Am, block_size, search_r, min_r);
	Ms = makeNewMeshSet();
	@time addMesh2MeshSet!(Am, Ms);
	@time addMesh2MeshSet!(Bm, Ms);
	@time addMatches2MeshSet!(Mab, Ms);
	@time addMatches2MeshSet!(Mba, Ms);
	@time solveMeshSet!(Ms, eta, n_steps, n_grad, true)
end