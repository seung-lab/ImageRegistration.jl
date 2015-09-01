type Matches
	src_index::Index					# source mesh index
	dst_index::Index					# destination mesh index

	n::Int64						# number of matches

	src_points_indices::Array{Int64, 1}			# index of points in src_mesh.points
	dst_points::Points					# location of points in the destination
	dst_triangles::Triangles				# index of the triangles
	dst_weights::Weights					# weights at each

	disp_vectors::Points					# displacement vector src->dest
end

function is_internal(A, pt::Point, disp::Point, d)
	return is_internal(A, pt - disp, d); 
end

function is_internal(A, pt::Point, d)
	if pt[1] > d && pt[1] <= (size(A, 1) - d) && pt[2] > d && pt[2] <= (size(A, 2) - d)
		return true
	end
	return false
end

function get_range(A, pt::Point, disp::Point, d)
	if !is_internal(A, pt, disp, d)
	  return NO_RANGE;
	end

	pt_local = pt - disp;
	Ai = round(Int64, ceil(pt_local[1]));
	Aj = round(Int64, ceil(pt_local[2]));

	Ai_range = Ai-d:Ai+d;
	Aj_range = Aj-d:Aj+d;

	return (Ai_range, Aj_range);
end

function get_max_xc_vector(A, B)

	xc = normxcorr2(A, B);
	r_max = maximum(xc); 
	rad = round(Int64, (size(xc, 1) - 1)/ 2)

	ind = findfirst(r_max .== xc);
	
	x1 = size(xc, 1);
	x2 = size(xc, 2);

	if ind == 0 return NO_MATCH; end
	(i_max, j_max) = (rem(ind, size(xc, 1)), cld(ind, size(xc, 1)));
	if i_max == 0 i_max = size(xc, 1); end
	return [i_max - 1 - rad; j_max - 1 - rad; r_max];
end

function Meshes2Matches(A, Am::Mesh, B, Bm::Mesh, params::Params)

	if (Am==Bm)
		return Void;
	end

	src_index = Am.index;
	dst_index = Bm.index;
	p1 = Am.name;
	p2 = Bm.name;

	n = 0;
	n_total = 0;
	n_low_r = 0;
	n_outlier = 0;
	n_no_triangle = 0;
	n_not_enough_dyn_range = 0;
	n_upperbound = Am.n;

	src_points_indices = Array(Int64, 0);
	dst_points = Points(0);
	dst_triangles = Triangles(0);
	dst_weights = Weights(0);
	disp_vectors = Points(0);
	disp_vectors_raw = Array{Array{Float64, 1}, 1}(n_upperbound);
	disp_vectors_mags = Array{Float64, 1}(0);

	src_ranges = Array{Tuple{UnitRange{Int64}, UnitRange{Int64}}, 1}(n_upperbound);
	dst_ranges = Array{Tuple{UnitRange{Int64}, UnitRange{Int64}}, 1}(n_upperbound);

	min_dyn_range_ratio = params.min_dyn_range_ratio;
	block_size = params.block_size;
	search_r = params.search_r;
	min_r = params.min_r;
	b_rad = block_size + search_r;

	
	# preprocessing
	for idx in 1:n_upperbound
	pt = Am.nodes[idx];
	src_ranges[idx] = get_range(A, pt, Am.disp, block_size);
	dst_ranges[idx] = get_range(B, pt, Bm.disp, b_rad);
	end
	

	inc_total() = (n_total += 1;)
	inc_not_enough_dyn_range() = (n_not_enough_dyn_range += 1;)

	#if (is_pre_aligned(Am.index) && is_pre_aligned(Bm.index))
	k = 1;
	nextidx() = (idx=k; k+=1; idx);
	@sync begin
	for p in 1:num_procs
		if p != myid() || num_procs == 1
			@async begin
				while true
					idx = nextidx();
						if idx > n_upperbound
							break
						end
					if src_ranges[idx] == NO_RANGE || dst_ranges[idx] == NO_RANGE
						disp_vectors_raw[idx] = NO_MATCH;
						continue;
					end
					inc_total();
					if maximum(A[src_ranges[idx][1], src_ranges[idx][2]]) / minimum(A[src_ranges[idx][1], src_ranges[idx][2]]) < min_dyn_range_ratio
						disp_vectors_raw[idx] = NO_MATCH;
						inc_not_enough_dyn_range();
						continue;
					end
					disp_vectors_raw[idx] = remotecall_fetch(p, get_max_xc_vector, A[src_ranges[idx][1], src_ranges[idx][2]],  B[dst_ranges[idx][1], dst_ranges[idx][2]]);
					println("$p: Matched point $idx, with displacement vector $(disp_vectors_raw[idx])");
				end
			end
		end
	end
      	end
	#=else
	for idx in 1:n_upperbound
		if src_ranges[idx] == NO_RANGE || src_range[idx] == NO_RANGE
			disp_vectors_raw[idx] = NO_MATCH;
			continue;
		end
		inc_total();
		if maximum(A[src_ranges[idx][1], src_ranges[idx][2]]) / minimum(A[src_ranges[idx][1], src_ranges[idx][2]]) < min_dyn_range_ratio
			disp_vectors_raw[idx] = NO_MATCH;
			inc_not_enough_dyn_range();
			continue;
		end
		disp_vectors_raw[idx] = get_max_xc_vector(A[src_ranges[idx][1], src_ranges[idx][2]],  B[dst_ranges[idx][1], dst_ranges[idx][2]]);
	end
	end=#

	for idx in 1:n_upperbound
	  	v = disp_vectors_raw[idx];
		if v != NO_MATCH && v[3] >= min_r
			push!(disp_vectors_mags, norm(v[1:2]));
		end
	end

	if length(disp_vectors_mags) == 0
	return Void;
	end


	mu = mean(disp_vectors_mags);
	sigma = std(disp_vectors_mags);
	max = maximum(disp_vectors_mags);

	for idx in 1:n_upperbound
		v = disp_vectors_raw[idx];	
		if v == NO_MATCH continue; end	
		if v[3] < min_r; n_low_r +=1; continue; end
		disp_vector = v[1:2];
		if norm(disp_vector) > mu + 2.5 * sigma; n_outlier +=1; continue; end
		dst_point = Am.nodes[idx] + disp_vector;
		dst_triangle = findMeshTriangle(Bm, dst_point[1], dst_point[2]); 
		if dst_triangle == NO_TRIANGLE n_no_triangle +=1; continue; end
		n += 1;
		push!(src_points_indices, idx);
		push!(disp_vectors, disp_vector);
		push!(dst_points, dst_point);
		push!(dst_triangles, dst_triangle);
		push!(dst_weights, getTriangleWeights(Bm, dst_triangle, dst_point[1], dst_point[2]));
	end

	if n == 0
	return Void;
	end

	println("###\n$p1 -> $p2: $n_upperbound in mesh, $n_total in overlap, $n accepted.\nRejections: $n_not_enough_dyn_range (low dynamic range), $n_low_r (low r), $n_outlier (outliers), $n_no_triangle (outside triangles).\nDisplacement mean = $mu, sigma = $sigma, max = $max\n###");

	matches = Matches(src_index, dst_index, n, src_points_indices, dst_points, dst_triangles, dst_weights, disp_vectors);
	return matches;
end

# multiple dispatch for when called remotely - A, B are remote references to the mesh
function Meshes2Matches(A, A_rr::RemoteRef, B, B_rr::RemoteRef, params_rr::RemoteRef)
Am = take!(A_rr);
Bm = take!(B_rr);
params = take!(params_rr);
return Meshes2Matches(A, Am, B, Bm, params);
end
