type Matches
	src_index::Index					# source mesh index
	dst_index::Index					# destination mesh index

	n::Int64						# number of matches

	src_points_indices::Array{Int64, 1}			# index of points in src_mesh.points
	dst_points::Points					# location of points in the destination
	dst_triangles::Triangles				# index of the triangles
	dst_weights::Weights 	# barycentric weights for respective triangle index

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
	return [i_max - 1 - rad; j_max - 1 - rad; r_max], xc;
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
	disp_vectors_mags_i = Array{Float64, 1}(0);
	disp_vectors_mags_j = Array{Float64, 1}(0);
	disp_vectors_mags_f = Array{Float64, 1}(0);
#=
	A_im_array = Array{Array{UInt8, 2}, 1}(n_upperbound);
	B_im_array = Array{Array{UInt8, 2}, 1}(n_upperbound);
	xc_im_array = Array{Array{Float64, 2}, 1}(n_upperbound);
=#
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
#=
	if is_pre_aligned(Am.index)
		blockmatch_impath = joinpath(ALIGNED_DIR, "blockmatches", string(Am.name, "-", Bm.name));
        else
		blockmatch_impath = joinpath(MONTAGED_DIR, "blockmatches", string(Am.name, "-", Bm.name));
        end
	if !isdir(blockmatch_impath)
		mkdir(blockmatch_impath);
	end
=#

	inc_total() = (n_total += 1;)
	inc_not_enough_dyn_range() = (n_not_enough_dyn_range += 1;)

	#if (is_prealigned(Am.index) && is_prealigned(Bm.index))
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
					#A_im = A[src_ranges[idx][1], src_ranges[idx][2]];
				#	B_im = B[dst_ranges[idx][1], dst_ranges[idx][2]];
					max_vect_xc = remotecall_fetch(p, get_max_xc_vector, A[src_ranges[idx][1], src_ranges[idx][2]],  B[dst_ranges[idx][1], dst_ranges[idx][2]]);
				#	max_vect_xc = remotecall_fetch(p, get_max_xc_vector, A_im, B_im); 
					disp_vectors_raw[idx] = max_vect_xc[1];
				#	xc_im_array[idx] = (max_vect_xc[2] .+ 1)./ 2;
	#				A_im_array[idx] = A_im;	
	#				B_im_array[idx] = B_im;	

					#println("$p: Matched point $idx, with displacement vector $(disp_vectors_raw[idx])");
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


		if v[3] < min_r; 
		  
		n_low_r +=1; 
	#= 	imwrite(grayim((A_im_array[idx]/255)'), joinpath(blockmatch_impath, string("bad_low_r_", n_low_r,"_src.jpg")));
	 	imwrite(grayim((B_im_array[idx]/255)'), joinpath(blockmatch_impath, string("bad_low_r_", n_low_r,"_dst.jpg")));
if (!isnan(sum(xc_im_array[idx])))	 
  imwrite(grayim(xc_im_array[idx]'), joinpath(blockmatch_impath, string("bad_low_r_", n_low_r,"_xc.jpg")));
end=#
		continue; end

		disp_vector = v[1:2];
		if norm(disp_vector) > mu + 3.5 * sigma; 
		n_outlier +=1; 
	      
#=	 	imwrite(grayim((A_im_array[idx]/255)'), joinpath(blockmatch_impath, string("bad_outlier_", n_outlier,"_src.jpg")));
	 	imwrite(grayim((B_im_array[idx]/255)'), joinpath(blockmatch_impath, string("bad_outlier_", n_outlier,"_dst.jpg")));
if (!isnan(sum(xc_im_array[idx])))	 
	 	imwrite(grayim((xc_im_array[idx])'), joinpath(blockmatch_impath, string("bad_outlier_", n_outlier,"_xc.jpg")));
	      end=#
	      	continue; end


		dst_point = Am.nodes[idx] + disp_vector;
		dst_triangle = findMeshTriangle(Bm, dst_point[1], dst_point[2]); 
		if dst_triangle == NO_TRIANGLE n_no_triangle +=1; 
	#=	
	 	imwrite(grayim((A_im_array[idx]/255)'), joinpath(blockmatch_impath, string("bad_triangle_", n_no_triangle,"_src.jpg")));
	 	imwrite(grayim((B_im_array[idx]/255)'), joinpath(blockmatch_impath, string("bad_triangle_", n_no_triangle,"_dst.jpg")));
if (!isnan(sum(xc_im_array[idx])))	 
	 	imwrite(grayim((xc_im_array[idx])'), joinpath(blockmatch_impath, string("bad_triangle_", n_no_triangle,"_xc.jpg")));
		
	      end		=#
		continue; end
		n += 1;#=
	 	imwrite(grayim((A_im_array[idx]/255)'), joinpath(blockmatch_impath, string("accepted_", n,"_src.jpg")));
	 	imwrite(grayim((B_im_array[idx]/255)'), joinpath(blockmatch_impath, string("accepted_", n,"_dst.jpg")));
if (!isnan(sum(xc_im_array[idx])))	 
	 	imwrite(grayim((xc_im_array[idx])'), joinpath(blockmatch_impath, string("accepted_", n,"_xc.jpg")));
	      end=#
		push!(src_points_indices, idx);
		push!(disp_vectors, disp_vector);
		push!(disp_vectors_mags_f, norm(disp_vector));
		push!(disp_vectors_mags_i, disp_vector[1]);
		push!(disp_vectors_mags_j, disp_vector[2]);
		push!(dst_points, dst_point);
		push!(dst_triangles, dst_triangle);
		push!(dst_weights, getTriangleWeights(Bm, dst_triangle, dst_point[1], dst_point[2]));
	end

	mu_f = mean(disp_vectors_mags_f);
	sigma_f = std(disp_vectors_mags_f);
	max_f = maximum(disp_vectors_mags_f);
	mu_i = mean(disp_vectors_mags_i);
	sigma_i = std(disp_vectors_mags_i);
	max_i = maximum(disp_vectors_mags_i);
	mu_j = mean(disp_vectors_mags_j);
	sigma_j = std(disp_vectors_mags_j);
	max_j = maximum(disp_vectors_mags_j);

	if n == 0
	return Void;
	end

	println("###\n$p1 -> $p2: $n_upperbound in mesh, $n_total in overlap, $n accepted.\nRejections: $n_not_enough_dyn_range (low dynamic range), $n_low_r (low r), $n_outlier (outliers), $n_no_triangle (outside triangles).\nDisplacement statistics:\nNorms, before filtering: mean = $mu, sigma = $sigma, max = $max\nNorms, after filtering: mean = $mu_f, sigma = $sigma_f, max = $max_f\ni-coord. after filtering: mean = $mu_i, sigma = $sigma_i, max = $max_i\nj-coord. after filtering: mean = $mu_j, sigma = $sigma_j, max = $max_j\n###");

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
