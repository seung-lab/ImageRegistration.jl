using Julimaps
using Params
using MeshModule

################################# SCRIPT FOR TESTING ###################################
tic();

max_tile_size = 0;
Ms = MeshModule.makeNewMeshSet();
sr = readdlm("W001_sec21_offsets.txt");

#sr = readdlm(joinpath("input_images", "W001_sec20", "W001_sec20_offsets.txt"))


@time for i in 1:num_tiles
	path = string("./W001_sec21/", sr[i, 1]);
	row = sr[i, 2];
	col = sr[i, 3];
	di = sr[i, 4];
	dj = sr[i, 5];
	MeshModule.addMesh2MeshSet!(MeshModule.Tile2Mesh(path, (1, 21, row, col), di, dj, false, mesh_length, mesh_coeff), Ms);
	meshImage = MeshModule.getMeshImage(Ms.meshes[i]);
	max_size = max(size(meshImage, 1), size(meshImage, 2));
	if max_tile_size < max_size max_tile_size = max_size; end
end
	imageArray = SharedArray(Float64, max_tile_size, max_tile_size, num_tiles);

@time for k in 0:num_procs:num_tiles
	@sync @parallel for l in 1:num_procs
	i = k+l;
	if i > num_tiles return; end;
	meshImage = MeshModule.getMeshImage(Ms.meshes[i]);
	imageArray[1:size(meshImage, 1), 1:size(meshImage, 2), i] = meshImage;
	end
end

print("Initialisation, "); toc(); println();

tic();
adjacent_pairs = Pairings(0);
diagonal_pairs = Pairings(0);

for i in 1:Ms.N, j in 1:Ms.N
	if MeshModule.isAdjacent(Ms.meshes[i], Ms.meshes[j]) push!(adjacent_pairs, (i, j)); end
	if MeshModule.isDiagonal(Ms.meshes[i], Ms.meshes[j]) push!(diagonal_pairs, (i, j)); end
end

pairs = vcat(adjacent_pairs, diagonal_pairs);

@time for k in 0:num_procs:length(pairs)
	toFetch = @sync @parallel for l in 1:num_procs
	ind = k + l;
	if ind > length(pairs) return; end
	(i, j) = pairs[ind];
	return MeshModule.Meshes2Matches(imageArray[:, :, i], Ms.meshes[i], imageArray[:, :, j], Ms.meshes[j], block_size, search_r, min_r);
	end
	for i = 1:length(toFetch)
		M = fetch(toFetch[i])
		if typeof(M) == Void continue; end
		MeshModule.addMatches2MeshSet!(M, Ms);
	end
end

print("Blockmatching, "); toc(); println();



@time MeshModule.solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);

residuals = Points(0);
residuals_t = Points(0);


@time for k in 60:70
src_p = Points(0);
dst_p = Points(0);
for i in 1:Ms.matches[k].n
		w = Ms.matches[k].dst_weights[i];
		t = Ms.matches[k].dst_triangles[i];
		p = Ms.matches[k].src_pointIndices[i];
		src = Ms.meshes[MeshModule.findIndex(Ms, Ms.matches[k].src_index)]
		dst = Ms.meshes[MeshModule.findIndex(Ms, Ms.matches[k].dst_index)]
		p1 = src.nodes[p];
		p2 = dst.nodes[t[1]] * w[1] + dst.nodes[t[2]] * w[2] + dst.nodes[t[3]] * w[3]
		push!(src_p, p1);
		push!(dst_p, p2);
end

src_im = imageArray[:, :, Ms.matches_pairs[k][1]];
dst_im = imageArray[:, :, Ms.matches_pairs[k][2]];

square_size = 2*block_size + 1;
cols = round(Int64, sqrt(Ms.matches[k].n));


sbs_rows = Array{Array{Float64, 2}, 1}(0);
for r in 1:cols:Ms.matches[k].n
sbs_row = zeros(Float64, 3*square_size, 0);
	for c in 1:cols
		i = r+c;
		if i > Ms.matches[k].n continue; end;
		iv = round(Int64, src_p[i][1]-Ms.meshes[Ms.matches_pairs[k][1]].disp[1])
		jv = round(Int64, src_p[i][2]-Ms.meshes[Ms.matches_pairs[k][1]].disp[2])
		src_block = src_im[iv-block_size:iv+block_size, jv-block_size:jv+block_size];
		iv = round(Int64, dst_p[i][1]-Ms.meshes[Ms.matches_pairs[k][2]].disp[1])
		jv = round(Int64, dst_p[i][2]-Ms.meshes[Ms.matches_pairs[k][2]].disp[2])
		dst_block = dst_im[iv-block_size:iv+block_size, jv-block_size:jv+block_size];
		sbs = vcat(src_block, dst_block, 1 - abs(src_block-dst_block));
		sbs_row = hcat(sbs_row, zeros(Float64, 3*square_size, 2), sbs);
	end
push!(sbs_rows, sbs_row)
end

row_size = size(sbs_rows[1], 2);

sbs_total = zeros(Float64, 0, row_size);

for i in 1:size(sbs_rows,1)-1

sbs_total = vcat(sbs_total, sbs_rows[i], zeros(Float64, 10, row_size));

end

lastrow = sbs_rows[size(sbs_rows, 1)];
lastrow = hcat(lastrow, zeros(Float64, 3 * square_size, row_size-size(lastrow, 2)));

sbs_total = vcat(sbs_total, lastrow);

view(sbs_total);

end




#=
for k in 1:Ms.M
	for i in 1:Ms.matches[k].n
		w = Ms.matches[k].dst_weights[i];
		t = Ms.matches[k].dst_triangles[i];
		p = Ms.matches[k].src_pointIndices[i];
		src = Ms.meshes[MeshModule.findIndex(Ms, Ms.matches[k].src_index)]
		dst = Ms.meshes[MeshModule.findIndex(Ms, Ms.matches[k].dst_index)]
		p1 = src.nodes_t[p];
		p2 = dst.nodes_t[t[1]] * w[1] + dst.nodes_t[t[2]] * w[2] + dst.nodes_t[t[3]] * w[3]
		push!(residuals_t, p2-p1);
		p1 = src.nodes[p];
		p2 = dst.nodes[t[1]] * w[1] + dst.nodes[t[2]] * w[2] + dst.nodes[t[3]] * w[3]
		push!(residuals, p2-p1);
	end
end
=#
@time MeshModule.MeshSet2JLD("solvedMesh(1,21,0).jld", Ms);


####### LEGACY CODE FOR PAIR TESTING #############

#Ap = "./EM_images/Tile_r4-c2_S2-W001_sec20.tif";
#dAi = 21906;
#dAj = 36429;

#Bp = "./EM_images/Tile_r4-c3_S2-W001_sec20.tif";
#dBi = 10000#29090; # 2908.6;
#dBj = 10000#36251; # 3624.3;


#=

Ms = makeNewMeshSet();
@time Am = MeshModule.Tile2Mesh(Ap, (1, 2, 42), (4, 2), dAi, dAj, false, mesh_length, mesh_coeff);
@time Bm = MeshModule.Tile2Mesh(Bp, (1, 2, 43), (4, 3), dBi, dBj, false, mesh_length, mesh_coeff);
@time A = MeshModule.getMeshImage(Am);
@time B = MeshModule.getMeshImage(Bm);
@time Mab = MeshModule.Meshes2Matches(A, Am, B, Bm, block_size, search_r, min_r);
@time Mba = MeshModule.Meshes2Matches(B, Bm, A, Am, block_size, search_r, min_r);

@time MeshModule.addMesh2MeshSet!(Am, Ms);
@time MeshModule.addMesh2MeshSet!(Bm, Ms);
@time MeshModule.addMatches2MeshSet!(Mab, Ms);
@time MeshModule.addMatches2MeshSet!(Mba, Ms);
=#
