################################# SCRIPT FOR TESTING ###################################


Ms = MeshModule.makeNewMeshSet();
sr = readdlm("section.txt");

#params
block_size = 20;
search_r = 80; #20;
min_r = 0.65; #0.50;
mesh_length = 200;
mesh_coeff = 1;
match_coeff = 10;
eta = 0.01;
n_steps = 2000;
n_grad = 1950;
show_plot = false;
num_procs = 12;
grad_threshold = 1/1000;
n_newton = 50;
max_tile_size = 10000;
num_tiles = size(sr, 1);

imageArray = SharedArray(Float64, 10000, 10000, num_tiles);

@time for i in 1:num_tiles
	path = string("./EM_images/", sr[i, 1]);
	println(path);
	di = sr[i, 2];
	dj = sr[i, 3];
	MeshModule.addMesh2MeshSet!(MeshModule.Tile2Mesh(path, (1, 20, i), (0, 0), di, dj, false, mesh_length, mesh_coeff), Ms);
	meshImage = MeshModule.getMeshImage(Ms.meshes[i]);
	imageArray[1:size(meshImage, 1), 1:size(meshImage, 2), i] = meshImage;
end

@time for k in 0:num_procs:Ms.N^2
	toFetch = @sync @parallel for l in num_procs
	(i, j) = (rem(k+l, Ms.N), cld(k+l, Ms.N));
	if i == 0 i = Ms.N; end
	if i == j || i > Ms.N || j > Ms.N continue end;
		return MeshModule.Meshes2Matches(imageArray[:, :, i], Ms.meshes[i], imageArray[:, :, j], Ms.meshes[j], block_size, search_r, min_r);
	end
	for i = 1:length(toFetch)
		MeshModule.addMatches2MeshSet!(fetch(toFetch[1]), Ms);
	end
end


@time for i in 1:Ms.N
	for j in 1:Ms.N
		if j == i continue; end
	@time M = MeshModule.Meshes2Matches(imageArray[i], Ms.meshes[i], imageArray[j], Ms.meshes[j], block_size, search_r, min_r);
	if M.n == 0 continue; end
	MeshModule.addMatches2MeshSet!(M, Ms);
	end
end

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
@time MeshModule.solveMeshSet!(Ms, match_coeff, eta, grad_threshold, n_newton);
@time MeshModule.MeshSet2JLD("solvedMesh(1,20,0).jld", Ms);
