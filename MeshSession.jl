################################# SCRIPT FOR TESTING ###################################

#params
block_size = 40;
search_r = 80;
min_r = 0.50;
mesh_length = 200;
mesh_coeff = 0.5;
match_coeff = 1.0;
eta = 0.10
n_steps = 120;
n_grad = 80;
show_plot = false;
num_procs = 4;
num_tiles = 16;

Ms = MeshModule.makeNewMeshSet();
sr = readdlm("section.txt");

#meshArray = Array{Mesh, 1}(0);
imageArray = Array{Array{Float64, 2}, 1}(0);

@time for i in 1:size(sr, 1)
	path = string("./EM_images/", sr[i, 1]);
	println(path);
	di = sr[i, 2];
	dj = sr[i, 3];
	MeshModule.addMesh2MeshSet!(MeshModule.Tile2Mesh(path, (1, 20, i), (0, 0), di, dj, false, mesh_length, mesh_coeff), Ms);
	meshImage = MeshModule.getMeshImage(Ms.meshes[i]);
	push!(imageArray, meshImage);
end
																																																																																																	


#=
PLoop = @sync @parallel for i=1:Ms.N, j in 1:Ms.N
		if j == i continue; end
	MeshModule.Meshes2Matches(imageArray[i], Ms.meshes[i], imageArray[j], Ms.meshes[j], block_size, search_r, min_r);
end
	for i in 1:size(PLoop)
	MeshModule.addMatches2MeshSet!(fech(PLoop[i]), Ms);
	end
=#

@time for i in 1:Ms.N
	for j in 1:Ms.N
		if j == i continue; end
	M = MeshModule.Meshes2Matches(imageArray[i], Ms.meshes[i], imageArray[j], Ms.meshes[j], block_size, search_r, min_r);
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
@time MeshModule.solveMeshSet!(Ms, match_coeff, eta, n_steps, n_grad, show_plot);
@time MeshModule.MeshSet2JLD("solvedMesh(1,20,0).jld", Ms);
