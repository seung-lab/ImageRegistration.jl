#=Tile_r1-c1_S2-W001_sec20.tif 14628 14628
Tile_r1-c2_S2-W001_sec20.tif 21815 14420
Tile_r1-c3_S2-W001_sec20.tif 28980 14292
Tile_r1-c4_S2-W001_sec20.tif 35871 14145
Tile_r2-c1_S2-W001_sec20.tif 14659 21795
Tile_r2-c2_S2-W001_sec20.tif 21880 21706
Tile_r2-c3_S2-W001_sec20.tif 28916 21564
Tile_r2-c4_S2-W001_sec20.tif 36075 21565
Tile_r3-c1_S2-W001_sec20.tif 14621 29135
Tile_r3-c2_S2-W001_sec20.tif 21851 28960
Tile_r3-c3_S2-W001_sec20.tif 29125 28892
Tile_r3-c4_S2-W001_sec20.tif 36176 28852
Tile_r4-c1_S2-W001_sec20.tif 14602 36500
Tile_r4-c2_S2-W001_sec20.tif 21907 36430
Tile_r4-c3_S2-W001_sec20.tif 29091 36252
Tile_r4-c4_S2-W001_sec20.tif 36259 36183
=#
################################# SCRIPT FOR TESTING ###################################

Ap = "./EM_images/Tile_r4-c2_S2-W001_sec20.tif";
dAi = 21906;
dAj = 36429;

Bp = "./EM_images/Tile_r4-c3_S2-W001_sec20.tif";
dBi = 29090; # 2908.6;
dBj = 36251; # 3624.3;

block_size = 40;
search_r = 100;
min_r = 0.55;
mesh_length = 200;
mesh_coeff = 0.5;
match_coeff = 1.0;
eta = 0.10
n_steps = 120;
n_grad = 80;

@time Am = MeshModule.Tile2Mesh(Ap, (1, 2, 42), (4, 2), dAi, dAj, false, mesh_length, mesh_coeff);
@time Bm = MeshModule.Tile2Mesh(Bp, (1, 2, 43), (4, 3), dBi, dBj, false, mesh_length, mesh_coeff);
@time A = MeshModule.getMeshImage(Am);
@time B = MeshModule.getMeshImage(Bm);
@time Mab = MeshModule.Meshes2Matches(A, Am, B, Bm, block_size, search_r, min_r);
@time Mba = MeshModule.Meshes2Matches(B, Bm, A, Am, block_size, search_r, min_r);
Ms = MeshModule.makeNewMeshSet();
@time MeshModule.addMesh2MeshSet!(Am, Ms);
@time MeshModule.addMesh2MeshSet!(Bm, Ms);
@time MeshModule.addMatches2MeshSet!(Mab, Ms);
@time MeshModule.addMatches2MeshSet!(Mba, Ms);
@time MeshModule.solveMeshSet!(Ms, match_coeff, eta, n_steps, n_grad, false);
@time MeshModule.MeshSet2JLD("solvedMesh.jld", Ms);
