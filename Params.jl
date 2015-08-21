module Params

function parseName(name::String)
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", name)
	return parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2]) # wafer, section, row, column
end


function waferpaths2dict(wafer_path_file)
	warray = readdlm(wafer_path_file)
	wdict = Dict()
	for (idx, path) in zip(warray[:,1], warray[:,2])
		wdict[idx] = path
	end
	return wdict
end

function parseRoughAlign(info_path::String)
	file = readdlm(info_path);
	session = cell(size(file, 1), 4); # name, index, dx, dy
	for i in 1:size(file, 1)
	m = match(r"(Tile\S+).tif", file[i, 1]);
	session[i, 1] = m[1];
	session[i, 2] = parseName(session[i, 1]);
	session[i, 3] = file[i, 2];
	session[i, 4] = file[i, 3];
	end
	return session;
end


bucket_dir_path = "/usr/people/dih/seungmount";
affine_dir_path = "~";
dataset_dir_path = "./datasets/piriform"
wafer_filename = "wafer_paths.txt";
rough_align_filename = "rough_align_test1-10.txt";

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, SESSION

global BUCKET = bucket_dir_path;
global AFFINE_DIR = affine_dir_path;
global DATASET_DIR = dataset_dir_path;
global WAFER_DIR_DICT = waferpaths2dict(joinpath(dataset_dir_path, wafer_filename));
global SESSION = parseRoughAlign(joinpath(dataset_dir_path, rough_align_filename));

export tile_size, block_size, search_r, min_r, mesh_length, mesh_coeff, match_coeff, eta_grad, eta_newton, show_plot, num_procs, ftol_grad, ftol_newton, num_tiles, num_rows, num_cols;

tile_size = 8000;
block_size = 60;
search_r = 80;
min_r = 0.80;
mesh_length = 140;
mesh_coeff = 1;
match_coeff = 10;
eta_grad = 0.01;
eta_newton = .5; 
show_plot = false;
num_procs = length(procs());
ftol_grad = 1/1000;
ftol_newton = 1/1000000;
num_tiles = 16;
num_rows = 4;
num_cols = 4;


end


