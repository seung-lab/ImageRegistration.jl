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
datasets_dir_path = "./datasets";
cur_dataset = "piriform";
montage_dir_path = "meshsets_montage";
alignment_dir_path = "meshsets_alignment";
wafer_filename = "wafer_paths.txt";
rough_align_filename = "rough_align.txt";

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, SESSION, MONTAGE_DIR, ALIGNMENT_DIR

global BUCKET = bucket_dir_path;
global AFFINE_DIR = affine_dir_path;
global DATASET_DIR = joinpath(datasets_dir_path, cur_dataset);
global MONTAGE_DIR = joinpath(datasets_dir_path, cur_dataset, montage_dir_path);
global ALIGNMENT_DIR = joinpath(datasets_dir_path, cur_dataset, alignment_dir_path);
global WAFER_DIR_DICT = waferpaths2dict(joinpath(datasets_dir_path, cur_dataset, wafer_filename));
global SESSION = parseRoughAlign(joinpath(datasets_dir_path, cur_dataset, rough_align_filename));

export tile_size, block_size, search_r, min_r, mesh_length, mesh_coeff, match_coeff, eta_grad, eta_newton, show_plot, num_procs, ftol_grad, ftol_newton, num_tiles, num_rows, num_cols, num_concurrent, num_procs_total;

tile_size = 8000;
block_size = 40;
search_r = 80;
min_r = 0.75;
mesh_length = 200;
block_size_alignment = 200;
search_r_alignment = 1200;
min_r_alignment = 0.25;
mesh_length_alignment = 4000;
mesh_coeff = 1;
match_coeff = 20;
eta_grad = 0.01;
eta_newton = .5; 
show_plot = false;
num_procs_total = maximum([length(procs())-1; 1]);
num_concurrent = 2;
num_procs = maximum([1, fld(num_procs_total, num_concurrent)]);
ftol_grad = 1/500;
ftol_newton = 1/1000000;
num_tiles = 16;
num_rows = 4;
num_cols = 4;


end


