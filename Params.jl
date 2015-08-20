module Params

importall IO

bucket_dir_path = "~/seungmount";
affine_dir_path = "~";
dataset_dir_path = "./datasets/piriform"
wafer_path_file = "wafer_paths.txt";
info_path_file = "import_info.txt";

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR

global BUCKET = bucket_dir_path;
global AFFINE_DIR = affine_dir_path;
global DATASET_DIR = dataset_dir_path;
global WAFER_DIR = waferpaths2dict(wafer_path_file);

export block_size, search_r, min_r, mesh_length, mesh_coeff, match_coeff, eta_grad, eta_newton, show_plot, num_procs, ftol_grad, ftol_newton, num_tiles, num_rows, num_cols;

block_size = 50;
search_r = 120;
min_r = 0.80;
mesh_length = 100;
mesh_coeff = 1;
match_coeff = 100;
eta_grad = 0.001;
eta_newton = .5; 
show_plot = false;
num_procs = length(procs());
ftol_grad = 1/1000;
ftol_newton = 1/1000000;
num_tiles = 16;
num_rows = 4;
num_cols = 4;

end



