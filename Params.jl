module Params

export BUCKET, AFFINE_DIR, WAFER_DIR

BUCKET = ".";
AFFINE_DIR = "~/150502_piriform/affine_transforms";
WAFER_DIR = "~/seungmount/research/GABA/data/atlas/MasterUTSLdirectory/07122012S2/S2-W001/HighResImages_ROI1_7nm_120apa"; #waferpath2dict("")

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



