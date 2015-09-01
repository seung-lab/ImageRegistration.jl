type Params
  	scaling_factor::Float64

  	# Mesh parameters
	mesh_length::Int64
	mesh_coeff::Float64

	# Matching parameters
	min_dyn_range_ratio::Float64
	block_size::Int64
	search_r::Int64
	min_r::Float64

	# Solver parameters
	match_coeff::Float64
	eta_gradient::Float64
	ftol_gradient::Float64
	eta_newton::Float64
	ftol_newton::Float64
end

function print_params(p::Params)
  	println("#Params:#######################")
	println("    scaling_factor: $(p.scaling_factor)");
	println("- Mesh parameters:");
	println("    mesh_length: $(p.mesh_length)");
	println("    mesh_coeff: $(p.mesh_coeff)");
	println("- Match parameters:");
	println("    min_dyn_range_ratio: $(p.min_dyn_range_ratio)");
	println("    block_size: $(p.block_size)");
	println("    search_r: $(p.search_r)");
	println("    min_r: $(p.min_r)");
	println("- Solver parameters:");
	println("    match_coeff: $(p.match_coeff)");
	println("    eta_gradient: $(p.eta_gradient)");
	println("    ftol_gradient: $(p.ftol_gradient)");
	println("    eta_newton: $(p.eta_newton)");
	println("    ftol_newton: $(p.ftol_newton)");
  	println("###############################")
end

SCALING_FACTOR_MONTAGE = 1.0;
MESH_LENGTH_MONTAGE = 200;
MESH_COEFF_MONTAGE = 1.0;
MIN_DYN_RANGE_RATIO_MONTAGE = 5;
BLOCK_SIZE_MONTAGE = 40;
SEARCH_R_MONTAGE = 80;
MIN_R_MONTAGE = 0.75;
MATCH_COEFF_MONTAGE = 20;
ETA_GRADIENT_MONTAGE = 0.01;
FTOL_GRADIENT_MONTAGE = 1/5000;
ETA_NEWTON_MONTAGE = 0.75;
FTOL_NEWTON_MONTAGE = 1/1000000;

SCALING_FACTOR_ALIGNMENT = 1.0;
MESH_LENGTH_ALIGNMENT = 1500;
MESH_COEFF_ALIGNMENT = 1.0;
MIN_DYN_RANGE_RATIO_ALIGNMENT = 5;
BLOCK_SIZE_ALIGNMENT = 1000;
SEARCH_R_ALIGNMENT = 500;
MIN_R_ALIGNMENT = 0.20;
MATCH_COEFF_ALIGNMENT = 20;
ETA_GRADIENT_ALIGNMENT = 0.01;
FTOL_GRADIENT_ALIGNMENT = 1/5000;
ETA_NEWTON_ALIGNMENT = 0.75;
FTOL_NEWTON_ALIGNMENT = 1/1000000;

global PARAMS_MONTAGE = Params(SCALING_FACTOR_MONTAGE, MESH_LENGTH_MONTAGE, MESH_COEFF_MONTAGE, MIN_DYN_RANGE_RATIO_MONTAGE, BLOCK_SIZE_MONTAGE, SEARCH_R_MONTAGE, MIN_R_MONTAGE, MATCH_COEFF_MONTAGE, ETA_GRADIENT_MONTAGE, FTOL_GRADIENT_MONTAGE, ETA_NEWTON_MONTAGE, FTOL_NEWTON_MONTAGE);

global PARAMS_ALIGNMENT = Params(SCALING_FACTOR_ALIGNMENT, MESH_LENGTH_ALIGNMENT, MESH_COEFF_ALIGNMENT, MIN_DYN_RANGE_RATIO_ALIGNMENT, BLOCK_SIZE_ALIGNMENT, SEARCH_R_ALIGNMENT, MIN_R_ALIGNMENT, MATCH_COEFF_ALIGNMENT, ETA_GRADIENT_ALIGNMENT, FTOL_GRADIENT_ALIGNMENT, ETA_NEWTON_ALIGNMENT, FTOL_NEWTON_ALIGNMENT);

function optimize_all_cores(params::Params)
  	img_d = 2 * (params.search_r + params.block_size) + 1;
	optimize_all_cores(img_d);
end
