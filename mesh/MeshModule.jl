module MeshModule

importall Julimaps
importall Params
importall IO

export montageSection, montageSections

using Images
using HDF5
using JLD
include("../convolve.jl")
include("Mesh.jl")
include("Matches.jl")
include("MeshSet.jl")
include("MeshSolve.jl")

function montageSection(n)
	@time Ms, imageArray = MeshModule.loadSection(SESSION, n);
	@time addAllMatches!(Ms, imageArray);
	@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
	printResidualStats(Ms);
	save(Ms);
end

function montageSections(k::UnitRange{Int64})
	imageArray = SharedArray(Float64, tile_size, tile_size, num_tiles);
	for i in k
		@time Ms, imageArray = loadSection(SESSION, i);
		@time addAllMatches!(Ms, imageArray);
		@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
		printResidualStats(Ms);
		save(Ms);
	end
end




################################### DIAGNOSTIC #########################################


end #END MESHMODULE
################################### SCRIPT ############################################

