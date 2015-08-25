module MeshModule

importall Julimaps
importall Params
importall IO

export montageSection, montageSections

using Images
using HDF5
using JLD
include("./convolve.jl")
include("Mesh.jl")
include("Matches.jl")
include("MeshSet.jl")
include("MeshSolve.jl")

function montageSection(n)
	@time Ms, imageArray = loadSection(SESSION, n);
	@time addAllMatches!(Ms, imageArray);
	@time solveMeshSet!(Ms, match_coeff, eta_grad, eta_newton, ftol_grad, ftol_newton);
	printResidualStats(Ms);
	save(Ms);
	imageArray = 0;
	gc();
end

function montageSections(k::UnitRange{Int64})
	for i in k
		@time montageSection(i);
	end
end




################################### DIAGNOSTIC #########################################


end #END MESHMODULE
################################### SCRIPT ############################################

