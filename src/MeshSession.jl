
function montage_section(n)
	@time Ms, images = load_section(PREMONTAGED_OFFSETS, n);
	@time add_all_matches!(Ms, images);
	@time solve_meshset!(Ms);
	printResidualStats(Ms);
	save(Ms);
	imageArray = 0;
	gc();
end

function montage_sections(wafer_num, k::UnitRange{Int64})
	optimize_all_cores(PARAMS_MONTAGE);
	for n in k
	  montage_section(n);
	end

end

function align_stack(wafer_num, k::UnitRange{Int64})
	@time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, k);
	for i in 1:Ms.N-1
	  a = Ms.meshes[i].index[2];
	  b = Ms.meshes[i+1].index[2];
	 add_pair_matches!(Ms, a, b); 
	end
	solve_meshset!(Ms);
#	printResidualStats(Ms);
	save(Ms);
end
