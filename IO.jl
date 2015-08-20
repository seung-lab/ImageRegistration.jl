module IO

export parseName, getPath, getFloatImage, loadSection, toJLD, parseRoughAlign, waferpaths2dict, loadSectionImages

using Julimaps
using Params
using Images

function parseName(name::String)
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", name)
	return parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2]) # wafer, section, row, column
end

# extensions:
# Mesh.jl: getPath(mesh::Mesh)
function getPath(name::String)
	index = parseName(name);
	section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage");
	return joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"));
end

# extensions:
# Mesh.jl getFloatImage(mesh::Mesh)
function getFloatImage(path::String)
	return convert(Array{Float64, 2}, convert(Array, imread(path)));
end

function loadAffine(path::String)
	affinePath = joinpath(AFFINE_DIR, string(path, ".csv"))
	return readcsv(path);
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
	session[i, 2] = m[1];
	session[i, 3] = parseName(array[i, 1]);
	session[i, 4] = file[i, 2];
	session[i, 5] = file[i, 3];
	end
	return session;
end

# extensions:
# MeshSet.jl loadSectionImages(Ms::MeshSet)
function loadSectionImages(session, section_num)
	indices = find(i -> session[i,2][2] == section_num, 1:size(session, 1))
	max_tile_size = 0;
	num_tiles = length(indices);
	paths = Array{String, 1}(num_tiles);
	for i in 1:num_tiles
		name = session[i, 1];
		paths[i] = getPath(name);
		image = getFloatImage(paths[i])
		max_size = max(size(image, 1), size(image, 2));
		if max_tile_size < max_size max_tile_size = max_size; end
	end
	imageArray = SharedArray(Float64, max_tile_size, max_tile_size, num_tiles);

	@time for k in 0:num_procs:num_tiles
		@sync @parallel for l in 1:num_procs
		i = k+l;
		if i > num_tiles return; end;
		image = getFloatImage(paths[i]);
		imageArray[1:size(image, 1), 1:size(image, 2), i] = image;
		end
	end

	return imageArray;

end

function toJLD()
	return;
end


end #module


#=
function 

function set_affine!(tile::Tile, affine::Array{Real,(3,3)})
	tile.affine = affine
end

"""
Extract wafer no, section no, row, and col for tile from filename
"""

"""
Load original image for tile, using the WAFER_DIR_DICT and filename
"""
function load_image(tile::Tile)
	section_folder = string("S2-W00", tile.id[1], "_Sec", tile.id[2], "_Montage")
	# path = joinpath(homedir(), WAFER_DIR_DICT[tile.id[1]], section_folder, string(tile.name, ".tif"))
	path = joinpath(".", "input_images", "W001_sec20", string(tile.name, ".tif"))
	return imread(path)
end
""" Load affine transformation matrix for provided tile """ function load_affine!(tile::Tile)
	path = joinpath(AFFINE_DIR, string(tile.name, ".csv"))
	tile.affine = readcsv(path)
end

"""
Set affine transformation matrix for tile
"""
function set_affine!(tile::Tile, affine::Array{Real,(3,3)})
	tile.affine = affine
end

"""
Set affine transformation matrix for tile
"""
function set_mesh!(tile::Tile, mesh::MeshModule.Mesh)
	tile.mesh = mesh
end

# """
# Remove image from a tile
# """
# function remove_img!(tile::Tile)
# 	tile.img = []
# 	tile.spatial_ref = SpatialRef()
# end

=#
