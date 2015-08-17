using MeshModule
using SpatialRef

type Tile
	img::Array{UInt8, 2}
	mask::Array{Binary, 2}
	spatial_ref::SpatialRef
	name::String
	affine::Array{FloatingPoint,2}
	mesh::Mesh
	id::Tuple{Integer, 4} # wafer, section, row, col
end

Tile(path::String) = Tile(nothing, nothing, path, nothing, nothing, parsename(path))

const BUCKET = "."
const AFFINE_DIR = "~/seungmount/research/150502_piriform/affine_transforms"
const WAFER_DIR = waferpaths2dict("piriform_wafer_paths.txt")

"""
Create dictionary of waferpaths in bucket from file
"""
function waferpaths2dict(wafer_path_file)
	warray = readdlm(wafer_path_file)
	wdict = Dict()
	for (idx, path) in zip(warray[:,1], warray[:,2])
		wdict[idx] = path
	end
	return wdict
end

"""
Extract wafer no, section no, row, and col for tile from filename
"""
function parsename(path)
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", path)
	return m[3], m[4], m[1], m[2] # wafer, section, row, column
end

"""
Load original image for tile, using the WAFER_DIR and filename
"""
function load_image!(tile::Tile)
	section_folder = string("S2-W00", tile.id[1], "_Sec", tile.id[2], "_Montage")
	path = joinpath(WAFER_DIR[tile.id[1]], section_folder, string(tile.name, ".tif")
	tile.img = imread(path)
end

"""
Load affine transformation matrix for provided tile
"""
function load_affine!(tile::Tile)
	path = joinpath(AFFINE_DIR, string(tile.name, ".csv")
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
function set_mesh!(tile::Tile, mesh::Mesh})
	tile.mesh = mesh
end

"""
Remove image from a tile
"""
function remove_img!(tile::Tile)
	tile.img = []
	tile.spatial_ref = SpatialRef()
end

function load_section(dir_path)
# Returns:
# 	Array{Tile, 1}
end
