module IO

export parsePath, loadImage

using Julimaps
using Params
using Images

function parsePath(path::String)
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", path)
	return parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2]) # wafer, section, row, column
end

function loadAffine(path::String)
	affinePath = joinpath(AFFINE_DIR, string(path, ".csv"))
	return readcsv(path);
end

function set_affine!(tile::Tile, affine::Array{Real,(3,3)})
	tile.affine = affine
end
end #module
#=
function 

"""
Extract wafer no, section no, row, and col for tile from filename
"""

"""
Load original image for tile, using the WAFER_DIR and filename
"""
function load_image(tile::Tile)
	section_folder = string("S2-W00", tile.id[1], "_Sec", tile.id[2], "_Montage")
	# path = joinpath(homedir(), WAFER_DIR[tile.id[1]], section_folder, string(tile.name, ".tif"))
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
