Using MeshModule

type Tile
	img::Array{UInt8, 2}
	path::String
	affine::Array{FloatingPoint,2}
	mesh::Mesh
	id::Tuple{Integer, 4} # wafer, section, row, col
end

Tile(path::String) = Tile(nothing, path, nothing, nothing, parsename(path))

function parsename(path)
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", path)
	return m[3], m[4], m[1], m[2]
end

function load_image(tile::Tile)
	tile.img = imread(tile.path)
end

function load_affine(tile::Tile)
	tile.affine = readcsv(tile.path)
end

function load_section(dir_path)
# Returns:
# 	Array{Tile, 1}
end
