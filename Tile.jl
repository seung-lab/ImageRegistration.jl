Using MeshModule

type Tile
	img::Array{Uint8, 2}
	path::String
	affine::Array{Float64,2}
	mesh::Mesh
	id::Tuple{Int64, 3}
end

Tile(path::String, affine_path::String) = 
function load_image(tile::Tile)
	tile.img = imread(tile.path)
end

function load_section(dir_path)
# Returns:
# 	Array{Tile, 1}
end
