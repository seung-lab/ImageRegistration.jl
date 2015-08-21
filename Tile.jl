using Julimaps
using Params

include("Mesh.jl")
using MeshModule


type Tile
	path::String
	affine
	mesh::Mesh
	id::Index # wafer, section, row, col
end

Tile(path::String) = Tile(path, nothing, nothing, parsepath(path))

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

const BUCKET = "."
const AFFINE_DIR = "~/seungmount/research/150502_piriform/affine_transforms"
const WAFER_DIR = waferpaths2dict("piriform_wafer_paths.txt")

"""
Extract wafer no, section no, row, and col for tile from filename
"""
function parsename(path)
	m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", path)
	return parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2]) # wafer, section, row, column
end

"""
Load original image for tile, using the WAFER_DIR and filename
"""
function load_image(tile::Tile)
	# section_folder = string("S2-W00", tile.id[1], "_Sec", tile.id[2], "_Montage")
	# path = joinpath(homedir(), WAFER_DIR[tile.id[1]], section_folder, string(tile.name, ".tif"))
	section_folder = string("W00", tile.id[1], "_Sec", tile.id[2])
	path = joinpath(".", "input_images", section_folder, string(tile.name, ".tif"))
	return rawdata(imread(path))
end

"""
Load affine transformation matrix for provided tile
"""
function load_affine!(tile::Tile)
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

"""
Parse Mesh object to retrieve nodes, edges, and spatial reference

Args:

* mesh: Mesh object

Returns:

* src_nodes: 2xN array of original mesh nodes
* dst_nodes: 2xN array of mesh nodes after elastic solving
* mesh.edges: MxN array of edges as edge-node incidence matrix 
* offset: global offset of the mesh
"""
function parse_mesh(mesh)
    src_nodes = hcat(mesh.nodes...)
    dst_nodes = hcat(mesh.nodes_t...)
    offset = convert(Array{Int64,1}, mesh.disp)

    src_nodes = xy2yx(src_nodes .- offset)
    dst_nodes = xy2yx(dst_nodes .- offset)
    return src_nodes, dst_nodes, mesh.edges, offset
end

function load_section(dir_path)
# Returns:
# 	Array{Tile, 1}
end

"""
Create tile array with meshes from mesh_set
"""
function load_tiles(mesh_set)
	tiles = []
	for mesh in mesh_set.meshes
		t = Tile(mesh.path[13:end-4])
		set_mesh!(t, mesh)
		push!(tiles, t)
	end
	return tiles
end
