using HDF5
using Images
using ImageView
using PiecewiseAffineTransforms
using Color
using FixedPointNumbers

include("incidence2triangles")

function xy2ij(shape, height)
	[height .- shape[:, 2] shape[:, 1]]
end

function demo_mesh()
	bucket = "/usr/people/tmacrina/seungmount/research/tommy/Julimaps"
	file = "11x11.h5"
	path = joinpath(bucket, file)
	v = h5read(path, "v") # nodes
	e = h5read(path, "e") # edges
	k = h5read(path, "k") # spring constants
	l = h5read(path, "l") # resting lengths
	n = h5read(path, "n") # number of edges
	m = h5read(path, "m") # number of nodes

	display_mesh(e', v')
end

function display_mesh(incidence, nodes)
	# Display white mesh on black background
	triangles = incidence2triangles(incidence)
	sz = round(Int,maximum(nodes[:,1]))+1, round(Int,maximum(nodes[:,2]))+1
	img = Image(zeros(sz)) # Needs to be Float
	shape = xy2ij(nodes, size(img, 1))

	println("Original images and shapes")
	triplot(img, shape, triangles)
end

function warp_piecewiseaffine(img, incidence, initial_nodes, final_nodes)
	triangles = incidence2triangles(incidence)
	return pa_warp(img, initial_nodes, final_nodes, triangles)
end