# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150810
#
# Create color plots from JLD meshes

using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "") # include current directory
@pyimport create_color_plots as ccp

include("piecewiseaffine_render.jl")

function create_hexplot(mesh_path, write_name, offset=[0,0])
    v, vt, e = load_mesh(mesh_path) # original verts, final verts, incidence
    v[1,:] = v[1,:] - offset[1]
    v[2,:] = v[2,:] - offset[2]
    vt[1,:] = vt[1,:] - offset[1]
    vt[2,:] = vt[2,:] - offset[2]
    ccp.write_image_from_points(v', vt', write_name)
end

function demo_color_plot()
    mesh_path = joinpath(BUCKET, "EM_Images", "r4c3_solved.jld")
    write_name = "Tile_r4-c3_S2-W001_sec20_hexplot2.png"
    offset = [29090, 36251]
    create_hexplot(mesh_path, write_name, offset)
end