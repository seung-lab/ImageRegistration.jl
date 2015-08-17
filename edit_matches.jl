using Images
using ImageView

include("render.jl")

"""
Find index of location in array with point closest to the given point (if there
exists such a unique point within a certain pixel limit).

Args:

* pts: 2xN array of coordinates
* pt: 2-element coordinate array of point in interest
* limit: number of pixels that nearest must be within

Returns:

* index of the nearest point in the pts array

  idx = find_idx_of_nearest_pt(pts, pt, limit)
"""
function find_idx_of_nearest_pt(pts, pt, limit)
    d = sum((pts.-pt).^2, 1).^(1/2)
    idx = eachindex(d)'[d .< limit]
    if length(idx) == 1
        return idx[1]
    else
        return 0
    end
end

"""
Provide bindings for GUI to right-click and remove matches from a mesh. End 
manual removal by exiting the image window, or by pressing enter while focused
on the image window.

Args:

* imgc: ImageCanvas object (from image with annotations)
* img2: ImageZoom object (from image with annotations)
* annotation: the annotation object from the image

Returns:

* pts_to_remove: array of indices for points to remove from mesh

  pts_to_remove = edit_matches(imgc, img2, annotation)
"""
function edit_matches(imgc, img2, annotation)
    e = Condition()

    pts_to_remove = Array{Integer,1}()
    pts = copy(annotation.ann.data.pts)

    c = canvas(imgc)
    win = Tk.toplevel(c)
    bind(c, "<Button-3>", (c, x, y)->right_click(parse(Int, x), parse(Int, y)))
    bind(win, "<Return>", path->end_edit())
    bind(win, "<KP_Enter>", path->end_edit())
    bind(win, "<Destroy>", path->end_edit())

    function right_click(x, y)
        xu,yu = Graphics.device_to_user(Graphics.getgc(imgc.c), x, y)
        xi, yi = floor(Integer, 1+xu), floor(Integer, 1+yu)
        # println(xi, ", ", yi)
        limit = (img2.zoombb.xmax - img2.zoombb.xmin) * 0.0125 # 100/8000
        annpts = annotation.ann.data.pts
        annidx = find_idx_of_nearest_pt(annpts, [xi, yi], limit)
        if annidx > 0
            idx = find_idx_of_nearest_pt(pts, [xi, yi], limit)
            println(idx, ": ", [xi, yi])
            annotation.ann.data.pts = hcat(annpts[:,1:annidx-1], 
                                                        annpts[:,annidx+1:end])
            ImageView.redraw(imgc)
            push!(pts_to_remove, idx)
        end
    end

    function end_edit()
        println("End edit")
        notify(e)
        bind(c, "<Button-3>", path->path)
        bind(win, "<Return>", path->path)
        bind(win, "<KP_Enter>", path->path)
        bind(win, "<Destroy>", path->path)
    end

    println("Right click to remove correspondences, then press enter.")
    wait(e)

    return pts_to_remove
end

function demo_edit_matches()
  mesh_set = load(joinpath(BUCKET, "EM_images", "Test.jld"))["MeshSet"]
  src_pts, dst_pts = load_matches(mesh_set.matches[1])
  tile_pathB = joinpath(BUCKET, "EM_images", "Tile_r4-c3_S2-W001_sec20.tif")
  # draw_points(make_isotropic(rawdata(imread(tile_pathA))), src_pts)
  imgB = make_isotropic(rawdata(imread(tile_pathB)))
  imgc, img2, annotation = draw_points(imgB, dst_pts)
  a = edit_matches(imgc, img2, annotation)
end

"""
Remove matches from mesh_set corresponding to the indices provided.

Args:

* indices_to_remove: array of indices of match points to remove
* match_index: the index of the matches object in the mesh_set
* mesh_set: the mesh_set containing all the matches

Returns:

* (updated mesh_set)

  remove_matches_from_meshset!(indices_to_remove, match_index, mesh_set)
"""
function remove_matches_from_meshset!(indices_to_remove, match_index, mesh_set)
  matches = mesh_set.matches[match_index]
  flag = trues(matches.n)
  flag[indices_to_remove] = false
  matches.src_pointIndices = matches.src_pointIndices[flag]
  matches.dst_points = matches.dst_points[flag]
  matches.dst_triangles = matches.dst_triangles[flag]
  matches.dst_weights = matches.dst_weights[flag]
  matches.dispVectors = matches.dispVectors[flag]
  no_pts_removed = length(indices_to_remove)
  matches.n -= no_pts_removed
  mesh_set.m -= no_pts_removed
  mesh_set.m_e -= no_pts_removed
end

"""
Count number of displacement vector lengths k-sigma from the group mean

Args:

* vectors: 4xN array of points (start and end points of displacement vectors)
* k: number of sigmas to set the threshold for counting

Returns:

* integer for number of displacement vectors above the k-sigma threshold

  num = count_outliers(vectors, k)
"""
function count_outliers(vectors, k)
  d = sum((vectors[1:2,:] - vectors[3:4,:]).^2, 1).^(1/2)
  d_mean = mean(d)
  d_std = mean((d.-d_mean).^2).^(1/2)
  return sum((d.-d_mean)./d_std .> k)
end

"""
Loop through mesh_set & manually intervene at bad seams (3-sigma displacements)

Args:

* mesh_set: MeshSet with Matches objects

Returns:
* (updated mesh_set)

  fix_bad_seams(mesh_set)
"""
function fix_bad_seams(mesh_set)
  for (idx, match) in enumerate(mesh_set.matches)
    src_index =  match.src_mesh.index[3]
    dst_index =  match.dst_mesh.index[3]
    vectors = load_vectors(match)
    outliers = count_outliers(vectors, 3)
    println(idx, ": ", outliers)
    if outliers > 0
      println((src_index, dst_index))
      c = canvasgrid(1,2)
      tile_pathA = joinpath(BUCKET, match.src_mesh.path)
      tile_pathB = joinpath(BUCKET, match.dst_mesh.path)
      src_pts, dst_pts = load_matches(match)
      # vectors = load_vectors(match)
      imgA = make_isotropic(rawdata(imread(tile_pathA)))
      imgc, img2 = view(c[1,1], imgA)
      draw_points(imgc, img2, src_pts, RGB(0,0,1))
      imgB = make_isotropic(rawdata(imread(tile_pathB)))
      imgc, img2 = view(c[1,2], imgB)
      imgc, img2, points, lines = draw_vectors(imgc, img2, vectors, RGB(1,0,0))
      pts_to_remove = edit_matches(imgc, img2, points)
      println(pts_to_remove)
      if length(pts_to_remove) > 0
        remove_matches_from_meshset!(pts_to_remove, idx, mesh_set)
      end
    end
  end
end

"""
Run mesh_set back through elastic solver and save
"""
function rerun_mesh_set(mesh_set)
  match_coeff = 10;
  eta = 0.01;
  show_plot = false;
  grad_threshold = 1/1000;
  n_newton = 50;
  @time MeshModule.solveMeshSet!(mesh_set, match_coeff, 
                                                eta, grad_threshold, n_newton)
  @time MeshModule.MeshSet2JLD(joinpath(BUCKET, "EM_images", 
                              "section20x5_blocksize20_FIXED.jld"), mesh_set)
  return mesh_set
end

"""
Resave imfuse images of fixed seams
"""
function redisplay_mesh_set(mesh_set, problem_seams)
  for match in mesh_set.matches
    src_index =  match.src_mesh.index[3]
    dst_index =  match.dst_mesh.index[3]
    src_mesh =  match.src_mesh
    dst_mesh =  match.dst_mesh
    if (src_index, dst_index) in problem_seams || 
                          reverse((src_index, dst_index)) in problem_seams
      tile_pathA = joinpath(BUCKET, src_mesh.path)
      A, SR_A = warp_tile(src_mesh, tile_pathA) 

      dst_mesh =  match.dst_mesh
      tile_pathB = joinpath(BUCKET, dst_mesh.path)
      B, SR_B = warp_tile(dst_mesh, tile_pathB)

      O, SR_O = imfuse(A, SR_A, B, SR_B)
      imwrite(O, joinpath("/usr/people/tmacrina/Desktop/", 
                string(tile_pathA[15:end-4],"_",tile_pathB[15:end-4],".jpg") )) 
    end
  end
end

function fix_mesh_set()
  mesh_set = load(joinpath(BUCKET, "EM_images", "section20x5.jld"))["MeshSet"]
  fix_bad_seams(mesh_set)
  rerun_mesh_set(mesh_set)
end
