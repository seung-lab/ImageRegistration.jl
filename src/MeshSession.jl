
function montage_section(wafer_num, n)
  @time Ms, images = load_section(PREMONTAGED_OFFSETS, wafer_num, n);
  @time add_all_matches!(Ms, images);
  @time solve_meshset!(Ms);
  save(Ms);
  images = 0;
end

function montage_sections(wafer_num, k::UnitRange{Int64})
  optimize_all_cores(PARAMS_MONTAGE);
  for n in k
    @time montage_section(wafer_num, n);
  end

end

function align_stack(wafer_num, k::UnitRange{Int64})
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, k);
  @time for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
    @time add_pair_matches!(Ms, a, b); 
  end
  solve_meshset!(Ms);
  save(Ms);
end

function align_stack(wafer_num, k::UnitRange{Int64}, fixed_interval)
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, k, fixed_interval);
  @time for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
    @time add_pair_matches!(Ms, a, b); 
  end
  solve_meshset!(Ms);
  save(Ms);
end

function align_batch_to_fixed(wafer_num, aligned, batch::UnitRange{Int64})
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, aligned, batch);
  @time for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
   @time add_pair_matches!(Ms, a, b); 
  end
  solve_meshset!(Ms);
  save(Ms);
end

function prealign(wafer_num, dst, src)# k::UnitRange{Int64})
  @time Ms = affine_make_stack(MONTAGED_OFFSETS, wafer_num, dst, src);
  @time affine_add_pair_matches!(Ms, src, dst);
  @time affine_solve_meshset!(Ms);
  save(Ms);
  return Ms;
end

function prealign(wafer_num_a, sec_num_a, wafer_num_b, sec_num_b)# k::UnitRange{Int64})
if wafer_num_a != wafer_num_b println("No support for different wafers yet."); return; end
for src in sec_num_a:sec_num_b
dst = src - 1;
  @time Ms = affine_make_stack(MONTAGED_OFFSETS, wafer_num_a, dst, src);
  @time affine_add_pair_matches!(Ms, src, dst);
  @time affine_solve_meshset!(Ms);
  save(Ms);
  end
end

function align_to_fixed(wafer_num, aligned, prealigned)
  @time Ms = make_stack(PREALIGNED_OFFSETS, wafer_num, aligned, prealigned);
  @time for i in 1:Ms.N-1
    @time a = Ms.meshes[i].index[2];
    @time b = Ms.meshes[i+1].index[2];
   @time add_pair_matches!(Ms, a, b); 
  end
  solve_meshset!(Ms);
  save(Ms);
end

function prealign(section_range::UnitRange{Int64}, is_batch_start=false)
  for k in section_range
    println("Prealignment for montaged image no. ", k)
    if k == 1
      # check that first image has been copied through
      first_img_fn = sort_dir(MONTAGED_DIR, "tif")[k]
      first_img = get_ufixed8_image(joinpath(MONTAGED_DIR, first_img_fn))
      @time first_img = rescopeimage(first_img, [0,0], GLOBAL_BB)
      index = parse_name(first_img_fn)

      # Save image to prealigned
      log_path = joinpath(PREALIGNED_DIR, "prealigned_offsets.txt")
      if !isfile(log_path)
        f = open(log_path, "w")
        close(f)
      end
      log_file = open(log_path, "a")
      new_fn = string(join(index[1:2], ","), "_prealigned.tif")
      println("Writing ", new_fn)
      @time imwrite(first_img, joinpath(PREALIGNED_DIR, new_fn))
      log_line = join((new_fn, 0, 0, 
                          size(first_img,1), size(first_img,2)), " ")
      write(log_file, log_line, "\n")
      close(log_file)

      # Save iamge to aligned
      log_path = joinpath(ALIGNED_DIR, "aligned_offsets.txt")
      if !isfile(log_path)
        f = open(log_path, "w")
        close(f)
      end
      log_file = open(log_path, "a")
      new_fn = string(join(index[1:2], ","), "_aligned.tif")
      println("Writing ", new_fn)
      @time imwrite(first_img, joinpath(ALIGNED_DIR, new_fn))
      log_line = join((new_fn, 0, 0, 
                          size(first_img,1), size(first_img,2)), " ")
      write(log_file, log_line, "\n")
      close(log_file)

    else
      if k == section_range[1] && is_batch_start # beginning of a batch
        moving_fn = sort_dir(MONTAGED_DIR, "tif")[k][1:end-4]
        fixed_fn = sort_dir(ALIGNED_DIR, "tif")[k-1][1:end-4]
      else
        moving_fn = sort_dir(MONTAGED_DIR, "tif")[k][1:end-4]
        fixed_fn = sort_dir(PREALIGNED_DIR, "tif")[k-1][1:end-4]
      end
      println("Prealigning ", moving_fn, " to ", fixed_fn)
      affine_align_sections(moving_fn, fixed_fn)
    end
  end
end

function premontage(section_range::UnitRange{Int64})
  for sec in section_range
    index = (1, sec, 0, 0)   # Wafer #1 for now
    overview_path = get_path(get_overview_index(index))
    dir,name = splitdir(overview_path)
    println(dir)
    tiles = sort_dir(dir, "tif");
    tiles = filter(x->contains(x,"Tile"), tiles)

    save_fused_img_to = name[1:end-4]"_fused.png"
    save_xcorr_img_to = name[1:end-4]"_xcorr.png"
    offsets, = tiles_to_overview(tiles, overview_path, 0.07; tile_img_dir = dir,
        save_fused_img_to = joinpath(PREMONTAGED_DIR, save_fused_img_to),
        save_xcorr_img_to = joinpath(PREMONTAGED_DIR, save_xcorr_img_to))

    offset_file = joinpath(PREMONTAGED_DIR, "premontaged_offsets_julimaps.txt")
    f = open(offset_file, "a")
    for pair in offsets
      line = join((pair[1], pair[2]...), " ")
      write(f, line, "\n")
    end
    close(f)
  end
end
