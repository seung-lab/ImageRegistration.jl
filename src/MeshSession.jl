
function montage_section(n)
  @time Ms, images = load_section(PREMONTAGED_OFFSETS, n);
  @time add_all_matches!(Ms, images);
  @time solve_meshset!(Ms);
  save(Ms);
  images = 0;
end

function montage_sections(wafer_num, k::UnitRange{Int64})
  optimize_all_cores(PARAMS_MONTAGE);
  for n in k
    @time montage_section(n);
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
