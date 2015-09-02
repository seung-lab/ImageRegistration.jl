function is_overview(index::Index)
    if index[3:4] == (OVERVIEW_INDEX, OVERVIEW_INDEX)   return true else return false end
end

function is_montaged(index::Index)
    if index[3:4] == (MONTAGED_INDEX, MONTAGED_INDEX)   return true else return false end
end

function is_prealigned(index::Index)
    if index[3:4] == (PREALIGNED_INDEX, PREALIGNED_INDEX)   return true else return false end
end

function is_aligned(index::Index)
    if index[3:4] == (ALIGNED_INDEX, ALIGNED_INDEX) return true else return false 
end
end


function parseName(name::String)

    ret = (0, 0, 0, 0)
    # singleton tile
    m = match(r"Tile_r(\d*)-c(\d*).*W(\d*)_sec(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[3]), parse(Int, m[4]), parse(Int, m[1]), parse(Int, m[2])
    end

    # overview image
    m = match(r"MontageOverviewImage_S2-W00(\d*)_sec(\d*)", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), OVERVIEW_INDEX, OVERVIEW_INDEX   
    end

    # montaged section
    m = match(r"(\d*),(\d*)_montaged", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), MONTAGED_INDEX, MONTAGED_INDEX   
    end

    # prealigned section
    m = match(r"(\d*),(\d*)_prealigned", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), PREALIGNED_INDEX, PREALIGNED_INDEX 
    end

    # aligned-section
    m = match(r"(\d*),(\d*)_aligned", name)
    if typeof(m) != Void
    ret = parse(Int, m[1]), parse(Int, m[2]), ALIGNED_INDEX, ALIGNED_INDEX 
    end

    return ret
    
end

function getName(index::Index)
    if is_overview(index)
    return string("MontageOverviewImage_S2-W00", index[1], "_sec", index[2])
    elseif is_montaged(index)
    return string(index[1], ",", index[2], "_montaged")
    elseif is_prealigned(index)
    return string(index[1], ",", index[2], "_prealigned")
    elseif is_aligned(index)
    return string(index[1], ",", index[2], "_aligned")
    else
    return string("Tile_r", index[3], "-c", index[4], "_S2-W00", index[1], "_sec", index[2])
    end
end

# function getPath()
# methods: 
#     
# extensions:
# Mesh.jl: getPath(mesh::Mesh)
function getPath(index::Index)
    name = getName(index)
    if is_overview(index)
        section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"))
    elseif is_montaged(index)
        path = joinpath(MONTAGED_DIR, string(name, ".tif"))
    elseif is_prealigned(index)
        path = joinpath(PREALIGNED_DIR, string(name, ".tif"))
    elseif is_aligned(index)
        path = joinpath(ALIGNED_DIR, string(name, ".tif"))
    else
        section_folder = string("S2-W00", index[1], "_Sec", index[2], "_Montage")
        path = joinpath(BUCKET, WAFER_DIR_DICT[index[1]], section_folder, string(name, ".tif"))
    end
    println(path)
    return path
end

function getPath(name::String)
    return getPath(parseName(name))
end



# extensions:
# Mesh.jl getImage(mesh::Mesh)
function getImage(path::String)
    img = imread(path)
    img = img[:, :, 1]
    img.properties["timedim"] = 0
    return convert(Array{UInt8, 2}, round(convert(Array, img)*255))
end

function getImage(index::Index)
    return getImage(getPath(index))
end


function waferpaths2dict(waferpath_filename)
    wdict = Dict()
    if isfile(waferpath_filename)
        warray = readdlm(waferpath_filename)
        for (idx, path) in zip(warray[:,1], warray[:,2])
            wdict[idx] = path
        end
    end
    return wdict
end

function parse_offsets(path::String)
    offsets = cell(0, 0)
    if isfile(path)
        file = readdlm(path)
        offsets = cell(size(file, 1), size(file, 2) + 1) # name, index, dx, dy
        for i in 1:size(offsets, 1)
            index = parseName(file[i, 1])
            offsets[i, 1] = getName(index)
            offsets[i, 2] = index
            for j in 3:size(offsets, 2)
                offsets[i, j] = file[i, j-1]
            end
        end
      end
    return offsets
end

bucket_dir_path = ""
if isfile("bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("bucket_dir_path.txt"), '\n')
elseif isfile("../bucket_dir_path.txt")
    bucket_dir_path = rstrip(readall("../bucket_dir_path.txt"), '\n')
end
datasets_dir_path = "research/Julimaps/datasets"
cur_dataset = "piriform"
affine_dir_path = "~"

premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"

wafer_filename = "wafer_paths.txt"
premontaged_offsets_filename = "premontaged_offsets.txt"
prealigned_offsets_filename = "prealigned_offsets.txt"

export BUCKET, DATASET_DIR, AFFINE_DIR, WAFER_DIR_DICT, PREMONTAGED_OFFSETS, PREMONTAGE_DIR, ALIGNMENT_DIR

global BUCKET = bucket_dir_path
global AFFINE_DIR = affine_dir_path
global DATASET_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset)
global PREMONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path)
global MONTAGED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, montaged_dir_path)
global PREALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path)
global ALIGNED_DIR = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, aligned_dir_path)

waferpath_filename = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, wafer_filename)
global WAFER_DIR_DICT = waferpaths2dict(waferpath_filename)
premontaged_offsets_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, premontaged_dir_path, premontaged_offsets_filename)
global PREMONTAGED_OFFSETS = parse_offsets(premontaged_offsets_path)
prealigned_offsets_path = joinpath(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path, prealigned_offsets_filename)
global PREALIGNED_OFFSETS = parse_offsets(prealigned_offsets_path)

show_plot = false
num_procs = nprocs()

