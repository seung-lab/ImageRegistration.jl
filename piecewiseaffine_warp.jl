function source_point_mat(X, Y)
    x1, x2, x3 = X
    y1, y2, y3 = Y    
    F = [x1 x2 x3; y1 y2 y3; 1 1 1]
    F
end

function target_point_mat(U, V)
    u1, u2, u3 = U
    v1, v2, v3 = V
    T = [u1 u2 u3; v1 v2 v3]
    T
end

function affine_params(X, Y, U, V)
    F = source_point_mat(X, Y)
    T = target_point_mat(U, V)
    M = T * inv(F)
    M
end

warp_pixel(M, x, y) = M * [x, y, 1]

function pa_warp2{N}(img::Array{Float64, N},
                    src::Matrix{Float64}, dst::Matrix{Float64},
                    trigs::Matrix{Int64}, interp = true)
    warped = zeros(eltype(img), size(img))
    if interp
        println("w interpolation")
    else
        println("w/o interpolation")
    end
    for t=1:size(trigs, 1)    
        tr = squeeze(trigs[t, :], 1)
        Y = src[tr, 1]
        X = src[tr, 2]
        
        V = dst[tr, 1]
        U = dst[tr, 2]
               
        # warp parameters from target (U, V) to source (X, Y)
        M = affine_params(U, V, X, Y)
        
        # mask = poly2mask2(U, V, size(img)[1:2]...)
        # vs, us = findn(mask)
        vs, us = poly2source(U, V)
        
        # for every pixel in target triangle we find corresponding pixel in source
        # and copy its value
        for i=1:length(vs)
            u, v = us[i], vs[i]
            x, y = M * [u, v, 1]
            if interp
                fx, fy = floor(Int64, x), floor(Int64, y)
                wx, wy = x-fx, y-fy
                if 1 <= y && y+1 <= size(img, 1) && 1 <= x && x+1 <= size(img, 2)
                    p = [1-wy wy] * img[fy:fy+1, fx:fx+1] * [1-wx; wx]
                    warped[v, u] = p[1]
                end
            else
                y = round(Int64,y)
                x = round(Int64,x)
                if 1 <= y && y <= size(img, 1) && 1 <= x && x <= size(img, 2)
                    warped[v, u] = img[y, x]
                end
            end
        end
        
    end
    warped
end

function poly2source(px, py)
# Run fillpoly2 and findn over the basic bounding box of a given triangle
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    top, bottom = floor(Int64,minimum(py)), ceil(Int64,maximum(py))
    mask = zeros(Bool, bottom-top+1, right-left+1)
    fill = fillpoly2!(mask, px-left+1, py-top+1, true)
    vs, us = findn(fill)
    vs += top-1
    us += left-1
    return vs, us
end


function fillpoly2!{T,P<:Number}(M::Matrix{T}, px::Vector{P}, py::Vector{P}, value::T)
    @assert length(px) == length(py)    
    left, right = floor(Int64,minimum(px)), ceil(Int64,maximum(px))
    for x=left:right     
        ys = Set{Int64}()
        j = length(px)
        for i=1:length(px)            
            if (px[i] <= x && x <= px[j]) || (px[j] <= x && x <= px[i])
                # special case: adding the whole cut to ys                            
                if px[i] == px[j]
                    push!(ys, py[i])
                    push!(ys, py[j])
                else
                    y = py[i] + (x - px[i]) / (px[j] - px[i]) * (py[j] - py[i])
                    push!(ys, ceil(Int64, y))
                end            
            end
            j = i
        end
        ys = sort([y for y in ys])
        # if there's an odd number of intersection points, add one imeginary point
        if length(ys) % 2 == 1
            push!(ys, ys[end])
        end
        for i=1:2:length(ys)           
            M[ys[i]:ys[i+1], x] = value  # <-- bounds error here!
        end
    end
    return M
end


function poly2mask2{P}(px::Vector{P}, py::Vector{P}, m::Int64, n::Int64)
    mask = zeros(Bool, m, n)
    fillpoly2!(mask, px, py, true)
end