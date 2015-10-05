function find_zero_indices(a)
  eachindex(a)'[a' .== 0]
end

function find_nonzero_indices(a)
  eachindex(a)'[a' .!= 0]
end

"""
`INCIDENCE2DICT` - Create dictionary of sets of nodes connected to indexing node
"""
function incidence2dict(D)
  D = abs(D)
  node_dict = Dict()

  for j = 1:size(D,2)
    i = 1
    while D[i,j] == 0
      i += 1
    end
    if !(i in keys(node_dict))
      node_dict[i] = Set{Int64}()
    end
    k = i+1
    while D[k,j] == 0
      k += 1
    end
    push!(node_dict[i], k)
  end
  return node_dict
end

"""
`DICT2TRIANGLES` - Convert node dict to 1x3 list of triangle node index vectors
"""
function dict2triangles(node_dict)
  triangles = Array(Int64, 0, 3)
  for a in sort(collect(keys(node_dict)))
    setA = node_dict[a]
    for b in sort(collect(setA))
      if !(b in keys(node_dict))
        continue
      end
      setB = node_dict[b]
      setC = intersect(setA, setB)
      for c in sort(collect(setC))
        triangles = vcat(triangles, [a b c])
      end
    end
  end
  return triangles
end

"""
`INCIDENCE2TRIANGLES` - Convert edge-node incidence matrix into triangle array 
(triples of node indices defining all triangles in the mesh). The triangle array
will be used when applying piecewise affine transforms through the meshwarp 
function.

See Mesh type documentation for definition of an incidence matrix.

```
triangles = incidence2triangles(D)
```

* D: edge-node incidence matrix
* triangles: Nx3 array of ints, each row represents the indices of three nodes
    in one triangle.
"""
function incidence2triangles(D)
  node_dict = incidence2dict(D)
  return dict2triangles(node_dict)
end
