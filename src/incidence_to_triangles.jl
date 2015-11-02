function find_zero_indices(a)
  eachindex(a)'[a .== 0]
end

function find_nonzero_indices(a)
  eachindex(a)'[a .!= 0]
end

"""
`INCIDENCE_TO_DICT` - Create dictionary of node sets connected to indexing node
"""
function incidence_to_dict(D)
  D = abs(D)
  node_dict = Dict()

  for i = 1:size(D,1)
    j = 1
    while D[i,j] == 0
      j += 1
    end
    if !(j in keys(node_dict))
      node_dict[j] = Set{Int64}()
    end
    k = j+1
    while D[i,k] == 0
      k += 1
    end
    push!(node_dict[j], k)
  end
  return node_dict
end

"""
`DICT_TO_TRIANGLES` - Convert node dict to Nx3 list triangle node index vectors
"""
function dict_to_triangles(node_dict)
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
`INCIDENCE_TO_TRIANGLES` - Convert edge-node incidence matrix into triangle 
array (triples of node indices defining all triangles in the mesh). The triangle 
array will be used when applying piecewise affine transforms through the 
meshwarp function.

See Mesh type documentation for definition of an incidence matrix.

```
triangles = incidence_to_triangles(D)
```

* D: edge-node incidence matrix
* triangles: Nx3 array of ints, each row represents the indices of three nodes
    in one triangle.
"""
function incidence_to_triangles(D)
  node_dict = incidence_to_dict(D)
  return dict_to_triangles(node_dict)
end
