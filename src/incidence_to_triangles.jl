function find_zero_indices(a)
  eachindex(a)'[a .== 0]
end

function find_nonzero_indices(a)
  eachindex(a)'[a .!= 0]
end


"""
`INCIDENCE_TO_INCIDENT_LIST_ARRAY` - takes in M edges on N nodes in NxM sparse format and returns an 3xN array where the N-th column includes the indices of the nodes connected (one-direction) to node N
 * Currently assumes that there's at most 3 points, and furthermore that the nonzeros alternate
 * 0 in the incidence_array means that there's less than 3 points
"""
function incidence_to_incident_list_array(E)
  incident_list_array = zeros(Int64, 3, size(E, 1))
  incident_list_counts = ones(Int64, size(E, 1))
  rows = rowvals(E)
  vals = nonzeros(E)

  ind = 1
  # for each edge
  for j in 1:size(E, 2)
    	node_first = rows[ind]
    	node_second = rows[ind+1]
	incident_list_count = incident_list_counts[node_first]
	incident_list_array[incident_list_count, node_first] = node_second
	incident_list_counts[node_first] = incident_list_count+1
	ind = ind+2
  end
  return incident_list_array
end

# checks if p -> p(pa,pb) is the same as p -> pc -> pca(b,c) and writes the triangle
function check_and_write_triangles!(triangles, cur_triangle, p, pa, pb, pc, pca, pcb, pcc)
  if p == 0 return cur_triangle end
  if pc == 0 return cur_triangle end
  if pa != 0
	if pa == pca || pa == pcb || pa == pcc
	triangles[1, cur_triangle] = p
	triangles[2, cur_triangle] = pc
	triangles[3, cur_triangle] = pa
	cur_triangle += 1
      end
  end
  if pb != 0
	if pb == pca || pb == pcb || pb == pcc
	triangles[1, cur_triangle] = p
	triangles[2, cur_triangle] = pc
	triangles[3, cur_triangle] = pb
	cur_triangle += 1
      end
  end
  return cur_triangle
end
"""
`INCIDENT_LIST_ARRAY_TO_TRIANGLES` - Convert incidence array to 3xN list triangle node index vectors
"""
function incident_list_array_to_triangles(incident_list_array)
  triangles = zeros(Int64, 3, size(incident_list_array, 2) * 10)
  max_edges = size(incident_list_array, 1)
  cur_triangle = 1
  # absolutely awful, but really fast
  # for each node at j
  for p in 1:size(incident_list_array, 2)
    # the up to 3 points that have edges from j
    	pa = incident_list_array[1, p]
    	pb = incident_list_array[2, p]
    	pc = incident_list_array[3, p]

      if pa != 0
	# the up to 9 points that have edges from pa, pb, pc
	paa = incident_list_array[1, pa]; pab = incident_list_array[2, pa]; pac = incident_list_array[3, pa]
      else
	paa = 0; pab = 0; pac = 0
      end

      if pb != 0
	pba = incident_list_array[1, pb]; pbb = incident_list_array[2, pb]; pbc = incident_list_array[3, pb]
      else
	pba = 0; pbb = 0; pbc = 0
      end

      if pc != 0
	pca = incident_list_array[1, pc]; pcb = incident_list_array[2, pc]; pcc = incident_list_array[3, pc]
      else
	pca = 0; pcb = 0; pcc = 0
      end

	cur_triangle = check_and_write_triangles!(triangles, cur_triangle, p, pa, pb, pc, pca, pcb, pcc)
	cur_triangle = check_and_write_triangles!(triangles, cur_triangle, p, pa, pc, pb, pba, pbb, pbc)
	cur_triangle = check_and_write_triangles!(triangles, cur_triangle, p, pb, pc, pa, paa, pab, pac)
  end
  return triangles[:, 1:(cur_triangle - 1)]
end

"""
`INCIDENCE_TO_TRIANGLES` - Convert node-edge incidence matrix into triangle 
array (triples of node indices defining all triangles in the mesh). The triangle 
array will be used when applying piecewise affine transforms through the 
meshwarp function.

See Mesh type documentation for definition of an incidence matrix.

```
triangles = incidence_to_triangles(E)
```

* E: node-edge incidence matrix
* triangles: 3xN array of ints, each column represents the indices of three nodes
    in one triangle.
"""
function incidence_to_triangles(E)
	incident_list_array = incidence_to_incident_list_array(E)
	return incident_list_array_to_triangles(incident_list_array)
end

#= THESE FUNCTIONS ARE SLOW
"""
`INCIDENCE_TO_DICT` - Create dictionary of node sets connected to indexing node
"""
function incidence_to_dict(D)
  D = abs.(D)
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
  triangles = Array{Int64}(0, 3)
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
=#
