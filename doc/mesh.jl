"""
`CREATE_MESH` - Create mesh for image with global offset, given point spacing

  ```mesh = create_mesh(img, offset, dist)```

* img: 2D array representing image
* offset: 2-element array for location of img's [0,0] pixel in global space
* dist: positive number for distance between two neighboring nodes in mesh

Created triangular mesh forms equilateral triangles, so all interior nodes
are equidistant from immediate neighbors.

```
   *-----*
  / \   / \
 /   \ /   \
*-----*-----* 
 \   / \   /
  \ /   \ /
   *-----*  

All edges are dist
```

"""
function create_mesh(img, offset, dist)
  sz = collect(size(img))
  spacing = [dist*sin(pi/3), dist]
  dims = round(Int64, sz./spacing + 1)
  margin = sz.%spacing / 2
  
  n = max(ij2index(dims, dims...), ij2index(dims, dims[1], dims[2]-1))
  m = 0
  nodes = Points(n)
  edges = spzeros(Float64, n, 3*n)

  for i in 1:dims[1], j in 1:dims[2]
    k = ij2index(dims, i, j)
    if k == 0 
      continue
    end
    nodes[k] = ij2coordinates(dims, i, j, offset+margin, spacing)
    if (j != 1)
      m += 1; edges[k, m] = -1; edges[ij2index(dims, i, j-1), m] = 1;
    end

    if (i != 1)
      if iseven(i) || j != dims[2]
        m += 1; edges[k, m] = -1; edges[ij2index(dims, i-1, j), m] = 1;
      end
      if iseven(i) && (j != dims[2])      
        m += 1; edges[k, m] = -1; edges[ij2index(dims, i-1, j+1), m] = 1;
      end
      if isodd(i) && (j != 1)
        m += 1; edges[k, m] = -1; edges[ij2index(dims, i-1, j-1), m] = 1;
      end
      if isodd(i) && ((j == 1) || (j == dims[2]))
        m += 1; edges[k, m] = -1; edges[ij2index(dims, i-2, j), m] = 1;
      end
    end
  end

  edges = edges[:, 1:m];
  return Mesh(offset, nodes, nodes, edges)
end

"""
Convert mesh row & column to node index
"""
function ij2index(dims, i, j)
  ind = 0
  if iseven(i) && (j == dims[2]) 
    return ind
  end
  if ((i < 1) || (j < 1) || (i > dims[1]) || (j > dims[2])) 
    return ind
  end
  ind += div(i-1, 2) * (dims[2] - 1)  #even rows
  ind += div(i, 2) * dims[2]          #odd rows
  ind += j
  ind = convert(Int64, ind)
  return ind
end

"""
Convert mesh row & column to coordinate in global space
"""
function ij2coordinates(dims, i, j, global_offset, spacing)
  if iseven(i) && (j == dims[2]) 
    return [0, 0]
  end
  pi = (i-1)*spacing[1] + global_offset[1]
  if iseven(i)  
    pj = (j-0.5)*spacing[2] + global_offset[2]
  else    
    pj = (j-1)*spacing[2] + global_offset[2]
  end
  return [pi, pj]
end
