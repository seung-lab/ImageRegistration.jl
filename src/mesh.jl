"""
`MESH` - Type defining a graph with nodes in two separate instances. 

* src_nodes: Nx2 array of points in global space (original instance)
* dst_nodes: Nx2 array of points in global space (updated instance)
* edges: NxM sparse incidence matrix, each column is zero except at two elements 
    with value -1 and 1
    i.e.
    edges    1  2  3  4
    nodes
        A    1  1  0  0  
        B   -1  0  1  1  
        C    0 -1 -1  0  
        D    0  0  0 -1  
        E    0  0  0  0  

```
            A - B   E
             \ / \ 
              C   D
```

"""
type Mesh
  src_nodes::Array 
  dst_nodes::Array
  edges::SparseMatrixCSC{Int64, Int64}
end

Mesh() = Mesh([], [], spzeros(0,0))

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
  dims = floor(Int64, sz./spacing + 1)
  margin = sz.%spacing / 2
  
  n = max(ij_to_index(dims, dims...), ij_to_index(dims, dims[1], dims[2]-1))
  m = 0
  nodes = zeros(Float64, (n,2))
  edges = spzeros(Float64, n, 3*n)

  for i in 1:dims[1], j in 1:dims[2]
    k = ij_to_index(dims, i, j)
    if k == 0 
      continue
    end
    nodes[k,:] = ij_to_coordinates(dims, i, j, offset+margin, spacing)
    if (j != 1)
      m += 1; edges[k, m] = -1; edges[ij_to_index(dims, i, j-1), m] = 1;
    end

    if (i != 1)
      if iseven(i) || j != dims[2]
        m += 1; edges[k, m] = -1; edges[ij_to_index(dims, i-1, j), m] = 1;
      end
      if iseven(i) && (j != dims[2])      
        m += 1; edges[k, m] = -1; edges[ij_to_index(dims, i-1, j+1), m] = 1;
      end
      if isodd(i) && (j != 1)
        m += 1; edges[k, m] = -1; edges[ij_to_index(dims, i-1, j-1), m] = 1;
      end
      if isodd(i) && ((j == 1) || (j == dims[2]))
        m += 1; edges[k, m] = -1; edges[ij_to_index(dims, i-2, j), m] = 1;
      end
    end
  end

  edges = edges[:, 1:m];
  return Mesh(nodes, nodes, edges)
end

"""
`IJ_TO_INDEX` - Convert mesh row & column to node index

```
ind = ij_to_index(dims, i, j)
```

* dims: tuple of dimension sizes in i & j
* i: row index
* j: column index
* ind: int representing index of the node in the mesh
"""
function ij_to_index(dims, i, j)
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
`IJ_TO_COORDINATES` - Convert mesh row & column to coordinate in global space

```
[pi, pj] = ij_to_coordinates(dims, i, j, offset, spacing)
```

* dims: tuple of dimension sizes in i & j
* i: row index
* j: column index
* offset: 2-element array for location of img's [0,0] pixel in global space
* spacing: 2-element array for distance between mesh rows and columns
* [pi, pj]: global coordinates of the mesh node
"""
function ij_to_coordinates(dims, i, j, offset, spacing)
  if iseven(i) && (j == dims[2]) 
    return [0, 0]
  end
  pi = (i-1)*spacing[1] + offset[1]
  if iseven(i)  
    pj = (j-0.5)*spacing[2] + offset[2]
  else    
    pj = (j-1)*spacing[2] + offset[2]
  end
  return [pi, pj]
end

"""
`GET_TRIANGLES` - Use mesh edge-node incidence matrix to form triangle array

```
triangles = get_triangles(mesh)
```

See `INCIDENCE_TO_TRIANGLES` documentation for more information.
"""
function get_triangles(mesh::Mesh)
  node_dict = incidence_to_dict(mesh.edges)
  return dict_to_triangles(node_dict)
end

function count_nodes(mesh::Mesh)
  return size(src_nodes, 1)
end