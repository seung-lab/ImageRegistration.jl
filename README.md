[![Build Status](https://travis-ci.org/seung-lab/ImageRegistration.svg?branch=master)](https://travis-ci.org/seung-lab/ImageRegistration)

# ImageRegistration.jl
An image registration package toolbox for Julia. 

* Create point set correspondences
* Calculate geometric transforms
* Render images with transforms

## Installation
Use the package manager:
```
Pkg.add("ImageRegistration")
```

## Dependencies
* FixedPointNumbers (to allow for Ufixed series of image types)
* Images 
* ImageView (for the visualization functions)

### Types
There are two basic types to derive and apply transforms in this pacakge, the **Mesh** type and the **Matches** type. 

The Mesh type includes two sets of nodes, source nodes (src_nodes) and the destination nodes that represent the transformation of the mesh (dst_nodes), as well as an incidence matrix the represents the edges of the mesh. See the `Mesh` documentation for more information.

The Matches type includes two set of points (keeping with mathematical terminology, though they are identical to nodes in implementation), as well as a list of indices representing the nodes from a given mesh with which the matches correspond (mesh_indices).

Mesh types exist to apply piecewise transforms. Matches types exist as a more general set of correspondences. Viewed another way, Mesh types apply to only one image, while Matches indicate that there is a correspondence between two images. In the near future, Matches will link Meshes together to enable the registration of systems with more than one mesh.

### Create point sets
The blockmatch method takes two images along with some parameters. It creates a mesh so that at each node it can blockmatch between one image and the other. The blockmatched nodes which pass the filter criteria (for now, just a threshold on the cross-correlation value), will create a Matches type. The initial mesh and the matches type are returned.

### Geometric transforms

* **affine** transforms
* **rigid** transforms (isometries)
* **translations** (displacements)
* **piecewise affine** transforms (diffeomorphisms)
 
All transforms are calculated and applied as right-hand matrix mulitplication. See the documentation in `IMWARP` and `MESHWARP`, as well as each of the transforms for more information.

### Render

* **imwarp**: apply a non-piecewise transform to an entire image
* **meshwarp**: apply a piecewise transform, as defined by a Mesh, to an image

### Indexing notation
Since Julia stores arrays in column-majored order, this package uses row-column indexing (row represented by "i", and column represented by "j"). For example, the coordinate [5,2] will access the element in the 5th row and 2nd column of an array.

### Offsets
Every image has an offset to translate between its intrinsic coordinate system and a global coordinate system that it might share with other images. The offset is a 2-element array (again, in i,j notation) that indicates where the upper left pixel in the image is located within the global space.

When a transform is applied to an image, the render method will return the transformed image along with a new offset, denoting the transformed image's new location in global space.

## Partial list of core methods
```
matches = blockmatch(imgA, imgB, offsetA, offsetB, params)
img, offset = imwarp(img, tform, offset)
img, offset = meshwarp(img, mesh, offset)
tform = calculate_translation(matches)
tform = calculate_rigid(matches)
tform = calculate_affine(matches)
tform = calculate_translation(mesh)
tform = calculate_rigid(mesh)
tform = calculate_affine(mesh)
new_mesh = matches2mesh(matches, old_mesh)
```

### Not currently included in this package
* Alpha masks
* Feature-based image registration (SURF or SIFT)
* Multi-mesh solvers (i.e. affine solvers for more than one mesh, elastic solvers)
* imfuse visualization function (see Overlay type in Images for a start)

### Examples
```
"""
Create a mesh with src_nodes and deformed dst_nodes, then use
meshwarp to deform the image via piecewise affine transforms.
"""
img = load_test_image() # load test/test_images/turtle.jpg
src_nodes = [20.0 20.0;
                620.0 20.0;
                620.0 560.0;
                20.0 560.0;
                320.0 290.0]
dst_nodes = [20.0 20.0;
                620.0 20.0;
                620.0 560.0;
                20.0 560.0;
                400.0 460.0]
edges = spzeros(Int64, 8, 5)
edges[1,1:2] = [1, -1]
edges[2,1] = 1
edges[2,4] = -1
edges[3,1] = 1
edges[3,5] = -1
edges[4,2:3] = [1, -1]
edges[5,2] = 1
edges[5,5] = -1
edges[6,3] = 1
edges[6,5] = -1
edges[7,3:4] = [1, -1]
edges[8,4:5] = [1, -1]
# 1  -1   0   0   0
# 1   0   0  -1   0
# 1   0   0   0  -1
# 0   1  -1   0   0
# 0   1   0   0  -1
# 0   0   1   0  -1
# 0   0   1  -1   0
# 0   0   0   1  -1

mesh = Mesh(src_nodes, dst_nodes, edges)
imgc, img2 = view(img, pixelspacing=[1,1])
draw_mesh(imgc, img2, mesh)

warped_img, offset = meshwarp(img, mesh)
wimgc, wimg2 = view(warped_img, pixelspacing=[1,1])
draw_mesh(wimgc, wimg2, mesh)
```
