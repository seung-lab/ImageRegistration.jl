[![Build Status](https://travis-ci.org/seung-lab/ImageRegistration.jl.svg?branch=master)](https://travis-ci.org/seung-lab/ImageRegistration.jl)

# ImageRegistration.jl
An image registration toolbox for Julia. 

* Create point set correspondences
* Calculate geometric transforms
* Render images with transforms

## Installation
Install via the package manager:
```
Pkg.add("ImageRegistration")
```

## Dependencies
* Images 
* Cairo (for visualizations)

## The Toolbox
### Types
There are two basic types to derive and apply transforms in this pacakge, the **Mesh** type and the **Matches** type. 

The Mesh type includes two sets of nodes, source nodes (src_nodes) and the destination nodes that represent the transformation of the mesh (dst_nodes), as well as an incidence matrix the represents the edges of the mesh. See the `Mesh` documentation for more information.

The Matches type includes two set of points (keeping with mathematical terminology, though they are identical to nodes in implementation), as well as a list of indices representing the nodes from a given mesh with which the matches correspond (mesh_indices).

Mesh types exist to apply piecewise transforms. Matches types exist as a more general set of correspondences. Viewed another way, Mesh types apply to only one image, while Matches indicate that there is a correspondence between two images. In the near future, Matches will link Meshes together to enable the registration of systems with more than one mesh.

### Create point sets
The blockmatch method takes two images along with some parameters. It creates a mesh so that at each node it can blockmatch between one image and the other. The blockmatched nodes which pass the filter criteria (for now, just a threshold on the cross-correlation value), will create a Matches type. The initial mesh and the matches type are returned.

Blockmatch is written to handle parallel computation. Starting Julia across multiple nodes will automatically take advantage of it:
```
julia -p <number of nodes>
```

### Geometric transforms

* **affine** transforms
* **rigid** transforms (isometries)
* **translations** (displacements)
* **piecewise affine** transforms (diffeomorphisms)
 
All transforms are calculated and applied as right-hand matrix mulitplication. See the documentation in `IMWARP` and `MESHWARP`, as well as each of the transforms for more information.

### Render

* **imwarp**: apply a non-piecewise transform to an entire image
* **meshwarp**: apply a piecewise transform, as defined by a Mesh, to an image

## Additional information
### Work with arrays only
For now, the blockmatch function only works on arrays, not on Image types. So if you read in an image using Images' imread, make sure to convert that type to just the data array.

### Indexing notation
Since Julia stores arrays in column-majored order, this package uses row-column indexing (row represented by "i", and column represented by "j"). For example, the coordinate [5,2] will access the element in the 5th row and 2nd column of an array.

### Offsets
Every image has an offset to translate between its intrinsic coordinate system and a global coordinate system that it might share with other images. The offset is a 2-element array (again, in i,j notation) that indicates where the upper left pixel in the image is located within the global space.

When a transform is applied to an image, the render method will return the transformed image along with a new offset, denoting the transformed image's new location in global space.

### Partial list of core methods
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
new_mesh = matches_to_mesh(matches, old_mesh)
```
### Visualizations
Working with the ImageView package, there are some functions included to help visualize function outputs.

* **draw_vectors**: display the starting point and line segment of correspondence points (matches or mesh) on an image
* **draw_mesh**: display the nodes and their edges on an image

### Not currently included in this package
* Alpha masks
* Feature-based image registration (SURF or SIFT)
* Multi-mesh solvers (i.e. affine solvers for more than one mesh, elastic solvers)
* imfuse visualization function (see Overlay type in Images for a start)

### Example
```
# Load the package and load two images
using ImageRegistration
imgA, imgB = load_test_images() # Loading two images from test/test_images

# Set the parameters for the blockmatch
params = default_params()
params["search_r"] = 1000

# Blockmatch the first image to the second imag (both with offsets of [0,0])
mesh, matches = blockmatch(imgA, imgB, [0,0], [0,0], params)

# Calculate rigid transform and render the first image with it
tform = calculate_rigid(matches)
rigid_imgA, rigid_offset = imwarp(imgA, tform)

# Convert matches to mesh & render piecewise affine transform of first image
warped_mesh = matches_to_mesh(matches, mesh)
warped_imgA, warped_offset = meshwarp(imgA, warped_mesh)

# Visually compare the images
using ImageView
using FixedPointNumbers
view(convert(Array{ufixed8}, imgA), pixelspacing=[1,1])
view(convert(Array{ufixed8}, imgB), pixelspacing=[1,1])
view(convert(Array{ufixed8}, rigid_imgA), pixelspacing=[1,1])
view(convert(Array{ufixed8}, warped_imgA), pixelspacing=[1,1])
```
