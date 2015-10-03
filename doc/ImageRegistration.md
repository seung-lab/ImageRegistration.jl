# Registration & Transforms
An image registration and image transformation package for Julia.

## Installation
Use the package manager:

```
Pkg.add("ImageRegistration")
```

## Who is this package for?
This package was designed for people migrating pipelines from other languages who are in need of similar image registration functionality. It was produced in a neuroscience lab, initially for aligning electron microscopy images, then abstracted for more general applications.

## Dependencies
Though not necessary, this package works well with the Images and ImageView 
packages.

## Methods
* matches = blockmatch(imgA, imgB, offsetA, offsetB, params)

//* blockmatch(meshA, meshB, params)

* img, offset = imwarp(img, tform, offset)
* img, offset = meshwarp(img, src_mesh, dst_mesh, triangles, offset)
* tform = solve_translation(matches)
* tform = solve_rigid(matches)
* tform = solve_affine(matches)
* src_mesh, dst_mesh = matches2mesh(matches)

//* meshA, meshB = solve_elastic(meshA, meshB, matches, params)

//* convolve(imgA, imgB, params)

### Objects
#### Mesh
* img_path
* index
* offset
* src_nodes
* dst_nodes
* fixed_nodes
* edges

#### Matches
* src_index
* dst_index
* src_points
* dst_points

#### MeshSet
* Meshes
* Matches
* Params

### Registration
#### convolution
#### blockmatch

### Transformation Methods

Why not geometric transforms?

#### geometric
##### affine
Solve a system of moving and fixed points to preserve parallelism.
##### rigid
Solve a system of moving and fixed points to preserve distances and angles.
##### similar
Solve a system of moving and fixed points to preserve distance ratios.
#### elastic
Treat meshes and their connections between each other (as contained in a
meshset) as a spring system. Then deform those meshes to match their
correspondences by relaxing the spring system.

### Visualizations
#### imfuse
Create a two-color composite of two images.
#### padimage
Pad an array.
#### rescopeimage
Render the image in the region of the global coordinate system defined by a 
bounding box.

### Rendering Methods
#### imwarp
Apply an affine transform to an entire image.
#### meshwarp
Apply a piecewise affine transform to an image, as defined by a mesh.

### Offsets
To translate between intrinsic and global coordinate systems, a 2-element array
notes where the upper left pixel in the image is located in global space. When
a transform is applied to an image, the render method will return the
transformed image along with an offset, denoting the transformed image's
location in global space.

Transform methods will work on images with offsets. Offsets also allow an image
to be rendered in the global coordinate system to align with other images.

### Not included in this package
* Alpha masks
* Feature-based image registration (SURF or SIFT)

http://gbayer.com/development/moving-files-from-one-git-repository-to-another-preserving-history/
