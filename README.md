# Julimaps
A set of tools for elastic image registration in Julia.

# To do
* Find max of cross correlation with subpixel accuracy
* Find eignevalue ratio for principal curvature
* Try VoronoiDelaunay for mesh generation
* Create render section method to stitch all tiles together
* Optimize pa_warp2 (rectify number of operations with time spent in method)
* Test better interpolation method in piecewise affine warping
* Detect overlapping tile pairs based on affine transforms
* Include unit tests in MeshSolve
* Include unit tests in Mesh
* Include unit tests in piecewiseaffine_warp
* Create downsampling method
