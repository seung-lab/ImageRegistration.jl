# Julimaps
A set of tools for elastic image registration in Julia.

# Milestones
* 8/7 Stitch two tiles (8/9)
* 8/14 Stitch one section (8/16)
* 8/21 Align one wafer
  * Render stitched section
  * Load tiles
  * Load affines
  * Apply imwarp to tiles with affines
  * Downsample 2d array to arbitrary scale (not just factor of 2)
  * Blockmatch with subpixel accuracy
  * Filter spurious matches that are close to image edges
  * Determine other necessary filters
  * Create list of overlapping tile pairs based on affine transforms
  * Make seam inspection easier
* 8/28 Align one stack
  * Parallelize
* 9/4 Pre-alignment
  * Create overview object
* 9/11 Optimize
  * AWS
