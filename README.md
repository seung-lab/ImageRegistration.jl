# Julimaps
A set of tools for elastic image registration in Julia.

# Milestones
* 8/7 Stitch two tiles (8/9)
* 8/14 Stitch one section (8/16)
* 8/21 Align one wafer
  * Render stitched section (Tommy)
  * Load tiles (Tommy + Dodam)
  * Load affines (Tommy + Dodam)
  * Apply imwarp to tiles with affines (Tommy + Dodam)
  * Downsample 2d array to arbitrary scale (not just factor of 2) (Shang)
  * Blockmatch with subpixel accuracy (Shang)
  * Filter spurious matches that are close to image edges (Dodam)
  * Determine other necessary match filters (Dodam)
  * Create list of overlapping tile pairs based on affine transforms (Tommy)
  * Make seam inspection easier (Tommy)
  * Update mesh code to handle stitched sections (Dodam)
* 8/28 Align one stack
  * Parallelize
* 9/4 Pre-alignment
  * Create overview object
* 9/11 Optimize
  * AWS
