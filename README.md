![Build Status](https://travis-ci.org/seung-lab/Julimaps.svg "travis")

# Julimaps
JULia IMAge Processing Suite:
A set of tools for elastic image registration in Julia.

# Process
(time in seconds per section)

| Step | Read | Match | Solve | Render | Write | Total | Review Method | Intervene Method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| premontage | 30 | 10 | 10 | 0 | 30 | 80 | (overlay tiles on overview) | NA |
| montage | 30 | 100 | 20 | 60 | 30 | 240 | (section overlay as checkerboard) | blockmatch image select |
| prealignment | 30 | 25 | 5 | 30 | 30 | 120 | overlay sections | blockmatch image select |
| alignment | 30 | 100 | 20 | 60 | 30 | 240 | movie of sections in FIJI | blockmatch image select |

# Milestones
* 8/7 Stitch two tiles (8/9)
* 8/14 Stitch one section (8/16)
* 8/21 Align one wafer
  * Render stitched section (Tommy) X
  * Load images (Dodam) X
  * Apply imwarp to tiles with affines (Tommy + Dodam) X
  * Downsample 2d array to arbitrary scale (not just factor of 2) (Shang) X
  * Implement rough filter for spurious matches that are close to image edges (Dodam) X
  * Determine other useful match filters (Sebastian)
  * Create list of overlapping tile pairs based on affine transforms (Dodam) X
  * Make seam inspection easier (Tommy)
  * Update mesh code to handle stitched sections (Dodam) X
  * Store cross correlation plots for the bad correspondences X
  * Pre-align two stitched section images (Shang) X
* 8/28 Elastically align one stack (piriform)
  * Parallelize
* 9/11 Pre-align & elastically align one stack (zebrafish)
  * AWS
* Other tasks:
  * Blockmatch with subpixel accuracy
  * Make rendering a section parallel
