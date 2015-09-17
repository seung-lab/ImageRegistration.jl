| Date | Step | Problem | Root Cause | Fix |
| --- | --- | --- | --- | --- |
| 09/15/15 | Prealignment | 1,8 had slightly poor prealignment with 1,7 | 2/25 blockmatches failed to find the best match. In both cases, a dark patch in both images drove the correspondence, instead of the neurites. | Manually removed the bad matches | 
| 09/16/15 | Montage | 1,40 failed to montage well | 1,40 is slightly out of focus, so cross correlation values were lower than the threshold in montage blockmatching | Reran with a lower threshold |
| 09/16/15 | Alignment | 1,2 failed locally in the interior aligning to 1,1 | There is a narrow scratch that runs through tiles 2 & 3  | Reran with a finer mesh |
| 09/17/15 | Prealignment | 1,33 failed to prealign well with 1,32 | 1,33 is translated by nearly 3000 px from 1,32 | Reran with a larger search radius |