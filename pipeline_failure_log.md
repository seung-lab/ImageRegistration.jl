| Date | Step | Problem | Root Cause | Fix |
| --- | --- | --- | --- | --- |
| 09/15/15 | Prealignment | 1,8 had slightly poor prealignment with 1,7 | 2/25 blockmatches failed to find the best match. In both cases, a dark patch in both images drove the correspondence, instead of the neurites. | Manually removed the bad matches | 
| 09/16/15 | Montage | 1,40 failed to montage well | 1,40 is slightly out of focus, so cross correlation values were lower than the threshold in montage blockmatching | Reran with a lower threshold |