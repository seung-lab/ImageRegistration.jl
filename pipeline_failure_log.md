| Date | Step | Problem | Root Cause | Fix | Rerendered |
| --- | --- | --- | --- | --- | --- |
| 09/15/15 | Prealignment | 1,8 had slightly poor prealignment with 1,7 | 2/25 blockmatches failed to find the best match. In both cases, a dark patch in both images drove the correspondence, instead of the neurites. | Manually removed the bad matches | Yes |
| 09/16/15 | Montage | 1,40 failed to montage well | 1,40 is slightly out of focus, so cross correlation values were lower than the threshold in montage blockmatching | Reran with a lower threshold | Yes |
| 09/16/15 | Alignment | 1,2 failed locally in the interior aligning to 1,1 | There is a narrow scratch that runs through tiles 2 & 3  | Reran with a finer mesh | Yes |
| 09/17/15 | Prealignment | 1,33 failed to prealign well with 1,32 | 1,33 is translated by nearly 3000 px from 1,32 | Reran with a larger search radius | Yes |
| 09/19/15 | Montage | 1,103 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,120 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,130 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,151 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,153 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,154 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,155 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,165 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,166 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,167 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/19/15 | Montage | 1,52 montaged with a blank tile | Blank tile from the microscope was not removed in premontage & noisy image | Removed blank tile from premontage offset file | Yes |
| 09/19/15 | Montage | 1,54 montaged with a blank tile | Blank tile from the microscope was not removed in premontage & noisy image | Removed blank tile from premontage offset file | Yes |
| 09/20/15 | Prealignment | 1,54 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/20/15 | Prealignment | 1,136 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/20/15 | Prealignment | 1,139 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/20/15 | Prealignment | 1,140 failed to prealign due to lack of points | Large offset | Manually entered offset | Yes |
| 09/22/15 | Montage | 1,56 failed to montage well | Noisy image | Reran with a lower xc threshold | Yes |
| 09/28/15 | Prealignment | 1,119-121 has bad prealignment due to lack of points | Wide mesh | Reran with fine mesh | No |
| 09/28/15 | Prealignment | 1,129-132 has bad prealignment due to lack of points | Wide mesh | Reran with fine mesh | No |
| 09/28/15 | Prealignment | 1,141-143 has bad prealignment due to lack of points | Wide mesh + large offset | Reran with fine mesh / manually entered offset | No |
| 09/28/15 | Prealignment | 1,155-158 has bad prealignment due to lack of points | Large offset | Manually entered offset | No |
| 09/28/15 | Prealignment | 1,167-168 has bad prealignment due to lack of points | Large offset | Manually entered offset | No |
| 09/28/15 | Premontage | 2,3 has bad premontage due to an extra tile | NA | NA | No |
| 09/28/15 | Premontage | 2,52 has bad premontage due to an extra tile | NA | NA | No |
| 09/29/15 | Premontage | 3,22 has bad premontage due to an extra tile | NA | NA | No |
| 09/29/15 | Premontage | 3,90 has bad premontage | NA | NA | No |
| 09/29/15 | Premontage | 3,133 has bad premontage | NA | NA | No |
| 09/30/15 | Premontage | 3,147 has bad premontage | NA | NA | No |

