![Build Status](https://travis-ci.org/seung-lab/Julimaps.svg "travis")

# Julimaps
JULia IMAge Processing Suite:
A set of tools for elastic image registration in Julia.

# Process
(time in seconds per section)

| Step | Read | Match | Solve | Render | Write | Total | Review Method | Intervene Method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| premontage | 30 | 10 | 10 | 0 | 30 | 80 | (overlay tiles on overview) | NA |
| montage | 30 | 40 | 20 | 60 | 30 | 180 | (section overlay as checkerboard) | blockmatch image select |
| prealignment | 30 | 25 | 5 | 30 | 30 | 120 | overlay sections | blockmatch image select |
| alignment | 30 | 480 | 20 | 80 | 60 | 670 | movie of sections in FIJI | blockmatch image select |

# Milestones
* Montage two tiles (8/9)
* Montage one section (8/16)
* Prealign two sections (9/1)
* Align two sections (9/3)
* Align five sections (9/5)
* Align one wafer
* Align one stack
* Premontage one section
* Run on AWS