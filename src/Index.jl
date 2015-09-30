function is_adjacent(A::Index, B::Index, cross_wafer_guess = false)
	if sum(abs([A...] - [B...])) == 1
		return true;
	elseif cross_wafer_guess && (A[1]-B[1]==1 && A[2]==1) || (A[1]-B[1] == -1 && B[2]==1)
		return true;
	end
	return false;
end

function is_adjacent(Am::Mesh, Bm::Mesh)
  if abs(Am.index[3] - Bm.index[3]) + abs(Am.index[4] - Bm.index[4]) + abs(Am.index[2] - Bm.index[2]) == 1 return true; end
  return false;
end

function is_diagonal(Am::Mesh, Bm::Mesh)
  if abs(Am.index[3] - Bm.index[3]) + abs(Am.index[4] - Bm.index[4]) == 2 && Am.index[3] != Bm.index[3] && Am.index[4] != Bm.index[4] return true; end
  return false;
end

function get_overview_index(index::Index)
  return (index[1:2]..., OVERVIEW_INDEX, OVERVIEW_INDEX)
end

"""
INCOMPLETE

Lazy function to generate list of indices between indexA and indexB

Fixes to include:
* switch between wafers
* check indexA[3:4] against indexB[3:4]
"""
function create_index_range(indexA, indexB)
  return [(indexA[1], i, indexA[3:4]...) for i in indexA[2]:indexB[2]]
end

"""
Return zip object of an index and the index that follows it
"""
function create_sequential_index_pairs(indexA, indexB)
  indices = create_index_range(indexA, indexB)
  return zip(indices[1:end-1], indices[2:end])
end

"""
Test if an index is the first section in the stack
"""
function is_first_section(index)
  return index[1] == 1 && index[2] == 1
end