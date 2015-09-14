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

