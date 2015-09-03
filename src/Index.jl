function isAdjacent(A::Index, B::Index, cross_wafer_guess = false)
	if sum(abs([A...] - [B...])) == 1
		return true;
	elseif cross_wafer_guess && (A[1]-B[1]==1 && A[2]==1) || (A[1]-B[1] == -1 && B[2]==1)
		return true;
	end
	return false;
end
