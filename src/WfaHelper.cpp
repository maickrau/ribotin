#include <limits>
#include <cassert>
#include "WfaHelper.h"

namespace Wfa
{
	const std::pair<size_t, size_t> uninitialized { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
	std::pair<size_t, size_t> getRealMatrixPosition(std::pair<size_t, size_t> wfaMatrixPosition, size_t offset)
	{
		size_t i = offset;
		size_t j = wfaMatrixPosition.second + offset - wfaMatrixPosition.first;
		return std::make_pair(i, j);
	}
	std::pair<size_t, size_t> getPredecessorDiagonalOffset(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& matrix)
	{
		if (edits > wfaMatrixPosition.first) return uninitialized;
		size_t newI = wfaMatrixPosition.first - edits;
		assert(newI < matrix.size());
		size_t newJ = wfaMatrixPosition.second - edits;
		if (verticalOffset > 0) newJ += verticalOffset;
		if (horizontalOffset > 0) newJ -= horizontalOffset;
		if (newJ >= matrix[newI].size()) return uninitialized;
		return matrix[newI][newJ];
	}
	bool canBacktrace(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& predecessorMatrix, const WfaMatrix& currentMatrix, size_t refLength, size_t queryLength)
	{
		assert(currentMatrix[wfaMatrixPosition.first][wfaMatrixPosition.second] != uninitialized);
		assert(wfaMatrixPosition.first < currentMatrix.size());
		assert(wfaMatrixPosition.second < currentMatrix[wfaMatrixPosition.first].size());
		if (verticalOffset >= wfaMatrixPosition.first) return false;
		auto value = getPredecessorDiagonalOffset(wfaMatrixPosition, edits, verticalOffset, horizontalOffset, predecessorMatrix);
		if (value == uninitialized) return false;
		if (value.second + verticalOffset == currentMatrix[wfaMatrixPosition.first][wfaMatrixPosition.second].first) return true;
		return false;
	}
	void updateMatrix(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& sourceMatrix, WfaMatrix& targetMatrix, size_t refLength, size_t queryLength)
	{
		auto value = getPredecessorDiagonalOffset(wfaMatrixPosition, edits, verticalOffset, horizontalOffset, sourceMatrix);
		if (value == uninitialized) return;
		value.second += verticalOffset;
		auto realPos = getRealMatrixPosition(wfaMatrixPosition, value.second);
		if (realPos.first >= refLength || realPos.second >= queryLength) return;
		if (targetMatrix[wfaMatrixPosition.first][wfaMatrixPosition.second] == uninitialized)
		{
			targetMatrix[wfaMatrixPosition.first][wfaMatrixPosition.second].first = value.second;
		}
		else if (targetMatrix[wfaMatrixPosition.first][wfaMatrixPosition.second].first < value.second)
		{
			targetMatrix[wfaMatrixPosition.first][wfaMatrixPosition.second].first = value.second;
		}
	}
	std::pair<size_t, size_t> getBacktracePosition(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& sourceMatrix, const WfaMatrix& targetMatrix)
	{
		assert(edits <= wfaMatrixPosition.first);
		size_t newI = wfaMatrixPosition.first - edits;
		assert(newI < targetMatrix.size());
		size_t newJ = wfaMatrixPosition.second - edits;
		if (verticalOffset > 0) newJ += verticalOffset;
		if (horizontalOffset > 0) newJ -= horizontalOffset;
		assert(newJ < targetMatrix[newI].size());
		return std::make_pair(newI, newJ);
	}
};
