#ifndef WfaHelper_h
#define WfaHelper_h

#include <vector>
#include <tuple>

namespace Wfa
{
	using WfaMatrix = std::vector<std::vector<std::pair<size_t, size_t>>>;
	std::pair<size_t, size_t> getRealMatrixPosition(std::pair<size_t, size_t> wfaMatrixPosition, size_t offset);
	std::pair<size_t, size_t> getPredecessorDiagonalOffset(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& matrix);
	bool canBacktrace(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& predecessorMatrix, const WfaMatrix& currentMatrix, size_t refLength, size_t queryLength);
	void updateMatrix(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& sourceMatrix, WfaMatrix& targetMatrix, size_t refLength, size_t queryLength);
	std::pair<size_t, size_t> getBacktracePosition(std::pair<size_t, size_t> wfaMatrixPosition, size_t edits, size_t verticalOffset, size_t horizontalOffset, const WfaMatrix& sourceMatrix, const WfaMatrix& targetMatrix);
};

#endif
