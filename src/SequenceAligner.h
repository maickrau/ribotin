#ifndef SequenceAligner_h
#define SequenceAligner_h

#include <vector>
#include <string>
#include <unordered_map>
#include "ClusterMisc.h"

size_t getEditDistanceWfa(const std::string& left, const std::string& right);
size_t getEditDistanceWfa(const PathSequenceView& left, const PathSequenceView& right, const size_t maxEdits);
size_t getEditDistance(const std::vector<Node>& left, const size_t leftIndex, const std::vector<Node>& right, const size_t rightIndex, const GfaGraph& graph, const std::unordered_map<Node, size_t>& pathStartClip, const std::unordered_map<Node, size_t>& pathEndClip, const size_t maxEdits, const std::unordered_set<size_t>& coreNodes, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodeCountIndex, const std::vector<phmap::flat_hash_map<Node, size_t>>& nodePosIndex, std::unordered_map<std::pair<std::vector<Node>, std::vector<Node>>, size_t>& memoizedEditDistances);

#endif
