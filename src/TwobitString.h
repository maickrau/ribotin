#ifndef TwobitString_h
#define TwobitString_h

#include <string>
#include <vector>

class TwobitString
{
public:
	TwobitString() = default;
	TwobitString(const std::string& str);
	__attribute__((always_inline)) inline std::pair<uint64_t, uint64_t> getBitvectorsPossiblyTruncated(size_t startIndex) const
	{
		size_t chunk = startIndex / 64;
		size_t off = startIndex % 64;
		uint64_t first = firstBit[chunk] >> off;
		uint64_t second = secondBit[chunk] >> off;
		return std::make_pair(first, second);
	}
	size_t size() const;
	char getChar(size_t index) const;
private:
	std::vector<uint64_t> firstBit;
	std::vector<uint64_t> secondBit;
	size_t realSize;
};

#endif
