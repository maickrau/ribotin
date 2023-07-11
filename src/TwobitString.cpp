#include <cassert>
#include "TwobitString.h"

TwobitString::TwobitString(const std::string& str)
{
	realSize = str.size();
	firstBit.resize((str.size()+63)/64, 0);
	secondBit.resize((str.size()+63)/64, 0);
	for (size_t i = 0; i < str.size(); i++)
	{
		size_t chunk = i/64;
		size_t off = i%64;
		switch(str[i])
		{
			case 'a':
			case 'A':
				break;
			case 'c':
			case 'C':
				firstBit[chunk] |= 1ull << off;
				break;
			case 'g':
			case 'G':
				secondBit[chunk] |= 1ull << off;
				break;
			case 't':
			case 'T':
				firstBit[chunk] |= 1ull << off;
				secondBit[chunk] |= 1ull << off;
				break;
			default:
				assert(false);
		}
	}
}

size_t TwobitString::size() const
{
	return realSize;
}

char TwobitString::getChar(size_t index) const
{
	size_t chunk = index / 64;
	size_t off = index % 64;
	return "ACGT"[2*((firstBit[chunk] >> off) & 1ull) + ((secondBit[chunk] >> off) & 1ull)];
}
