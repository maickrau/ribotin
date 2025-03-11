#ifndef RibotinUtils_h
#define RibotinUtils_h

#include <cstdint>
#include <vector>
#include <string>
#include <chrono>

bool fileExists(std::string filename);
std::vector<std::string> split(const std::string& str, char separator);
size_t getSequenceLength(const std::string& filename);
std::chrono::time_point<std::chrono::steady_clock> getTime();
std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end);
inline int popcount(uint64_t x)
{
	//https://gcc.gnu.org/onlinedocs/gcc-4.8.4/gcc/X86-Built-in-Functions.html
	// return __builtin_popcountll(x);
	//for some reason __builtin_popcount takes 21 instructions so call assembly directly
	__asm__("popcnt %0, %0" : "+r" (x));
	return x;
}

#endif
