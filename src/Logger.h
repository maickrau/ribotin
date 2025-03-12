#ifndef Logger_h
#define Logger_h

#include <cstddef>
#include <iostream>
#include <ostream>

class Logger
{
public:
	class EnabledOrDisabledLogger
	{
	public:
		EnabledOrDisabledLogger(bool enable, std::ostream& stream);
		EnabledOrDisabledLogger() = delete;
		EnabledOrDisabledLogger& operator=(const EnabledOrDisabledLogger&) = delete;
		EnabledOrDisabledLogger(const EnabledOrDisabledLogger&) = delete;
		template <typename T>
		EnabledOrDisabledLogger& operator<<(T info)
		{
			if (disabled) return *this;
			stream << info;
			return *this;
		}
		// for std::endl
		EnabledOrDisabledLogger& operator<<(decltype(std::cerr)& (*info)(decltype(std::cerr)&))
		{
			if (disabled) return *this;
			stream << info;
			return *this;
		}
	private:
		bool disabled;
		std::ostream& stream;
	};
	enum LogLevel : size_t
	{
		Always = 0,
		DebugInfo = 1,
		DetailedDebugInfo = 2
	};
	Logger(std::ostream& stream);
	Logger& operator=(const Logger&) = delete;
	Logger(const Logger&) = delete;
	EnabledOrDisabledLogger log(const LogLevel level);
	void setVerbosity(size_t newVerbosity);
	static Logger Log;
private:
	std::ostream& stream;
	size_t verbosity;
};

#endif
