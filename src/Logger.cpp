#include "Logger.h"

Logger Logger::Log { std::cerr };

Logger::Logger(std::ostream& stream) :
	stream(stream),
	verbosity(0)
{
}

void Logger::setVerbosity(const size_t newVerbosity)
{
	verbosity = newVerbosity;
}

Logger::EnabledOrDisabledLogger::EnabledOrDisabledLogger(bool disabled, std::ostream& stream) :
	disabled(disabled),
	stream(stream)
{
}

Logger::EnabledOrDisabledLogger Logger::log(const LogLevel level)
{
	if (level <= verbosity) return EnabledOrDisabledLogger { false, stream };
	return EnabledOrDisabledLogger { true, stream };
}

