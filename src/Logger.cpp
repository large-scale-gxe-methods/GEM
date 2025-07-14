#include "Logger.h"


int CoutRedirector::overflow(int c) 
{
    if (c == '\n') 
    {
        spdlog::info("{}", buffer_);
        buffer_.clear();
    } else if (c != EOF) 
    {
        buffer_ += static_cast<char>(c);
    }
    return c;
}

int CoutRedirector::sync() 
{
    if (!buffer_.empty()) 
    {
        spdlog::info("{}", buffer_);
        buffer_.clear();
    }
    return 0;
}

int CerrRedirector::overflow(int c) 
{
    if (c == '\n') 
    {
        sync();
    } else if (c != EOF) 
    {
        buffer_ += static_cast<char>(c);
    }
    return c;
}

int CerrRedirector::sync() 
{
    if (!buffer_.empty()) 
    {
        spdlog::error("{}", buffer_);
        buffer_.clear();
    }
    return 0;
}

void LoggerSetup::init(const std::string& filename) 
{
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, true);
    std::vector<spdlog::sink_ptr> sinks{console_sink, file_sink};

    auto logger = std::make_shared<spdlog::logger>("multi_logger", sinks.begin(), sinks.end());
    logger->set_level(spdlog::level::debug);
    spdlog::set_default_logger(logger);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");

    static CoutRedirector cout_redirector;
    std::cout.rdbuf(&cout_redirector);
    static CerrRedirector cerr_redirector;
    std::cerr.rdbuf(&cerr_redirector);
}
