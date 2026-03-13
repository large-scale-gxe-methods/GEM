#pragma once

#include <iostream>
#include <streambuf>
#include <memory>
#include <string>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <vector>

/**
 * @class CoutRedirector
 * @brief Redirects standard output (std::cout) into the spdlog logging system.
 *
 * This class inherits from std::streambuf and overrides the virtual functions
 * overflow() and sync() in order to intercept characters written to std::cout.
 *
 * Typical use:
 * - Installed via std::cout.rdbuf(&redirector)
 * - Enables capturing of all std::cout output into log files and console sinks.
 */
class CoutRedirector : public std::streambuf 
{
    protected:
        int overflow(int c) override;
        int sync() override;
    private:
        std::string buffer_;
};

/**
 * @class CerrRedirector
 * @brief Redirects standard error output (std::cerr) into the spdlog logging system.
 *
 * This class captures all output written to std::cerr and forwards it to spdlog
 * using the error-level logger.
 *
 * Characters are buffered until a newline or flush occurs.
 */
class CerrRedirector : public std::streambuf 
{
    protected:
        int overflow(int c) override;
        int sync() override;
    private:
        std::string buffer_;
};

/**
 * @class LoggerSetup
 * @brief Initializes a global spdlog logger and redirects std::cout/std::cerr.
 *
 * This utility class configures spdlog with:
 * - Console logging (colored output)
 * - File logging (persistent log file)
 * - Custom log formatting
 * - Automatic flushing behavior
 *
 * After initialization, all output written to std::cout and std::cerr is captured
 * and forwarded into spdlog, allowing unified logging of both C++ and library output.
 */
class LoggerSetup 
{
    public:
        static void init(const std::string& filename = "log.txt");
};
