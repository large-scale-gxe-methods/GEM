#pragma once

#include <iostream>
#include <streambuf>
#include <memory>
#include <string>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <vector>

class CoutRedirector : public std::streambuf 
{
    protected:
        int overflow(int c) override;
        int sync() override;
    private:
        std::string buffer_;
};

class CerrRedirector : public std::streambuf 
{
    protected:
        int overflow(int c) override;
        int sync() override;
    private:
        std::string buffer_;
};

class LoggerSetup 
{
    public:
        static void init(const std::string& filename = "log.txt");
};
