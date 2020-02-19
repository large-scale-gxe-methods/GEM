#pragma once

#ifndef TIMEUTILS_H
#define TIMEUTILS_H


#include <chrono>

class Time {

public:
};

void printExecutionTime(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point end_time);
#endif

