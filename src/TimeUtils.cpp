#include "declars.h"
#include "TimeUtils.h"


void printExecutionTime(std::chrono::high_resolution_clock::time_point start_time, std::chrono::high_resolution_clock::time_point end_time){

	auto execution_time_ms   = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	auto execution_time_sec  = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
	auto execution_time_min  = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();
	auto execution_time_hour = std::chrono::duration_cast<std::chrono::hours>(end_time - start_time).count();

	execution_time_ms   = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	execution_time_sec  = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
	execution_time_min  = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time).count();
	execution_time_hour = std::chrono::duration_cast<std::chrono::hours>(end_time - start_time).count();

	if (execution_time_hour > 0)
		cout << "" << execution_time_hour << "h, ";
	if (execution_time_min > 0)
		cout << "" << execution_time_min % 60 << "m, ";
	if (execution_time_sec > 0)
		cout << "" << execution_time_sec % 60 << "s, ";
	if (execution_time_ms > 0)
		cout << "" << execution_time_ms % long(1E+3) << " ms";
	cout << "\n";

}