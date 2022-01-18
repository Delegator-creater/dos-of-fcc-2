#ifndef _STOPWATCH_H
#define _STOPWATCH_H
#include <ctime>
#include <limits.h>
#include <iostream>

class Stopwatch {

	double max_time;
	double min_time;
	double average_time;
	double current_time;
	size_t size_meas;
	bool   is_start;
public:

	Stopwatch();

	void start();

	void end();

	using data = struct {
		double max_time;
		double min_time;
		double average_time;
		double last_current_time;
	};

	data get_data();

	template<typename OUT>
	void print( OUT & io_out , const char name[] )
	{
		auto data = get_data();
		io_out << "->" << name << data.max_time << " " << data.min_time << " " << data.average_time << "\n";
	}

	void reset();

};

#endif