#include "Stopwatch.h"



Stopwatch::Stopwatch() {
	size_meas = 0;
	max_time = 0;
	min_time = 1.7976931348623158e+308;
	is_start = false;
	average_time = 0;
}

void Stopwatch::start() {
	if (!is_start) {
		current_time = clock();
		is_start = true;
	}
}

void Stopwatch::end() {
	if (is_start) {
		is_start = false;
		++size_meas;
		current_time *= (-1);
		current_time += clock();
		if (current_time > max_time)
			max_time = current_time;
		if (current_time < min_time)
			min_time = current_time;
		average_time += current_time;
	}
}

Stopwatch::data Stopwatch::get_data() {
	data data_;
	data_.max_time = max_time / CLOCKS_PER_SEC;
	data_.min_time = min_time / CLOCKS_PER_SEC;
	data_.average_time = average_time / double(size_meas) / CLOCKS_PER_SEC;
	data_.last_current_time = current_time / CLOCKS_PER_SEC;

	return data_;
}


void Stopwatch::reset() {
	size_meas	 = 0;
	max_time	 = 0;
	min_time	 = 1.7976931348623158e+308;
	is_start	 = false;
	average_time = 0;
}