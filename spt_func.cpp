#include "spt_func.h"
#include <deque>
#include <string>

std::deque<std::string> split(std::string str, const char * delimiter) {
	std::deque<std::string> result;
	size_t index = 0;
	for (size_t i = 0; i < str.size(); ++i)
		if (str[i] == *delimiter) {
			result.push_back(str.substr(index, i - index));
			index = i + 1;
		}
	return result;
}

double sqr_(double x) {
	return x*x;
}

