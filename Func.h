#ifndef FUNC_H
#define FUNC_H
#include "Macros.h"
#include <fstream>
#include <string>
#include <set>
#include <fstream>
#include <iostream>

class Func {
	std::ofstream out;
	using _str = std::string;
public:
	Func(double t , std::string & prm_start) : out(_str("rho/logs/")+ prm_start + _str(" logs(t = ") + std::to_string(t) + _str(").txt"), std::ios::out) {}

	double operator()(ARG_3, double)  ;

	~Func() {
		out.close();
	}
};
#endif
