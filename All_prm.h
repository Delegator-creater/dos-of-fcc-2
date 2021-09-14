#ifndef ALL_PRM_H
#define ALL_PRM_H
#include "Func.h"
#include <string>

class all_prm {
	double e;

public:
	all_prm(double t) : func(t , std::string("")) {}

	double get_e() { return e; }
	void set_e(double new_e) { e = new_e; }

	double t;
	double s;
	double error_ext1;
	double error_ext2;
	double error_in1;
	double error_in2;



	double(*f)(double x, void* param);

	std::function<double(double, double, double)> up;
	std::function<double(double, double, double)> down;

	Func func;

	gsl_integration_workspace* memory1;
	gsl_integration_workspace* memory2;

	size_t size_memory_1;
	size_t size_memory_2;
	using vfi = struct
	{
		double s;
		double int_for_zeta;
		double zeta1;
		double zeta2;
	};
	std::deque<vfi> value_first_int;

};

#endif