#ifndef ALL_PRM_H
#define ALL_PRM_H
#include "Func.h"
#include <string>

class all_prm {
	

public:
	all_prm(double t) : func(t , std::string("dfhfh")) {}

	double get_e() { return e; }
	void set_e(double new_e) { e = new_e; }

	double e;
	double t;
	double s;
	double error_ext1;
	double error_ext2;
	double error_in1;
	double error_in2;

	UpBounds   * c_ub;
	DownBounds * c_db;

	double(*f)(double x, void* param);

	std::function<double(double, double, double)> up;
	std::function<double(double, double, double)> down;

	double v_ub;
	double v_db;

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

	// tmp_value:
	int  type_bounds_up;
	int  type_bounds_down;
};

#endif