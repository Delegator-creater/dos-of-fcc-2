#ifndef MY_INTEGRATION
#define MY_INTEGRATION

#include <functional>
#include "All_prm.h"
#include "out_err.h"
#include <string>


#ifdef WIN32

void gsl_set_error_handler_off() {
}

#include <string>
inline std::string gsl_strerror(int cod) {
	return std::string("");
}

static const size_t start_points = 100;

template<const int>//
static double integration( double(*f)(double, void*),  void*  prm, double a, double  b, const double error1, const double error2, int& cod_error) {
	a += error1;
	b -= error1;
	size_t n = start_points;
	all_prm* allp = (all_prm *)prm;
	auto const sumiring = [&]() {
		double h = (b - a) / n;
		double  sum = 0;
		for (double count = a; count <= b; count += h)
			sum += (f(count, prm) + f(count + h, prm)) / 2;
		return sum * h;
	};

	double sum1 = sumiring();
	n *= 2;
	double sum2 = sumiring();
	while ( (abs(sum1 - sum2) >= error2 ) && (n < allp->size_memory_1) ) {
		n *= 2;
		sum1 = sum2;
		sum2 = sumiring();
	}

	return sum2; //double integration(const std::function<double(double , void *)>, const void*, const double, const double, const double, const double);
}
#endif

///////////////////////////////////////////////////


#if defined(__linux)
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

gsl_function make_gsl_function(double (* function) (double,void *) , void* prm){
	gsl_function gsl_fun;
	gsl_fun.function = function;
	gsl_fun.params   = prm;
	return gsl_fun;
}

template<const int prm_integr>
double integration( double(*f)(double , void *), void* prm, double a, double  b, const double error1, const double error2, int& cod_error) {
	err_ << std::string("Unimplemented template parameter specified for function integration");
}

template<>
inline double integration<1>( double(*f)(double, void*), void* prm, double a, double  b, const double error1, const double error2, int& cod_error) {
	
	double bounds[2]({a , b});
	double abserr;
	double result;
	auto f_gsl = make_gsl_function(f , prm);
	all_prm* allp = (all_prm *)prm;
	static const auto memory = gsl_integration_workspace_alloc(allp->size_memory_2);
	cod_error = gsl_integration_qags(&f_gsl,
		bounds[0],
		bounds[1],
		error1,
		error2,
		allp->size_memory_2,
		memory,
		&result,
		&abserr);
	return result;
}

template<>
inline double integration<2>(double(*f)(double, void*), void* prm, double a, double  b, const double error1, const double error2, int& cod_error) {
	
	double abserr;
	double result;
	all_prm* allp = (all_prm *)prm;
	auto f_gsl = make_gsl_function(f , prm);
	static const auto memory = gsl_integration_workspace_alloc(allp->size_memory_2);
	cod_error = gsl_integration_qags(&f_gsl,
		a,
		b,
		error1,
		error2,
		allp->size_memory_2,
		memory,
		&result,
		&abserr);
	return result;
}

#endif

#endif