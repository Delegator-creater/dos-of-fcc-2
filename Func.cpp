#ifndef FUNC_CPP
#define FUNC_CPP
#include "Func.h"
#include "spt_func.h"
#include "phi.h"
#include <fstream>
#include <iostream>


double Func::operator()(ARG_3, double x)  {

	double R = (1. + t - 2. * sqr_(t))*sqr_(s) - 
	t * (e / 4. + (1. + 2.* t)* sqr_(1. - std::abs(s) )* sqr_(x));

	double M = sqr_(1.+ std::abs(s)) - sqr_(1.- std::abs(s))* sqr_(x) ;

	double arg = (sqrt(R) + s) / t;
	

	//if ( set_eps.size()==3 )
	//	if ( abs(s + 0.6) < 0.001 )
	//	std::cerr << "\t" << x << "\t" << s << "\t"
	//	<< t << "\t" << e << "\t" << arg << "\t" << M << "\t" << R << "\n";

	if ( (std::abs(arg) >= 1) || (R * M <= 0 )|| (R <= 0)) {


		//if (out.is_open())
		//{//"(xi,s,t,e,f,arg,M,R) 


		//}
		
		return 0;

	}
	double f_ = 1.;
	f_ *= phi(arg);
	f_ /= sqrt(M*R*(1 - sqr_(x)) );
	f_ = f_ > 10E5 ? 10E5 : f_;
	return f_;
}
#endif
