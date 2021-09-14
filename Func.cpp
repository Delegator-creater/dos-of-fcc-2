#ifndef FUNC_CPP
#define FUNC_CPP
#include "Func.h"
#include "spt_func.h"
#include "phi.h"
#include <fstream>
#include <iostream>


double Func::operator()(ARG_3, double x)  {

	double R_ = (1. + t - 2. * sqr_(t))*sqr_(s) - 
	t * (e / 4. + (1. + 2.* t)* sqr_(1. - std::abs(s) )* sqr_(x));

	double M = sqr_(1.+ std::abs(s)) - sqr_(1.- std::abs(s))* sqr_(x) ;

	double arg = (sqrt(R_) + s) / t;
	

	//if ( set_eps.size()==3 )
	//	if ( abs(s + 0.6) < 0.001 )
	//	std::cerr << "\t" << x << "\t" << s << "\t"
	//	<< t << "\t" << e << "\t" << arg << "\t" << M << "\t" << R << "\n";

	if ( (std::abs(arg) >= 1) || (R_ * M <= 0 )|| (R_ <= 0)) {


		//if (out.is_open())
		//{//"(xi,s,t,e,f,arg,M,R) 


		//}
		
		return 0;

	}
	double f_ = 1.;
	f_ *= phi(arg);
	f_ /= sqrt(M*R_*(1 - sqr_(x)) );
	return f_;
}
#endif
