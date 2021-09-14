#ifndef FUNC_H
#define FUNC_H
#include "Macros.h"
#include <fstream>
#include <string>
#include <set>
#include <fstream>
#include <iostream>
#include "spt_func.h"

class Func {
	std::ofstream out;
	using _str = std::string;


	template<const int case_>
	double Q(ARG_3, double d, double dir , double zeta_ast) {
		double z  = zeta_ast + dir*d;
		static double a1 = 1 / (2 *std::pow(3.14159, 3) );
		double a2 = case_ == 5 ? 1 / std::sqrt(d*(1-d)) : 1 / std::sqrt((1 - z) * z);
		double a3 = case_ == 5 ? std::sqrt(4 * std::abs(s) - sqr_(1 - std::abs(s))*d 
							   : std::sqrt( sqr_(1+std::abs(s))-sqr_(1-std::abs(s))*z );
		double R_class = (1. + t - 2. * sqr_(t)) * sqr_(s) -
			t * (e / 4. + (1. + 2. * t) * sqr_(1. - std::abs(s)) * sqr_(z));
		double R = case_ != 3 ? R_class : -t*(1+2*t)*sqr_(1-std::abs(s)) * dir *d;
		double phi_or_phi;
		switch (case_)
		{
			case 1: {
				phi_or_phi = psi( arg_psi<-1>(ARG_3, dir * d) );
				break;
			}
			case 2: {
				phi_or_phi = psi( arg_psi<+1>(ARG_3, dir * d) );
				break;
			}default: {
				phi_or_phi = phi( (std::sqrt(R) + s) / t );
				break;
			}
		}
		return a1 * a2 / std::sqrt(a3 * R) * phi_or_phi;

	}

public:
	Func(double t , std::string prm_start) : out(_str("rho/logs/")+ prm_start + _str(" logs(t = ") + std::to_string(t) + _str(").txt"), std::ios::out) {}

	double operator()(ARG_3, double)  ;

	~Func() {
		out.close();
	}
};
#endif
