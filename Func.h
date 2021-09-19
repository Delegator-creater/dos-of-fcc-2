#ifndef FUNC_H
#define FUNC_H
#include "Macros.h"
#include <fstream>
#include <string>
#include <set>
#include <fstream>
#include <iostream>
#include "spt_func.h"
#include "phi.h"

class Func {
	std::ofstream out;
	using _str = std::string;



public:
	Func(double t , std::string prm_start) : out(_str("rho/logs/")+ prm_start + _str(" logs(t = ") + std::to_string(t) + _str(").txt"), std::ios::out) {}

	double operator()(ARG_3, double)  ;

	template<const int case_>
	double Q(ARG_3, double d, double dir, double zeta_ast) {
		try {
			double z = zeta_ast + dir * d;
			static double a1 = 1 / (2 * std::pow(3.14159, 3));
			double a2 = case_ == 5 ? 1 / std::sqrt(d * (1 - d)) : 1 / std::sqrt((1 - z) * z);
			double a3 = case_ == 5 ? std::sqrt(4 * std::abs(s) - sqr_(1 - std::abs(s)) * d)
				: std::sqrt(sqr_(1 + std::abs(s)) - sqr_(1 - std::abs(s)) * z);
			double R_class = (1. + t - 2. * sqr_(t)) * sqr_(s) -
				t * (e / 4. + (1. + 2. * t) * sqr_(1. - std::abs(s)) * sqr_(z));
			double R = case_ != 3 ? R_class : R_f<0>(e, t, s, d * dir);
			double phi_or_phi;
			switch (case_)
			{
			case 1: {
				phi_or_phi = psi(arg_psi<-1>(e, t, s, dir * d));
				break;
			}
			case 2: {
				phi_or_phi = psi(arg_psi<+1>(e, t, s, dir * d));
				break;
			}
			default: {
				phi_or_phi = phi((std::sqrt(R) + s) / t);
				break;
			}
			}
			double result = a1 * a2 / std::sqrt(a3 * R) * phi_or_phi;
			//static bool flag = true;
			if ( ((a3 * R) <= 0) || ( (case_ > 2) && (R <= 0))) {
				//std::cerr << "in funcion Q<case_>, case_ = " << case_ << "; prm: e=" << e << " t =" << t << " s =" << s << " d=" << d << " dir =" << dir << " a1="
				//	<< a1 << " a2 =" << a2 << " a3=" << a3 << " R=" << R << " phi_or_phi=" << phi_or_phi << " zeta_ast=" << zeta_ast << " \n";
				result = 0;
			}
			if ((d * (1 - d) < 0) || (z * (1 - z) < 0)) {
				std::cerr << "z=" << z << " d=" << d << "\n";
				throw 1;
			}
			return result;
		}
		catch (int){
			throw 1;
			//return 0;
		}

	}

	~Func() {
		out.close();
	}
};
#endif
