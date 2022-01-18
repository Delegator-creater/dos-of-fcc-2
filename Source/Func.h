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
#include "out_err.h"
#include "Stopwatch.h"

static Stopwatch stch_Func_Q;
class Func {
	std::ofstream out;
	using _str = std::string;



public:
	Func(double t , std::string prm_start) : out(_str("rho/logs/")+ prm_start + _str(" logs(t = ") + std::to_string(t) + _str(").txt"), std::ios::out) {}

	double operator()(ARG_3, double)  ;

	template<const int case_>
	double Q(ARG_3, double d, double dir, double zeta_ast) {
		stch_Func_Q.start();
		try {
			double z = zeta_ast + dir * d;
			static double a1 = 1 / (2 * std::pow(3.14159, 3));
			double a2 = case_ == 5 ? 1 / std::sqrt(d * (1 - d)) : 1 / std::sqrt((1 - z) * z);
			double a3 = case_ == 5 ? 4 * std::abs(s) - sqr_(1 - std::abs(s)) * d
				: sqr_(1 + std::abs(s)) - sqr_(1 - std::abs(s)) * z;
			double R_class = (1. + t - 2. * sqr_(t)) * sqr_(s) -
				t * (e / 4. + (1. + 2. * t) * sqr_(1. - std::abs(s)) * z);
			double R = case_ != 3 ? R_class : R_f<0>(e, t, s, d * dir);
			double phi_or_phi;
			switch (case_)
			{
			case 1: {
				phi_or_phi = psi(arg_psi<-1>(e, t, s, dir * d));
				break;
			}
			case 2: {
				phi_or_phi = psi(arg_psi<1>(e, t, s, dir * d));
				break;
			}
			default: {
				phi_or_phi = phi((std::sqrt(R) + s) / t);
				break;
			}
			}
			double result = a1 * a2 / std::sqrt(a3 * R) * phi_or_phi;
			
			if ( ((a3 * R) <= 0) || 
					( 
						(case_ > 2) &&
							(
								(R <= 0) ||
								(std::abs((std::sqrt(R) + s) / t) >= 1)
							)
					)
				) {
				stch_Func_Q.end();
				//err << "in funcion Q<case_>, case_ = " << case_ << "; prm: e=" << e << " t =" << t << " s =" << s << " d=" << d << " dir =" << dir << " a1="
				//	<< a1 << " a2 =" << a2 << " a3=" << a3 << " R=" << R << " phi_or_phi=" << phi_or_phi << " zeta_ast=" << zeta_ast << " \n";
				result = 0;
			}
			if ((d * (1 - d) < 0) || (z * (1 - z) < 0)) {
				err_ << "z=" << z << " d=" << d << "\n";
				stch_Func_Q.end();
				throw 1;
			}
			stch_Func_Q.end();
			return result;
		}
		catch (int){
			stch_Func_Q.end();
			throw 1;
			//return 0;
		}

	}



	~Func() {
		out.close();
	}
};

using Regime = struct {
	char regime, bound;
};

using Integrand2_prms = struct {
	double e, tau, s;
	short sgn; // dir
	Regime R; // режим , граница
};

double f_zeta(double s, double e, double tau, char t);

double integrand2(double x, void* prms);
#endif
