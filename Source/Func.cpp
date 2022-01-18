#ifndef FUNC_CPP
#define FUNC_CPP
#include "Func.h"
#include "spt_func.h"
#include "phi.h"
#include <fstream>
#include <iostream>
#include "Bounds.h"


double Func::operator()(ARG_3, double x)  {

	double R_ = (1. + t - 2. * sqr_(t))*sqr_(s) - 
	t * (e / 4. + (1. + 2.* t)* sqr_(1. - std::abs(s) )*x);

	double M = sqr_(1.+ std::abs(s)) - sqr_(1.- std::abs(s))* x;

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
	double f_ = 1./(2*std::pow(3.14159,3));
	f_ *= phi(arg);
	f_ /= sqrt(M*R_*(1 -x)*x );
	return f_;
}

double f_zeta(double s, double e, double tau, char t) {
	double E_arg;
	switch (t) {
		case 'o': {
			E_arg = 4 * sqr_(s) / tau;
			break;
		}
		case '+': {
			E_arg = 8. *  s - 4. * tau;
			break;
		}
		case '-': {
			E_arg = 8. * (-1.) * s - 4. * tau;
			break; 
		}

		default: {
			fprintf(stderr, "Calculation::zeta : Wrong argument, t = %c\n", t);
			exit(1);
		}
	}
	return Bounds().Z( e, tau,s, E_arg);
}

double integrand2(double x, void* prms) {
	Integrand2_prms* p = (Integrand2_prms*)prms;
	//if ((p->sgn == 1 && p->R.bound == '>') || (p->sgn == -1 && p->R.bound == '<')) {
	//	fprintf(stderr, "Calculation::integrand2 : wrong regime argument. sgn = %hd, bound = %c\n", p->sgn, p->R.bound);
	//	exit(1);
	//}

	char regime = p->R.regime; // '-' , '+' , 'o' , '0' , '1'

	double zeta_ast;
	if (regime == '0' || regime == '1')
		zeta_ast = (regime == '0' ? 0e0 : 1e0);
	else
		zeta_ast = f_zeta(p->s, p->e, p->tau, regime);
	double zeta = zeta_ast + p->sgn * sqr_(x);
	double d_zeta_0 = zeta, d_zeta_1 = zeta - 1e0, d_zeta_p = zeta - f_zeta(p->s, p->e, p->tau, '+'),
		d_zeta_m = zeta - f_zeta(p->s, p->e, p->tau, '-'), d_zeta_o = zeta - f_zeta(p->s, p->e, p->tau, 'o');
	switch (regime) {
	case '0': {
		d_zeta_0 = p->sgn * x * x;
		break;
	}
	case '1': {
		d_zeta_1 = p->sgn * x * x;
		break;
	}
	case 'o': {
		d_zeta_o = p->sgn * x * x;
		break;
	}
	case '-': {
		d_zeta_m = p->sgn * x * x;
		break;
	}
	case '+': {
		d_zeta_p = p->sgn * x * x;
		break;
	}
	default:
		fprintf(stderr, "Calculation::integrand2 : wrong regime argument = %c%c\n", p->R.regime, p->R.bound);
		exit(1);
		break;
	}
	double a = -p->tau * (1 + 2 * p->tau) * sqr_(1e0 - fabs(p->s)),
		Kp = sqrt(fabs(p->tau - p->s) + sqrt(sqr_(p->tau - p->s) + a * d_zeta_p)),
		Km = sqrt(fabs(p->tau + p->s) + sqrt(sqr_(p->tau + p->s) + a * d_zeta_m)),
		sqrt_root = sqrt(sqr_(1e0 + fabs(p->s)) - sqr_(1e0 - fabs(p->s)) * zeta),
		tmp;
	unsigned the_case = (p->s < -fabs(p->tau) ? 0 : (p->tau < 0e0 ? 1 : 2));
	switch (the_case) {
	case 0: { // tau < -abs(s)
		tmp = Kp * Km / sqrt(-d_zeta_m * d_zeta_p) / sqrt(a * d_zeta_o);
		//		printf("a = %f, d_zeta_p = %f, d_zeta_m = %f, d_zeta_0 = %f\n", a, d_zeta_p, d_zeta_m, d_zeta_o);
		//		printf("0: %e, 1: %e, 2: %e\n",Kp*Km, sqrt(-d_zeta_m*d_zeta_p), sqrt(a*d_zeta_o) );
		break;
	}
	case 1: { // -abs(s) < tau < 0
		tmp = Km / Kp / sqrt(-d_zeta_o * d_zeta_m);
		break;
	}
	case 2: { // 0 < tau
		tmp = Kp / Km / sqrt(-d_zeta_o * d_zeta_p);
	}
	}
	double result = tmp * fabs(p->tau) / fabs(a) / sqrt_root / sqrt(-d_zeta_0 * d_zeta_1) / pow(3.14159, 3e0);
	/*if (result != result) {
		static unsigned n_err = 0;
		FILE* file = fopen("integrand2.err", "at");
		fprintf(file, "n_err = %u: regime = %c, sgn = %hd\n", ++n_err, regime, p->sgn);
		fprintf(file, "e = %+.10f\n", p->e);
		fprintf(file, "\tzeta = %+.12f\n", zeta);
		fprintf(file, "\td_zeta_0 = %e\n", d_zeta_0);
		fprintf(file, "\td_zeta_1 = %e\n", d_zeta_1);
		fprintf(file, "\td_zeta_p = %e\n", d_zeta_p);
		fprintf(file, "\td_zeta_m = %e\n", d_zeta_m);
		fprintf(file, "\td_zeta_o = %e\n", d_zeta_o);
		fprintf(file, "\ttmp = %e\n", tmp);
		fprintf(file, "\tKp = %e, Km = %e\n", Kp, Km);
		fprintf(file, "\tsqrt_root = %e\n", sqrt_root);
		fclose(file);
	}*/
	double J = x;
	/*
		double R_val = R(p->e, p->tau, p->s, zeta),
		phi_val = phi((sqrt(R_val) + p->s)/p->tau);
		double result_old = phi_val/sqrt(R_val)/sqrt_root/sqrt(zeta*(1e0 - zeta))/pow(M_PI, 3e0);
		printf("s = %f, zeta = %e, result = %e, result_old = %e\n", p->s, zeta, result, result_old);
		// getchar();
	*/
	return J * result;
}


#endif
