#include <cstring>
#include <fstream>
#include <deque>
#include <string>
#include <iostream>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>
#include <ctime>
#include <thread>
#include <functional>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include "Func.h"
#include "spt_func.h"
#include <array>
#include "Integr.h"
std::string prm_start;

using namespace std::placeholders;

class all_prm  {
	double e;

public:
	all_prm(double t) : func(t , prm_start) {}

	double get_e() { return e; }
	void set_e(double new_e) { e = new_e; }

	double t;
	double s;
	double error;
	double error_;



	double(*f)(double x, void * param);

	std::function<double(double , double , double)> up;
	std::function<double(double, double, double)> down;

	Func func;
	
	gsl_integration_workspace * memory1;
	gsl_integration_workspace * memory2;

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

using fd_2x = std::function<double(double, double)>;
using fd_3x = std::function<double(double, double, double)>;

class Bounds{
protected:

	inline double Z(double e, double t, double s, double f) {
		return ((f - e) / 4. + (1. - 2. * t)* sqr_(s))
			/ ((1. + 2. * t)* sqr_(1. - std::abs(s)));
	}

	inline double s1(ARG_3, double sign) {
		double a = e /
			(1. - 2.*t + 1./t);
		if (a > 0.)
			return sign * 0.50* std::sqrt(a);
		else
			return 0.;
	}

	inline double s2(ARG_3, double sign, double sigm) {
		auto a = 1. - (1. - 2. * t)*(t - e / 4.);
		if (a > 0)
			return -sigm / (1. - 2. * t) + sign * std::sqrt(a)
			/ std::abs(1. - 2. * t);
		else
			return 0;
	}

	inline double u1(ARG_3, double sign) {
		auto a = -4 * t * t + 1 / 4 / t * (e + 4) - t * e + 3;
		if (a > 0.)
			return -t * (1. + 2. * t) / (1. - 4. * sqr_(t)) + sign*std::abs(t) / (1. - 4. * t*t)* std::sqrt(a);
		else
			return 0;
	}

	inline double u2(ARG_3, double sign, double sigm) {
		double a = 2. - e * t + 2. * sigm  * (1. + 2. * t);
		if (a > 0)
			return (1. + 2. * t + sigm) / (4. * t) + sign * std::sqrt(a)
			/ std::abs(4. * t);
		else
			return 0;
	}


};

class UpBounds final : public Bounds {
public:
	UpBounds() : Bounds() {}

	

	//*******************************//

	/// Z(e_1 , t , s ,sigma_ * 8s +4t)
	double bounds1(ARG_3 , double sigm_) {
		double sigm = double(SIGN(s)) * sigm_;
		auto return_value1 = std::function<double(void)>( [=]() {return double(1); });
		auto return_value2 = std::function<double(void)>( [=]() {return Z(e, t, s, sigm_*8. * s + 4. * t); });
		if (t > -0.5 && t < 0.)
			std::swap(return_value1 , return_value2);
		if (inter(std::abs(s), u2(e, t, s, -1., sigm), u2(e, t, s, 1., sigm)))
			return return_value1();
		else 
			return return_value2();
		
	}
	//*******************************//

	///Z(e_1 , t , s , -4s^2/t)
	double bounds2(ARG_3) {
		auto return_value1 = std::function<double(void)>([=]() {return double(1); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, 4.* sqr_(s)/t); });
		if ( (0.< t) && (t < 0.5))
			std::swap(return_value1, return_value2);
		if (inter(std::abs(s), u1(e, t, s, -1.), u1(e, t, s, 1.)))
			return return_value1();
		else
			return return_value2();
	}

	//*******************************//
	fd_3x get_bounds(const double t_) {
		if (t_ >= 0) {
			return [=](ARG_3) {return (std::abs(s) > std::abs(t)) ? this->bounds1(e, t, s, -1.) : this->bounds2(e, t, s); };
			// return Z(... , -8s - 4t) or Z(... , 4s^2/t)
		}
		else if (t_ >= -0.5) {
			return [=](ARG_3) {return  this->bounds1(e, t, s, -1.); };
			// return Z(... , -8s - 4t)
		}
		else {
			return [=](ARG_3) {return (std::abs(s) > std::abs(t)) ? this->bounds1(e, t, s, 1.) : this->bounds2(e, t, s); };
			// return Z(... , 8s - 4t) or Z(... , 4s^2/t)
		}
	}



	
};

class DownBounds final : public Bounds {
public:
	DownBounds() : Bounds() {}

	//*******************************//

	/// Z(e_1 , t , s ,sigma * 8s +4t)
	double bounds1(ARG_3 , double sigm) {
		auto return_value1 = std::function<double(void)>([=]() {return double(0); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, sigm*8. * s + 4. * t); });
		if ( std::abs(t) > 0.5 )
			std::swap(return_value1, return_value2);
		if ( (s2(e, t, s, -1., sigm) < s) && (s < s2(e, t, s, 1., sigm)))
			return return_value2();
		else
			return return_value1();
	}

	//*******************************//

	///Z(e_1 , t , s , -4s^2/t)
	double bounds2(ARG_3) {
		auto return_value1 = std::function<double(void)>([=]() {return double(0); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, 4. * sqr_(s) / t); });
		if ((0. < t) && (t < 1.) )
			std::swap(return_value1, return_value2);
		if (inter(s, s1(e, t, s, -1.), s1(e, t, s, 1.)))
			return return_value2();
		else
			return return_value1();
	}

	//*******************************//


	

	fd_3x get_bounds(const double t_) {
		if (t_ >= 0.) {
			return [=](ARG_3) {return this->bounds1(e, t, s, 1.); };
		}// return Z(... , 8s - 4t)
		else if (t_ >= -0.5) {
			return [=](ARG_3) {return ((std::abs(s) > std::abs(t)) ? this->bounds1(e, t, s, 1.) : this->bounds2(e, t, s)); };
		}// return Z(... , 8s - 4t) or Z(... , 4s^2/t)
		else {
			return [=](ARG_3) {return this->bounds1(e, t, s, -1.); };
		}// return Z(... , -8s - 4t) 
	}

};

double f(double x , void * param) {
	auto prm = (all_prm *)param;
	return prm->func(prm->get_e(), prm->t, prm->s, x);
}

bool flagl = false;
bool flagu = false;

template<const int prm_comp>
double first_integral(double s , void * parametr_ ) {

	auto parametr = (all_prm*)parametr_;

	static std::ofstream file_logs_f_int("logs_f_int");
	file_logs_f_int.precision(3);
	const double eps_xi = 0.01;

	parametr->s = s;
	double t = parametr->t;
	double e = parametr->get_e();

	double bl[2];
	double bu[2];


	bl[0] = parametr->down(e, t, s);
	bu[0] = parametr->up(e, t, s);
	if (bu[0] > 1)
		bu[0] = 1;
	if (bl[0] < 0)
		bl[0] = 0;
	if (bu[0] < bl[0])
		return 0;

	bl[0] = std::sqrt(bl[0]);
	bu[0] = std::sqrt(bu[0]);

	bl[1] = 0.0;
	bu[1] = 1.0;

	double** arr_bounds = new double* [3];
	for (size_t i = 0; i < 3; ++i)
		arr_bounds[i] = new double[2];
	arr_bounds[0][0] = bl[1];
	arr_bounds[0][1] = bl[0];
	arr_bounds[1][0] = bl[0];
	arr_bounds[1][1] = bu[0];
	arr_bounds[2][0] = bu[0];
	arr_bounds[2][1] = bu[1];

	flagu = true;
	flagl = false;


	double result[3] = {0,0,0};
	double abserr;

	gsl_function f_integral;
	f_integral.function = parametr->f;
	f_integral.params = parametr_;

	bool coundit_bounds[3];
	coundit_bounds[0] = std::abs(bl[0]- 0) > 0.001;
	coundit_bounds[1] = true;
	coundit_bounds[0] = std::abs(bu[0]- 1 ) >0.001;

	for (size_t i = 0 ; i < 3 ; ++i)
		if (coundit_bounds[i])
			gsl_integration_qag(&(f_integral),
				arr_bounds[i][0], 
				arr_bounds[i][1],
				parametr->error,
				parametr->error_,
				parametr->size_memory_2 ,
				6,
				parametr->memory2,
				&result[i],
				&abserr);

	if ((result[0] != 0) || (result[2] != 0)) {
		file_logs_f_int
			<< e + 6. * parametr->t << '\t'
			<< parametr->t << '\t'
			<< s << '\t'
			<< bl[0] << '\t'
			<< bu[0] << '\t'
			<< result[0] << '\t'
			<< result[1] << '\t'
			<< result[2] << '\n';
	}

	//file_logs_f_int.close();


	double all_res = 0;
	for (size_t i = 0; i < 3; ++i) {
		delete[] arr_bounds[i];
		all_res += result[i];
	}
	delete[] arr_bounds;

	all_prm::vfi data;
	data.int_for_zeta = all_res;
	data.s = s;
	data.zeta1 = bl[0];
	data.zeta2 = bu[0];
	parametr->value_first_int.push_back(data);
	return all_res;
}
template<const int prm_int>
double second_integral(all_prm & a_prm  ) {

	

	gsl_function f_integral;
	f_integral.function = first_integral<prm_int>;
	f_integral.params = &a_prm;
	
	double t = a_prm.t;

	double abserr;
	double result;
	gsl_integration_qag(&f_integral,
			-1.,
			std::min(1., std::abs(t)),
			a_prm.error,
			a_prm.error_,
			a_prm.size_memory_1,
			6,
			a_prm.memory1,
			&result,
			&abserr);



	//result = integration(first_integral<prm_int>, &a_prm, -1, std::min(std::std::abs(t), 1.), a_prm.error, a_prm.error_);

	return result;
}

double rho(double t, double e , all_prm & parametr , int prm_ ) {
	try{


		parametr.f = f;
		parametr.t = t;
		parametr.set_e(e);
		parametr.size_memory_1 = 100000;
		parametr.size_memory_2 = 100000;
		parametr.memory1 = gsl_integration_workspace_alloc(parametr.size_memory_1);
		parametr.memory2 = gsl_integration_workspace_alloc(parametr.size_memory_2);
		double result;
		switch (prm_){
			case 0 :{
				result = second_integral<0>(parametr);
				break;}
			case 1 :{
				result = second_integral<1>(parametr);
				break;} 
			default:{
				throw -1;
				}
			
		}
		//mtx3.lock();
		//std::cout <<t <<" " << e << " "<<  result << " || t = " <<(clock() - start) / CLOCKS_PER_SEC << " sec\n";
		//mtx3.unlock();
		return result;
	}
	catch(...){
		std::cerr << "error in function \"rho\"";
	}	 
}

template<const int prm_progr = 0>
void progr(const double t , const double error1 , const double error2 , const double e = 0  , int type_bounds = 0 
	, double step=0.1) {

	using namespace std;
	
	all_prm prm(t);
	prm.error  = error1;
	prm.error_ = error2;
	   	 

	UpBounds ub;
	DownBounds db;

	prm.up = ub.get_bounds(t); 
	prm.down = db.get_bounds(t);


	double min_arg;
	double max_arg ;
	std::function<double(double, double, all_prm & , int)> main_function;
	double step_for_second_arg;
	switch (prm_progr)
	{
	case 0 : {
		step_for_second_arg = step;
		main_function = std::function<double(double, double, all_prm & , int)>(rho);
		min_arg = min(-4.-2.*t , min(-4.-6.*t, 6.0*t ));
		max_arg = max(12.0 - 6.0*t, 6.0*t);
		break; }
	case 1 : {
		step_for_second_arg = step;
		main_function = std::function<double(double, double s ,  all_prm & , int)>([&](double, double s, all_prm & prm , int) {
			return first_integral<0>(s, &prm);
			});
		prm.f = f;
		prm.t = t;
		prm.set_e(e);
		prm.size_memory_1 = 100000;
		prm.size_memory_2 = 100000;
		prm.memory1 = gsl_integration_workspace_alloc(prm.size_memory_1);
		prm.memory2 = gsl_integration_workspace_alloc(prm.size_memory_2);
		min_arg = -1;
		max_arg = std::min(std::abs(t), 1.);
		break;
		}
	}

	std::cout << std::string("start for t =") + std::to_string(t)  +  std::string("\n") ;

	prm.value_first_int;

	double start = clock();
	cout << "step_for_second_arg" << step_for_second_arg << '\n';
	for (double second_arg =min_arg; second_arg < max_arg; second_arg += step_for_second_arg) {
		flagu = flagl = false;
		cout << second_arg << " \t";
		double res = main_function(t, second_arg - 6. * t, prm, type_bounds);
		cout << res << " \t";
		cout << flagu ? 1 : 0;
		cout << " \t";
		cout << flagl ? -1 : 0;
		cout << " \t";
		for (all_prm::vfi& i : prm.value_first_int)
			cout << i.s << " \t" << i.int_for_zeta << " \t" << i.zeta1 << " \t" << i.zeta2 << " \t";
		cout << " \n";
		prm.value_first_int.clear();
		
	}
}

// ===========================================================================
// Формат запуска файла name_programm parametr_programm name_file_parametr
// ===========================================================================
// name_programm -- название программы
// ===========================================================================
// parametr_programm -- параметр запуска программы. Определяют поведение программы при запуске.
// "-с" -- расчета плотности состояний
// "-f" -- расчет первого интеграла плотности состояния
// ===========================================================================
// name_file_parametr -- название текстового файла с нужными данными для работы программы.
// Данные определяются параметром запуска программы.
// ===========================================================================
// При запуске программы данные будут записываться в файл, в названии которого будет соответсвующий префикс.
// ===========================================================================
int main(int argc, char *argv[]) {

	using namespace std;
	cout.setf(ios::scientific);
	cout.precision(3);
	cout << string("start ") 
		<< to_string(argc)
		<< string(" ")
		<< string( argv[0] )
		<< string(" ")
		<< string(argv[1])
		<< string(" ")
		<< string( argv[2] )
		<< string( "\n" ) ;

	
	gsl_set_error_handler_off();

	if (argc == 3 ) {
		

		std::deque<double> arr_double;

		std::ifstream param_file(std::string(argv[2]) + std::string(".txt"), std::ios::in);
		std::string line;
		if (param_file.is_open())
			while (getline(param_file, line))
				arr_double.push_back(std::stod(line));
		else {
			std::cout << "Error, not found file!\n" ;
			return -2;
		}

		param_file.close();

		prm_start = std::string(argv[1]);

		if ( std::strcmp( argv[1] , "-c") == 0 ) {
			std::cout << "start -c\n";
			progr(arr_double[0], arr_double[1], arr_double[2] , arr_double[3] , arr_double[4] ,arr_double[5] );
		}
		if ( std::strcmp( argv[1] , "-f" ) == 0 ) {
			std::cout << "start -f\n";
			progr<1>(arr_double[0], arr_double[1], arr_double[2], arr_double[3] , arr_double[4] ,arr_double[5]);
		}
		
	}


	return 0;

}

