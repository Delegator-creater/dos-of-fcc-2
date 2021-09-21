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
#include <map>
#include "Bounds.h"
#include "All_prm.h"
#include <set>

std::string prm_start;

using namespace std::placeholders;



double f(double x , void * param) {
	auto prm = (all_prm *)param;
	

	double Q;
	try {

		if (prm->version == 0) {

			Q = prm->func(prm->e, prm->t, prm->s, prm->v_db + x);
			Q += prm->func(prm->e, prm->t, prm->s, prm->v_ub - x);
		}
		else {
			switch (prm->type_bounds_down)
			{
				case 1 : {
					Q = prm->func.Q<1>(prm->e, prm->t, prm->s, x, +1., prm->v_db);
					break;
				}
				case 2: {
					Q = prm->func.Q<2>(prm->e, prm->t, prm->s, x, +1., prm->v_db);
					break;
				}
				case 3 : {
					Q = prm->func.Q<3>(prm->e, prm->t, prm->s, x, +1., prm->v_db);
					break;
				}
				case 4: {
					Q = prm->func.Q<4>(prm->e, prm->t, prm->s, x, +1., prm->v_db);
					break;
				}
			}

			switch (prm->type_bounds_up)
			{
				case 1: {
					Q += prm->func.Q<1>(prm->e, prm->t, prm->s, x, -1., prm->v_ub);
					break;
				}
				case 2: {
					Q += prm->func.Q<2>(prm->e, prm->t, prm->s, x, -1., prm->v_ub);
					break;
				}
				case 3: {
					Q += prm->func.Q<3>(prm->e, prm->t, prm->s, x, -1., prm->v_ub);
					break;
				}
				case 5: {
					Q += prm->func.Q<5>(prm->e, prm->t, prm->s, x, -1., prm->v_ub);
					break;
			}
			}

		}

	}
	catch (int) {
		std::cerr << "Error in function f(doube , void*)\n";
		exit(1);
		//return 0.; 
	}

	return Q;
}

bool flagl = false;
bool flagu = false;

using data_for_cod_status = struct {
	double e;
	double t;
	double s;
	int num_interg;
};

std::map<int, std::deque<data_for_cod_status*> > list_cods_status;

template<const int prm_comp>
double first_integral(double s , void * parametr_ ) {

	auto parametr = (all_prm*)parametr_;



	parametr->s = s;
	double t = parametr->t;
	double e = parametr->get_e();

	double bl;
	double bu;

	switch (prm_comp)
	{
	case 0: {
		bl = parametr->down(e, t, s);
		bu = parametr->up(e, t, s);
		if (bu > 1) {
			bu = 1;
			parametr->c_ub->type_bounds = 5;
		}
		if (bl < 0) {
			bl = 0;
			parametr->c_db->type_bounds = 4;
		}
		if (bu < bl)
			return 0;
		break;
	}
	case 1: {
		bl = 0.0;
		bu = 1.0;
	}
	}




	flagu = false;
	flagl = false;


	double result;
	double abserr;

	gsl_function f_integral;
	f_integral.function = parametr->f;
	f_integral.params = parametr_;

	double bounds[2] = {0 , (bu - bl) / 2. };

	parametr->v_db = bl;
	parametr->v_ub = bu;


	parametr->type_bounds_down = parametr->c_db->type_bounds;
	parametr->type_bounds_up = parametr->c_ub->type_bounds;

	auto status = gsl_integration_qagp(&(f_integral),
				bounds,
				2,
				parametr->error_in1,
				parametr->error_in2,
				parametr->size_memory_2 ,
				parametr->memory2,
				&result,
				&abserr);


	
	if (status != 0) {
		if (list_cods_status.find(status) == list_cods_status.end()) 
			list_cods_status.insert(std::pair<int , std::deque<data_for_cod_status*> >(status, std::deque<data_for_cod_status*>()));
		auto elem_list = list_cods_status.find(status);
		data_for_cod_status* data = new data_for_cod_status;
		data->e = e;
		data->t = t;
		data->s = s;
		data->num_interg = 2;
		elem_list->second.push_back(data);

	}

	//file_logs_f_int.close();


	/*all_prm::vfi data;
	data.int_for_zeta = all_res;
	data.s = s;
	data.zeta1 = bl[0];
	data.zeta2 = bu[0];
	parametr->value_first_int.push_back(data);*/
	return result;
}
template<const int prm_int>
double second_integral(all_prm & a_prm  ) {

	

	gsl_function f_integral;
	f_integral.function = first_integral<prm_int>;
	f_integral.params = &a_prm;
	
	double t = a_prm.t;

	double abserr;
	double result;
	auto status = gsl_integration_qag(&f_integral,
			-1.,
			std::min(1., std::abs(t)),
			a_prm.error_ext1,
			a_prm.error_ext2,
			a_prm.size_memory_1,
			6,
			a_prm.memory1,
			&result,
			&abserr);
	
	if (status != 0) {
		if (list_cods_status.find(status) == list_cods_status.end())
			list_cods_status.insert(std::pair<int, std::deque<data_for_cod_status*> >(status, std::deque<data_for_cod_status*>() ));
		auto elem_list = list_cods_status.find(status);
		data_for_cod_status* data = new data_for_cod_status;
		data->e = a_prm.get_e();
		data->t = t;
		data->num_interg = 1;
		elem_list->second.push_back(data);

	}
	//if (status != 0)
	//	std::cerr << "GSL status ext_int: " << gsl_strerror(status) <<'\n';


	//result = integration(first_integral<prm_int>, &a_prm, -1, std::min(std::std::abs(t), 1.), a_prm.error, a_prm.error_);

	return result;
}

double rho(double t, double e , all_prm & parametr , int prm_ ) {
	try{


		parametr.f = f;
		parametr.t = t;
		parametr.e = e;
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


// version: 1 - new, 0 - old
void progr(	const double t , 
			const double error_ext1 , // abs
			const double error_ext2 , // rel
			const double error_in1,   // abs
			const double error_in2,   // rel
			const double e = 0  ,
			const int type_bounds = 0	,
			const double step=0.1 ,
			const int prm_progr = 0,
			const int version = 1,
			const bool error_write = false) {

	using namespace std;
	
	all_prm prm(t);
	prm.error_ext1  = error_ext1;
	prm.error_ext2  = error_ext2;
	prm.error_in1   = error_in1;
	prm.error_in2   = error_in2;
	prm.version     = version;

	UpBounds ub;
	DownBounds db;

	prm.c_ub = &ub;
	prm.c_db = &db;
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
		prm.e = e;
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
	cout << "step_for_second_arg " << step_for_second_arg << '\n';
	for (double second_arg =min_arg; second_arg < max_arg; second_arg += step_for_second_arg) {
		//flagu = flagl = false;
		cout << second_arg << " \t";
		double res = main_function(t, second_arg - 6. * t, prm, type_bounds);
		cout << res  << '\t';
		//cout << flagu ? 1 : 0;
		//cout << " \t";
		//cout << flagl ? -1 : 0;
		//cout << " \t";
		//for (all_prm::vfi& i : prm.value_first_int)
		//	cout << i.s << " \t" << i.int_for_zeta << " \t" << i.zeta1 << " \t" << i.zeta2 << " \t";
		if (error_write) {
			std::map<int, std::set<int> > type_integration;
			for (auto& cods : list_cods_status)
				for (auto status : cods.second)
					if (type_integration.find(cods.first) == end(type_integration))
						type_integration.insert(std::pair<int, std::set<int>>(cods.first, std::set<int>({ status->num_interg })));
					else
						type_integration.find(cods.first)->second.insert(status->num_interg);
			for (auto& cods : list_cods_status) {
				std::cout << gsl_strerror(cods.first) << " ";
				const auto f_elem = type_integration.find(cods.first);
				if (f_elem != type_integration.end())
					for (auto i : f_elem->second)
						std::cout << i << " ";
				std::cout << "\t";
			}
		}
		cout << " \n";
		list_cods_status.clear();
		prm.value_first_int.clear();
		std::cerr << (second_arg - min_arg) / (max_arg - min_arg) * 100 << "\% " << (clock() - start) / CLOCKS_PER_SEC << " sec\n";
	}

}

// ===========================================================================
// Формат запуска файла name_programm parametr_programm name_file_parametr
// ===========================================================================
// name_programm -- название программы
// ===========================================================================
// parametr_programm -- параметр запуска программы. Определяют поведение программы при запуске.
// "-с"   -- расчета плотности состояний
// "-f"	  -- расчет первого интеграла плотности состояния
// "-old" --
// "-new" --
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
		<< string(" ")
		<< string(argv[3])
		<< string( "\n" ) ;

	
	gsl_set_error_handler_off();

	if (argc == 4 ) {
		

		std::deque<double> arr_double;

		std::ifstream param_file(std::string(argv[3]), std::ios::in);
		std::string line;
		if (param_file.is_open())
			while (getline(param_file, line))
				arr_double.push_back(std::stod(line));
		else {
			std::cout << "Error, not found file!\n" ;
			return -2;
		}

		param_file.close();

		const auto progr_bind = std::bind(progr, arr_double[0], arr_double[1], arr_double[2], arr_double[3], arr_double[4], arr_double[5], arr_double[6], arr_double[7], _1 , _2 , bool(arr_double[8]));

		prm_start = std::string(argv[1]);

		int type_bounds;
		int ver;

		if ( std::strcmp( argv[1] , "-c") == 0 ) {
			std::cout << "start -c\n";
			type_bounds = 0;
			
		}
		if ( std::strcmp( argv[1] , "-f" ) == 0 ) {
			std::cout << "start -f\n";
			type_bounds = 1;
			
		}
		if (std::strcmp(argv[2], "-old") == 0) {
			std::cout << "start -old\n";
			ver = 0;
		}
		if (std::strcmp(argv[2], "-new") == 0) {
			std::cout << "start -new\n";
			ver = 1;
		}
		progr_bind(type_bounds , ver);
		
	}


	return 0;

}

