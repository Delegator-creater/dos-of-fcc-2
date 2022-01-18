


#include "out_err.h"
#include "Macros.h"

#include "Stopwatch.h"
#include "All_prm.h"
#include "Integr.h"
#include <cstring>
#include <deque>
#include <string>
#include <ctime>
#include <thread>
#include <functional>
#include <algorithm>
#include "Func.h"
#include "spt_func.h"
#include <array>
#include <map>
#include "Bounds.h"
#include <set>
#include <thread>
#include "Interval.h"

std::string prm_start;

using namespace std::placeholders;

Stopwatch stopwatch_f;
double f(double x , void * param) {

	stopwatch_f.start();

	auto prm = (all_prm *)param;
	

	double Q;


		if (prm->version == 0) {

			Q = prm->func(prm->e, prm->t, prm->s, prm->v_db + x);
			Q += prm->func(prm->e, prm->t, prm->s, prm->v_ub - x);
		}
		else {
			Integrand2_prms new_prm;
			new_prm.e = prm->e + 6. * prm->t;
			new_prm.tau = prm->t;
			new_prm.s = prm->s;

			new_prm.sgn = +1;
			new_prm.R.regime = prm->c_db->type_bounds;
			new_prm.R.bound = prm->v_db;
			Q = integrand2(x, &new_prm);

			new_prm.sgn = -1;
			new_prm.R.regime = prm->c_ub->type_bounds;
			new_prm.R.bound = prm->v_ub;
			Q += integrand2(x, &new_prm);

		}


	stopwatch_f.end();
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

Stopwatch stch_first_integral;
template<const int prm_comp>
double first_integral(double s , void * parametr_ ) {

	stch_first_integral.start();

	auto parametr = (all_prm*)parametr_;



	parametr->s = s;
	double t = parametr->t;
	double e = parametr->get_e();

	double bl;
	double bu;

	switch (prm_comp){
		case 0: {
			bl = parametr->down(e, t, s);
			bu = parametr->up(e, t, s);
			if (bu > 1.) {
				bu = 1.;
				parametr->c_ub->type_bounds = '1';
			}
			if (bl < 0.) {
				bl = 0.;
				parametr->c_db->type_bounds = '0';
			}
			if (bu <= bl) {
				stch_first_integral.end();
				return 0.;
			}
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


	double bounds[2] = {0 , (bu - bl) / 2. };

	parametr->v_db = bl;
	parametr->v_ub = bu;


	parametr->type_bounds_down = parametr->c_db->type_bounds;
	parametr->type_bounds_up = parametr->c_ub->type_bounds;

	int status;
	result = integration<1>(&f,
								parametr_,
								bounds[0],
								bounds[1],
								parametr->error_in1,
								parametr->error_in2,
								status
								);


	
	if (status != 0) {
		if (list_cods_status.find(status) == list_cods_status.end()) {
			list_cods_status.insert(std::pair<int , std::deque<data_for_cod_status*> >(status, std::deque<data_for_cod_status*>()));
		auto elem_list = list_cods_status.find(status);
		data_for_cod_status* data = new data_for_cod_status;
		data->e = e;
		data->t = t;
		data->s = s;
		data->num_interg = 2;
		elem_list->second.push_back(data);}

	}

	//file_logs_f_int.close();


	/*all_prm::vfi data;
	data.int_for_zeta = all_res;
	data.s = s;
	data.zeta1 = bl[0];
	data.zeta2 = bu[0];
	parametr->value_first_int.push_back(data);*/
	stch_first_integral.end();
	return result;
}

Stopwatch stch_second_integral;
template<const int prm_int>
double second_integral(all_prm & a_prm  ) {

	stch_second_integral.start();

	
	double t = a_prm.t;

	double abserr;
	double result=0;
	int status;

	double bounds_s_for_tau = std::min(std::abs(t), 1.);
	
	auto inter_s_for_tau_less = Interval(-bounds_s_for_tau,  bounds_s_for_tau);
	auto inter_s_for_tau_more = Interval(-1.			  , -bounds_s_for_tau);

	auto s_bounds_d_old_less   = a_prm.c_db->get_s_bounds(t, a_prm.e , false)        * inter_s_for_tau_less;
	auto s_bounds_d_old_more   = a_prm.c_db->get_s_bounds(t ,a_prm.e , true )        * inter_s_for_tau_more;
	auto s_bounds_u_less       = a_prm.c_ub->get_s_bounds(t, a_prm.e , false)        * inter_s_for_tau_less;
	auto s_bounds_u_more	   = a_prm.c_ub->get_s_bounds(t, a_prm.e , true )        * inter_s_for_tau_more;

	static auto operator_inter = [](const Interval & inter1 , const Interval & inter2 , const Interval & total_set) {
		auto result = (inter1 / inter2) + (inter2 / inter1) + (inter1 * inter2);
		return total_set / result + result;
	};
	Interval s_bounds[2];
	s_bounds[0] = operator_inter(s_bounds_d_old_less, s_bounds_u_less , inter_s_for_tau_less);
	s_bounds[1] = operator_inter(s_bounds_d_old_more, s_bounds_u_more,  inter_s_for_tau_more);

	/*
	Interval s_bounds_d_less_p;
	Interval s_bounds_d_less_m;
	Interval s_bounds_d_more;
	for (auto& i : (*s_bounds_d_old.get_one_ints()))
	{
		s_bounds_d.add_O_I(i);
		Interval tmp_inter( -i.get_bounds()[1] , -i.get_bounds()[0] );
		s_bounds_d.add_O_I(tmp_inter.get_one_ints()->at(0));
	}*/

	//for (size_t i = 0; i < 2; ++i) {
	//	if (!s_bounds->empty())
	//		for (auto & bounds : s_bounds->get_bounds())
				result += integration<2>(&first_integral<prm_int>,
										&a_prm,
										-1,
										bounds_s_for_tau,
										a_prm.error_ext1,
										a_prm.error_ext2,
										status
										);
	//}
	
	if (status != 0) {
		if (list_cods_status.find(status) == list_cods_status.end()){
			list_cods_status.insert(std::pair<int, std::deque<data_for_cod_status*> >(status, std::deque<data_for_cod_status*>() ));
		auto elem_list = list_cods_status.find(status);
		data_for_cod_status* data = new data_for_cod_status;
		data->e = a_prm.get_e();
		data->t = t;
		data->num_interg = 1;
		elem_list->second.push_back(data);
		}
	}
	//if (status != 0)
	//	err_ << "GSL status ext_int: " << gsl_strerror(status) <<'\n';


	//result = integration(first_integral<prm_int>, &a_prm, -1, std::min(std::std::abs(t), 1.), a_prm.error, a_prm.error_);
	stch_second_integral.end();
	return result;
}

double rho(double t, double e , all_prm & parametr , int prm_ ) {
	try{


		parametr.f = f;
		parametr.t = t;
		parametr.e = e;
		parametr.size_memory_1 = 100000;
		parametr.size_memory_2 = 10000;
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
		//out <<t <<" " << e << " "<<  result << " || t = " <<(clock() - start) / CLOCKS_PER_SEC << " sec\n";
		//mtx3.unlock();
		return result;
	}
	catch(...){
		err_ << "error in function \"rho\"";
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
			const bool error_write = false,
			const int s = 0) {

	using namespace std;
	
	all_prm prm(t);
	prm.error_ext1  = error_ext1;
	prm.error_ext2  = error_ext2;
	prm.error_in1   = error_in1;
	prm.error_in2   = error_in2;
	prm.version     = version;
	prm.prm_progr   = prm_progr;
	
	const int prm_progr_ = (prm_progr == 1) ? 1 : 
						   (prm_progr == 3) ? 3 : 0;


	UpBounds ub;
	DownBounds db;

	prm.c_ub = &ub;
	prm.c_db = &db;
	prm.up = ub.get_bounds(t); 
	prm.down = db.get_bounds(t);


	double min_arg;
	double max_arg ;

	using type_main_func = std::function<double(double, double, all_prm&, int)>;
	type_main_func main_function;
	double step_for_second_arg;
	switch (prm_progr_)
	{
	case 0 : {
		step_for_second_arg = step;
		main_function = type_main_func(rho);
		min_arg = (prm_progr == 2 ) ? e : min(-4.-2.*t , min(-4.-6.*t, 6.0*t ));
		max_arg = max(12.0 - 6.0*t, 6.0*t);
		break; }
	case 1: {


		step_for_second_arg = step;
		main_function = type_main_func
			([=](double, double s, all_prm& prm_, int type_bounds_) {
			if (type_bounds_ == 0)
				return first_integral<0>(s, &prm_);
			else if (type_bounds_ == 1)
				return first_integral<1>(s, &prm_);
				});
		prm.f = f;
		prm.t = t;
		prm.e = e - 6*t;
		prm.size_memory_1 = 1000000;
		prm.size_memory_2 = 1000000;
		min_arg = -1;
		max_arg = std::min(std::abs(t), 1.);
		break;
	}
	case 3: {
		step_for_second_arg = step;
		main_function = type_main_func([&](double, double zeta, all_prm& prm_, int type_bounds_){
			double bl = prm_.down(prm_.e, t, s);
			double bu = prm_.up(prm_.e, t, s);

			prm_.type_bounds_down = prm_.c_db->type_bounds;
			prm_.type_bounds_up   = prm_.c_ub->type_bounds;
			if (bu > 1.) {
				bu = 1.;
				prm_.c_ub->type_bounds = 5;
			}
			if (bl < 0.) {
				bl = 0.;
				prm_.c_db->type_bounds = 4;
			}
			if (bu <= bl) {
				
				return 0.;
			}
			prm_.v_db = bl;
			prm_.v_ub = bu;
			return f(zeta , (void*)(&prm_) );
			});
		prm.f = f;
		prm.t = t;
		prm.e = e - 6 * t;
		prm.s = s;
		prm.size_memory_1 = 1000000;
		prm.size_memory_2 = 1000000;
		min_arg = 0.;
		max_arg = 1.;

		break;
	}
		
	}

	out << std::string("start for t =") + std::to_string(t)  +  std::string("\n") ;

	prm.value_first_int;

	double start = clock();
	out << "step_for_second_arg " << step_for_second_arg << '\n';
	for (double second_arg =min_arg; second_arg < max_arg; second_arg += step_for_second_arg) {
		//flagu = flagl = false;
		double chang_second_arg = prm_progr_ == 0 ? second_arg - 6. * t : second_arg;
		double res = main_function(t, chang_second_arg, prm, type_bounds);
		out << second_arg << " \t";
			out << res  << '\t';

			if (error_write) {
				std::map<int, std::set<int> > type_integration;
				for (auto& cods : list_cods_status)
					for (auto status : cods.second)
						if (type_integration.find(cods.first) == end(type_integration))
							type_integration.insert(std::pair<int, std::set<int>>(cods.first, std::set<int>({ status->num_interg })));
						else
							type_integration.find(cods.first)->second.insert(status->num_interg);
				for (auto& cods : list_cods_status) {
					out << gsl_strerror(cods.first) << " ";
					const auto f_elem = type_integration.find(cods.first);
					if (f_elem != type_integration.end())
						for (auto i : f_elem->second)
							out << i << " ";
					out << "\t";
				}
			}
			out << " \n";
			list_cods_status.clear();
			prm.value_first_int.clear();

			err_ << (second_arg - min_arg) / (max_arg - min_arg) * 100 << "% " << (clock() - start) / CLOCKS_PER_SEC << " sec\n";
			stopwatch_f.print(err_, "stopwatch_f");
			stch_Func_Q.print(err_, "stch_Func_Q");
			stch_first_integral.print(err_, "stch_first_integral");
			stch_second_integral.print(err_, "stch_second_integral");
			
			stopwatch_f.reset();
			stch_Func_Q.reset();
			stch_first_integral.reset();
			stch_second_integral.reset();

			out.close();
			out.open(name_file_out, std::ios_base::app);

			err_.close();
			err_.open(name_file_err_, std::ios_base::app);

		if (prm_progr == 2)
			break;
	}

}

// ===========================================================================
// Формат запуска файла name_programm parametr_programm name_file_parametr
// ===========================================================================
// name_programm -- название программы
// ===========================================================================
// parametr_programm -- параметр запуска программы. Определяют поведение программы при запуске.
// "-с"    -- расчета плотности состояний
// "-f"	   -- расчет первого интеграла плотности состояния
// "-c0"   --
// "-func" --
// "-old"  --
// "-new"  --
//"
// ===========================================================================
// name_file_parametr -- название текстового файла с нужными данными для работы программы.
// Данные определяются параметром запуска программы.
// ===========================================================================
// При запуске программы данные будут записываться в файл, в названии которого будет соответсвующий префикс.
// ===========================================================================
#include "Interval.h"
int main(int argc, char *argv[]) {

	//Interval i1 = Interval(-1,1) / Interval(-0.8,0.8);
	//for (auto i : i1.get_bounds())
	//	std::cout << i[0] << " " << i[1] << '\n';



	using namespace std;
 	std::cout.setf(ios::scientific);
	std::cout.precision(3);
	/*std::cout << string("start ")
				<< to_string(argc)
				<< string(" ")
				<< string( argv[0] ) // path file programm
				<< string(" ")
				<< string(argv[1]) // parametrs calculate
				<< string(" ")
				<< string( argv[2] ) // version
				<< string(" ")
				<< string(argv[3])  // name file parametrs
				<< string(" ")
				<< string(argv[4])//path OUT
				<< string( "\n" ) ;*/
	
	std::set<std::string> set_prm_start;
	for (size_t i = 1; i < argc; ++i) {
		set_prm_start.insert(std::string(argv[i]));
	}

	int prm_progr = 0;
	int ver       = 0;
	
	gsl_set_error_handler_off();

	
		std::deque<double> arr_double;

		std::ifstream param_file(std::string("./") + std::string(argv[3]), std::ios::in);
		std::string line;
		if (param_file.is_open())
			while (getline(param_file, line))
				arr_double.push_back(std::stod(line));
		else {
			std::cerr << "Error, not found file!\n" ;
			return -2;
		}

		param_file.close();
		std::string postfix_name_file("");
		
		for (auto i : arr_double)
			postfix_name_file += std::to_string(i) + std::string("_");
		postfix_name_file += std::string(argv[1]);
		postfix_name_file += std::string(argv[2]);

		name_file_out = std::string(argv[4]) + std::string("//") + std::string("OUT_") + postfix_name_file;
		out.open(name_file_out);

		name_file_err_ = std::string(argv[4]) + std::string("//") + std::string("ERROR_") + postfix_name_file;
		err_.open(name_file_err_);


		out.setf(ios::scientific);
		//out.precision(3);
		err_.setf(ios::scientific);
		err_.precision(3);


		const auto progr_bind = std::bind(progr, 
			arr_double[0], // tau
			arr_double[1], // abs1
			arr_double[2], // rel1
			arr_double[3], // abs2
			arr_double[4], // rel2
			arr_double[5], // epsilon (for -f)
			arr_double[6], // type bounds
			arr_double[7], // step for second arg
			_1 ,		
			_2 , 
			bool(arr_double[8]),// logs gsl error (true, false)
			arr_double[9] // s (for -fun)
		);

		prm_start = std::string(argv[1]);

		const auto find_in_set = [&](const char c[]) {
			bool res = set_prm_start.find(string(c)) != end(set_prm_start);
			if (res)
				out << "start " << c << "\n";
			return res;
		};

		if ( find_in_set("-c") ) 
			prm_progr = 0;

		if ( find_in_set("-f") )
			prm_progr = 1;
			
		if ( find_in_set("-c0"))
			prm_progr = 2;
		
		if ( find_in_set("-fun"))
			prm_progr = 3;

		if ( find_in_set("-old")) 
			ver = 0;
		
		if ( find_in_set("-new"))
			ver = 1;
		

		progr_bind(prm_progr , ver);
		out.close();
		err_.close();




	return 0;

}

