#ifndef SPT_FUNC_H
#define SPT_FUNC_H

#include <functional>
#include <string>
#include <deque>
#include <cmath>
#include "Macros.h"
#include "out_err.h"

//auto const sqr = std::bind(std::pow<double, double>, std::placeholders::_1, 2);
double sqr_(double x);

template<class Num>
inline int SIGN(Num num);

template<typename Num1 , typename Num2 , typename Num3>
inline bool cinter(Num1 x, Num2 x1, Num3 x2);

template<typename Num1 , typename Num2 , typename Num3>
inline bool inter(Num1 x, Num2 x1, Num3 x2);

template<typename T>
inline bool min(T & x1, T & x2) {
	return (x1 > x2) ? x2 : x1;
}

template<typename T>
inline bool max(T & x1, T & x2) {
	return (x1 < x2) ? x2 : x1;
}

std::deque<std::string> split(std::string  str, const char * delimiter);

template<class Num>
inline int SIGN(Num num)
{

	return num >= 0 ? 1 : -1;

}

template<typename Num1 , typename Num2 , typename Num3 >
inline bool cinter(Num1 x, Num2 x1, Num3 x2)
{

	return (x < x1) || (x > x2);

}

template<typename Num1 , typename Num2 , typename Num3 >
inline bool inter(Num1 x, Num2 x1, Num3 x2)
{
	return (x1 < x) && (x < x2);
}

template<const int type>// 0 , -1 , 1
double R_f(ARG_3, double d) {
	//err << "Unimplemented template parameter specified for function R_f<" << type << ">";
	exit(1);
}

template<>
inline double R_f<0>(ARG_3, double d) {
	return -t * (1. + 2. * t) * sqr_(1. - std::abs(s)) * d;
}

template<>
inline double R_f<-1>(ARG_3, double d) {
	return sqr_(t+s)-t * (1. + 2. * t) * sqr_(1. - std::abs(s)) * d;
}


template<>
inline double R_f<1>(ARG_3, double d) {
	return sqr_(t - s) - t * (1. + 2. * t) * sqr_(1. - std::abs(s)) * d;
}


template<const int bounds>
static double arg_psi(ARG_3, double d ) {
	double a1 = (1.+2.*t)*sqr_(1.-std::abs(s))*d;
	double a2 = t - double(bounds) * s + std::sqrt(R_f<bounds>(e,t,s,d));
	double result = a1 / a2;
	if ((R_f<bounds>(e, t, s, d) <= 0) || (result <=0. ) || (result >= 1. )) {
		err_ << "Error in function arg_psi, bounds =" << bounds << " e=" << e << " t=" << t << " s=" << s << " d=" << d <<" sqr_(t+s)="<< sqr_(t + s) <<
			" t * (1. + 2. * t)=" << t * (1. + 2. * t) << " result=" << result << " R_f<bounds>(e, t, s, d)="<< R_f<bounds>(e, t, s, d) << '\n';
		throw 1;
	}
	return result;
}

#include <set>
#include <memory>
#include <vector>
template <typename T>
using Matr_T = std::vector<std::vector<T>>;

std::set<size_t> * depth_search(const Matr_T<bool> & matrix_adj, size_t elem, std::set<size_t>* stack = new std::set<size_t>);

#endif