#ifndef SPT_FUNC_H
#define SPT_FUNC_H

#include <functional>
#include <string>
#include <deque>
#include<cmath>

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

#endif
