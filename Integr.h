#ifndef MY_INTEGRATION
#define MY_INTEGRATION
#include <functional>

double integration(std::function<double(double , void *)>, void*, double, double, double, double);

#endif

