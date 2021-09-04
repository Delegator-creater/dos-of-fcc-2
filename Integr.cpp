#include "Integr.h"

const size_t start_points = 100;
double integration(std::function<double(double , void*)> f , void * const prm, double a, double  b, const double error1, const double error2) {
	a += error1;
	b -= error1;
	size_t n = start_points;
	
	auto const sumiring = [&]() {
		double h = (b - a) / n;
		double  sum = 0;
		for (double count = a; count <= b; count += h)
			sum += (f(count, prm) + f(count + h, prm)) / 2;
		return sum * h;
	};

	double sum1 = sumiring();
	n *= 2;
	double sum2 = sumiring();
	while (abs(sum1 - sum2) >= error2) {
		n *= 2;
		sum1 = sum2;
		sum2 = sumiring();
	}

	return sum2;
}