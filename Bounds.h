#ifndef BOUNDS_H
#define BOUNDS_H
class Bounds {
protected:

	inline double Z(double e, double t, double s, double f); /*{
		return ((f - e) / 4. + (1. - 2. * t) * sqr_(s))
			/ ((1. + 2. * t) * sqr_(1. - std::abs(s)));
	}*/ 

	inline double s1(ARG_3, double sign);/* {
		double a = e /
			(1. - 2. * t + 1. / t);
		if (a > 0.)
			return sign * 0.50 * std::sqrt(a);
		else
			return 0.;
	}*/

	inline double s2(ARG_3, double sign, double sigm); /*{
		auto a = 1. - (1. - 2. * t) * (t - e / 4.);
		if (a > 0)
			return -sigm / (1. - 2. * t) + sign * std::sqrt(a)
			/ std::abs(1. - 2. * t);
		else
			return 0;
	}*/ 

	inline double u1(ARG_3, double sign); /*{
		auto a = -4 * t * t + 1 / 4 / t * (e + 4) - t * e + 3;
		if (a > 0.)
			return -t * (1. + 2. * t) / (1. - 4. * sqr_(t)) + sign * std::abs(t) / (1. - 4. * t * t) * std::sqrt(a);
		else
			return 0;
	}*/ 

	inline double u2(ARG_3, double sign, double sigm); /*{
		double a = 2. - e * t + 2. * sigm * (1. + 2. * t);
		if (a > 0)
			return (1. + 2. * t + sigm) / (4. * t) + sign * std::sqrt(a)
			/ std::abs(4. * t);
		else
			return 0;
	}*/ 


};
#endif