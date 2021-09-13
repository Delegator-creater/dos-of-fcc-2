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

using fd_2x = std::function<double(double, double)>;
using fd_3x = std::function<double(double, double, double)>;

class UpBounds final : public Bounds {
public:
	UpBounds() : Bounds() {}
	//*******************************//

	/// Z(e_1 , t , s ,sigma_ * 8s +4t)
	double bounds1(ARG_3, double sigm_);/* {
		double sigm = double(SIGN(s)) * sigm_;
		auto return_value1 = std::function<double(void)>([=]() {return double(1); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, sigm_ * 8. * s + 4. * t); });
		if (t > -0.5 && t < 0.)
			std::swap(return_value1, return_value2);
		if (inter(std::abs(s), u2(e, t, s, -1., sigm), u2(e, t, s, 1., sigm)))
			return return_value1();
		else
			return return_value2();

	}*/
	//*******************************//

	///Z(e_1 , t , s , -4s^2/t)
	double bounds2(ARG_3); /* {
		auto return_value1 = std::function<double(void)>([=]() {return double(1); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, 4. * sqr_(s) / t); });
		if ((0. < t) && (t < 0.5))
			std::swap(return_value1, return_value2);
		if (inter(std::abs(s), u1(e, t, s, -1.), u1(e, t, s, 1.)))
			return return_value1();
		else
			return return_value2();
	}*/

	//*******************************//
	fd_3x get_bounds(const double t_); /* {
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
	}*/




};

class DownBounds final : public Bounds {
public:
	DownBounds() : Bounds() {}

	//*******************************//

	/// Z(e_1 , t , s ,sigma * 8s +4t)
	double bounds1(ARG_3, double sigm);/* {
		auto return_value1 = std::function<double(void)>([=]() {return double(0); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, sigm * 8. * s + 4. * t); });
		if (std::abs(t) > 0.5)
			std::swap(return_value1, return_value2);
		if ((s2(e, t, s, -1., sigm) < s) && (s < s2(e, t, s, 1., sigm)))
			return return_value2();
		else
			return return_value1();
	}*/

	//*******************************//

	///Z(e_1 , t , s , -4s^2/t)
	double bounds2(ARG_3);/* {
		auto return_value1 = std::function<double(void)>([=]() {return double(0); });
		auto return_value2 = std::function<double(void)>([=]() {return Z(e, t, s, 4. * sqr_(s) / t); });
		if ((0. < t) && (t < 1.))
			std::swap(return_value1, return_value2);
		if (inter(s, s1(e, t, s, -1.), s1(e, t, s, 1.)))
			return return_value2();
		else
			return return_value1();
	}*/

	//*******************************//




	fd_3x get_bounds(const double t_);/* {
		if (t_ >= 0.) {
			return [=](ARG_3) {return this->bounds1(e, t, s, 1.); };
		}// return Z(... , 8s - 4t)
		else if (t_ >= -0.5) {
			return [=](ARG_3) {return ((std::abs(s) > std::abs(t)) ? this->bounds1(e, t, s, 1.) : this->bounds2(e, t, s)); };
		}// return Z(... , 8s - 4t) or Z(... , 4s^2/t)
		else {
			return [=](ARG_3) {return this->bounds1(e, t, s, -1.); };
		}// return Z(... , -8s - 4t) 
	}*/

};



#endif