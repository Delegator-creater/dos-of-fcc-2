#ifndef BOUNDS_H
#define BOUNDS_H
#include "Macros.h"
#include <functional>
#include "Interval.h"
class Bounds {
public:

	double Z(double e, double t, double s, double f);

	double s1(ARG_3, double sign);

	double s2(ARG_3, double sign, double sigm);

	double u1(ARG_3, double sign); 

	double u2(ARG_3, double sign, double sigm); 




	char type_bounds;
	void edit_fild(int x) {
		type_bounds = x;
	}


};

using fd_2x = std::function<double(double, double)>;
using fd_3x = std::function<double(double, double, double)>;

class UpBounds final : public Bounds {
public:
	UpBounds() : Bounds() {}
	//*******************************//

	/// Z(e_1 , t , s ,sigma_ * 8s +4t)
	double bounds1(ARG_3, double sigm_);
	//*******************************//

	///Z(e_1 , t , s , -4s^2/t)
	double bounds2(ARG_3);

	//*******************************//
	fd_3x get_bounds(const double t_); 

	Interval get_s_bounds(const double t, const double e, const bool is_s_more_t);



};

class DownBounds final : public Bounds {
public:
	DownBounds() : Bounds() {}

	//*******************************//

	/// Z(e_1 , t , s ,sigma * 8s +4t)
	double bounds1(ARG_3, double sigm);

	//*******************************//

	///Z(e_1 , t , s , -4s^2/t)
	double bounds2(ARG_3);

	//*******************************//

	fd_3x get_bounds(const double t_);

	Interval get_s_bounds(const double t, const double e, const bool is_s_more_t);
};



#endif