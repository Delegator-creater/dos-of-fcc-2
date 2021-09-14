#ifndef BOUNDS_H
#define BOUNDS_H
#include "Macros.h"
#include <functional>
class Bounds {
protected:

	double Z(double e, double t, double s, double f);

	double s1(ARG_3, double sign);

	double s2(ARG_3, double sign, double sigm);

	double u1(ARG_3, double sign); 

	double u2(ARG_3, double sign, double sigm); 

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

};



#endif