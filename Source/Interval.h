#ifndef _INTERVAL_H_
#define _INTERVAL_H_
#include <vector>
#include <deque>
#include <memory>


class Interval {

	class One_Interval {
		double bounds[2];
		friend Interval;
	public:

		One_Interval() {
			bounds[0] = 1.;
			bounds[1] = -1.;
		}

		One_Interval(const One_Interval & o_i) {
			for (size_t i = 0; i < 2; ++i)
				bounds[i] = o_i.bounds[i];
			
		}

		One_Interval(One_Interval&& o_i) {
			for (size_t i = 0; i < 2; ++i)
				bounds[i] = o_i.bounds[i];
		}

		One_Interval(double left_bounds, double right_bounds) {
			bounds[0] = left_bounds;
			bounds[1] = right_bounds;
		}
		

		const double* get_bounds() const { return bounds; }

		bool in(const double) const;

		bool empty() const;

		Interval operator+(const One_Interval) const;

		One_Interval operator*(const One_Interval &) const;

		Interval operator/(const One_Interval&)const;

		const One_Interval & operator=( const One_Interval & o_i) {

			for (size_t i = 0; i < 2; ++i)
				bounds[i] = o_i.bounds[i];

			return *this;
		}


	};

	std::vector<One_Interval> bounds;

	friend One_Interval;

	std::unique_ptr<std::vector<One_Interval>> get_one_ints()const;

	void add_O_I(One_Interval&&);

	void add_O_I(const One_Interval&);

public:

	Interval(double left_bounds , double right_bounds) {
		bounds.push_back(One_Interval(left_bounds , right_bounds));
	}

	Interval(Interval&& interval) {
		for (const auto & i : interval.bounds)
			bounds.push_back(i);
		interval.bounds.clear();
	}

	Interval(const Interval& interval) {
		
		bounds.clear();
		for (const auto & i : interval.bounds) 
			bounds.push_back(i);
		
	}

	Interval(const One_Interval& o_i) {
		bounds.push_back(o_i);
	}

	Interval(){
		bounds.push_back(One_Interval(1,-1));
	}


	Interval operator+(const Interval& other_inter);
	
	//Interval operator+(Interval&& other_inter);

	bool empty() const;

	Interval abs() const ;

	Interval operator*( const Interval & other_inter) const;

	Interval operator/(const Interval& other_inter) const;

	Interval operator*(Interval&& ) const;

	Interval & operator=(const Interval&);

	Interval & operator=(Interval&&);

	const std::vector<One_Interval>& get_bounds() const { return bounds; }

	using O_I = One_Interval;
};

#endif