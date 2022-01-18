#include "Interval.h"
#include <initializer_list>
#include <memory>


bool Interval::One_Interval::in(const double pun) const
{
	return (bounds[0] < pun ) && (pun < bounds[1]);
}

Interval::One_Interval Interval::One_Interval::operator*(const Interval::One_Interval & other_interval) const {
	
	if (this->empty() || other_interval.empty())
		return One_Interval(1,-1);
	if ( (bounds[0] >= other_interval.bounds[1]) || (bounds[1] <= other_interval.bounds[0]))
		return One_Interval(1, -1);

	return Interval::One_Interval(std::max(bounds[0] , other_interval.bounds[0]), std::min(bounds[1] , other_interval.bounds[1]));

}

Interval Interval::One_Interval::operator/(const One_Interval& other_interval)const {
	Interval result;
	result.add_O_I(*this);
	if (!other_interval.empty()) {

		bool cond_matrix[2][2];
		for (size_t i = 0; i < 2; ++i)
			for (size_t j = 0; j < 2; ++j)
				cond_matrix[i][j] = bounds[i] > other_interval.bounds[j];

		const auto sum_matrix = [&]() {
			return cond_matrix[0][0] + cond_matrix[0][1] + cond_matrix[1][0] + cond_matrix[1][1];
		};

		if ((sum_matrix() == 0) || (sum_matrix() == 4)) {
			return result;
		}

		else {
			Interval other_result;
			Interval::One_Interval first_interval(this->bounds[0] , other_interval.bounds[0]);
			Interval::One_Interval second_interval( other_interval.bounds[1] , this->bounds[1]);
			other_result.add_O_I(first_interval);
			other_result.add_O_I(second_interval);
			return other_result;
		}
	}
	else
		return result;
}

bool Interval::One_Interval::empty() const
{
	return bounds[0] >= bounds[1];
}

void Interval::add_O_I(Interval::One_Interval&& o_i)
{
	if (!o_i.empty()) {
		if (bounds.size() == 1)
			if (bounds[0].empty())
				bounds.clear();
		bounds.push_back(o_i);
	}
}

void Interval::add_O_I(const One_Interval& o_i)
{
	if (!o_i.empty()) {

		if (bounds.size() == 1)
			if (bounds[0].empty())
				bounds.clear();
		bounds.push_back(o_i);
	}
}

Interval Interval::One_Interval::operator+(const Interval::One_Interval other_interval) const {
	if (this->empty() && other_interval.empty())
		return Interval(One_Interval(1, -1));
	if (this->empty())
		return Interval(other_interval);
	if (other_interval.empty())
		return Interval(*this);
	if ((bounds[0] >= other_interval.bounds[1]) || (bounds[1] <= other_interval.bounds[0])) {
		Interval ret_inter(*this);
		ret_inter.add_O_I(other_interval);
		return ret_inter;
	}

	return Interval::One_Interval( std::min(bounds[0], other_interval.bounds[0]), std::max(bounds[1], other_interval.bounds[1]));

}



Interval Interval::operator*(const Interval& other_inter) const {

	

	auto v_O_I_this = &bounds;
	auto v_O_I_other = &other_inter.bounds;

	//size_t len_intersects = bounds.size() + other_inter.bounds.size();
	//Matr_T<bool> intersects;
	//intersects.resu

	Interval new_interval;
	for (size_t i = 0; i < bounds.size(); ++i)
		for (size_t j = 0; j < other_inter.bounds.size(); ++j)
			if (!((*v_O_I_this)[i] * (*v_O_I_other)[j]).empty()) 
				new_interval.add_O_I(  (*v_O_I_this)[i] * (*v_O_I_other)[j]  );


	return new_interval;

}
Interval Interval::operator/(const Interval& other_inter) const {
	if ((bounds.size() == 1) && (other_inter.bounds.size() == 1))
		return Interval::One_Interval(bounds[0].bounds[0], bounds[0].bounds[1]) / Interval::One_Interval(other_inter.bounds[0].bounds[0], other_inter.bounds[0].bounds[1]);
	Interval result;
	for (auto & i : bounds ) {
		Interval edit_i; 
		edit_i.add_O_I(i);
		for (auto & j : other_inter.bounds) {
			edit_i = edit_i / j;
		}
		for (auto & o_i : edit_i.bounds)
			result.add_O_I(o_i);
	}
	return result;
}

Interval Interval::operator*(Interval&& other_inter) const {
	Interval new_interval(other_inter);
	return (*this) * new_interval;
}

Interval & Interval::operator=(const Interval& interval) {
	if (&interval != this) {
		for (const auto & i : interval.bounds) {
			bounds.push_back(i);
		}
	}
	return (*this);

}

Interval & Interval::operator=(Interval&& interval) {
	if (&interval != this) {

		bounds.clear();
		for (const auto & i : interval.bounds)
			bounds.push_back(i);
	}
	return (*this);
	
}


#include <map>
#include "spt_func.h"
Interval Interval::operator+(const Interval & other_inter) {
	
	size_t len_intersect = bounds.size() + other_inter.bounds.size();
	std::vector<std::vector<bool>> intersects;
	intersects.resize(len_intersect);
	for (auto& i : intersects)
		i.resize(len_intersect);

	Interval new_interval;
	auto v_O_I_this  = &bounds;
	auto v_O_I_other = &other_inter.get_bounds();

	for (size_t i = 0; i < bounds.size(); ++i)
		for (size_t j = 0; j < other_inter.bounds.size(); ++j){
			bool is_empty = !((*v_O_I_this)[i] * (*v_O_I_other)[j]).empty();
			if (is_empty) {
				intersects[i][other_inter.bounds.size() + j] = is_empty;
				intersects[other_inter.bounds.size() + j][i] = is_empty;
			}
		}

	std::map<size_t, Interval::One_Interval> * map_o_i = new std::map<size_t , One_Interval>;
	
	size_t number = 0;
	for (const auto & i : bounds)
		map_o_i->insert(std::pair<size_t , Interval::One_Interval>(number++, One_Interval(i.bounds[0] , i.bounds[1]) ));
	for (const auto& i : other_inter.bounds)
		map_o_i->insert(std::pair<size_t, Interval::One_Interval>(number++, One_Interval(i.bounds[0], i.bounds[1])  ));

	while(!map_o_i->empty()) {
		auto elem = *map_o_i->begin();
		auto set_o_i = depth_search(intersects, elem.first);
		Interval::One_Interval first_elem(elem.second);

		for (auto i : *set_o_i) {
			first_elem = (first_elem + map_o_i->find(i)->second).bounds[0];
			//map_o_i.emplace(i);
		}
		new_interval.add_O_I(first_elem);
		auto* new_map_o_i = new std::map<size_t, One_Interval>;
		for (const auto& i : *map_o_i)
			if (set_o_i->find(i.first) == set_o_i->end())
				new_map_o_i->insert(i);
		std::swap(map_o_i, new_map_o_i);
		
		delete new_map_o_i;
		delete set_o_i;
	}
	return new_interval;
}

bool Interval::empty() const {
	if (bounds.size() == 1)
		return bounds[0].bounds[0] > bounds[0].bounds[1];
	else
		return false;
}
Interval Interval::abs() const {
	Interval new_interval;
	for (auto& o_i : bounds) {
		new_interval.add_O_I(One_Interval(-o_i.bounds[1] , -o_i.bounds[0]));

	}
	return  new_interval+(*this);
}