#include "spt_func.h"
#include <deque>
#include <string>

std::deque<std::string> split(std::string str, const char * delimiter) {
	std::deque<std::string> result;
	size_t index = 0;
	for (size_t i = 0; i < str.size(); ++i)
		if (str[i] == *delimiter) {
			result.push_back(str.substr(index, i - index));
			index = i + 1;
		}
	return result;
}

std::set<size_t> *  depth_search(const Matr_T<bool> & matrix_adj, size_t elem, std::set<size_t>* stack)
{
	stack->insert(elem);
	for (size_t i = 0; i < matrix_adj[elem].size(); ++i) {//i in range(len(matrix_adj[elem - 1])) //: # [0, 1, 2, 3]
		if (matrix_adj[elem][i])
			if (stack->find(i) == stack->end())
				depth_search(matrix_adj, i, stack);
	}
	return stack;
}

double sqr_(double x) {
	return x*x;
}

