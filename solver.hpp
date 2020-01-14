#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"

class solver_edp
{
public:
	
	solver_edp(model pde_model);
	std::vector<double> solve_pde();
	
	std::vector<double> product_inverse(std::vector<std::vector<double>> trig_mat, std::vector<double> d);
	std::vector<double> trig_matmul(std::vector<std::vector<double>> trig_mat, std::vector<double> x);
private:
	model s_pde_model;
	
};

#endif