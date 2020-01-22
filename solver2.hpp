#ifndef SOLVER2_HPP
#define SOLVER2_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"
#include "boundaries2.hpp"

class solver_edp2
{
public:
	
	solver_edp2(model pde_model, mesh grille, bound2 boundary, payoff f, double theta);
	
	void solve_pde(const bool& vega_bool = 0);
		
	std::vector<double> solution;
	std::vector<double> delta;
	std::vector<double> gamma;
	std::vector<double> vega;

private:
	model s_pde_model;
	mesh s_mesh;
	bound2 s_bound;
	payoff s_f;
	
	double s_theta;
	
	void pde_matrix(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, const std::vector<double>& sigma, const std::vector<double>& sigma_plus, double r, double theta, double dt, double dx, int nx, int i);
	void product_inverse(std::vector<double>& x, std::vector<std::vector<double>> trig_mat, std::vector<double> d);
	void trig_matmul(std::vector<double>& res, std::vector<std::vector<double>> trig_mat, std::vector<double> x);
	
	std::vector<std::vector<double>> s_cdt;
	

};

#endif