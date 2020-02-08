#ifndef SOLVER_HPP
#define SOLVER_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"
#include "boundaries.hpp"

class solver_edp
{
public:
	
	solver_edp(const model& pde_model, const mesh& grille, const bound& boundary, const payoff& f, double theta);
	
	void solve_pde(bool vega_bool = 0);
	void export_csv(std::string f_name = "output.csv") const;
	void print_results() const;
	
        // Design: Why public members instead of private ones
        // with getters (and not setters) ?
        // Where is the theta (greek)?
	std::vector<double> solution;
	std::vector<double> delta;
	std::vector<double> gamma;
	std::vector<double> vega;

private:
	const model& s_pde_model;
	const mesh& s_mesh;
	const bound& s_bound;
	const payoff& s_f;
	
	double s_theta;
	
	void pde_matrix(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, const std::vector<double>& sigma, const std::vector<double>& sigma_plus, double r, double r_plus, double theta, double dt, double dx, int nx, int i) const;
	void product_inverse(std::vector<double>& x, std::vector<std::vector<double>>& trig_mat, std::vector<double>& d);
	void trig_matmul(std::vector<double>& res, std::vector<std::vector<double>>& trig_mat, std::vector<double>& x);
	
        // Implementation: unused member
	std::vector<std::vector<double>> s_cdt;
	
};

#endif
