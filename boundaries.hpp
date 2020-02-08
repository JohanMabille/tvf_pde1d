#ifndef BOUND_HPP
#define BOUND_HPP
#include <vector>
#include "model.hpp"
#include "mesh.hpp"

// Design: if you want to add new boundary conditions (gamma null on boundaries, Robin conditions),
// you need to change the implementation, which is not a good thing.
//
// A better approach would be to have a bound base class with pure virtual methods:
// class bound
// {
// public:
//
//     virtual ~bound();
//
//     virtual void adapt_mat(...) const = 0;
//     virtual void get_boundaries(...) = const = 0;
//
// protected:
//     
//     // constructors
//     bound(....);
//
//     // members
//     .....
// };
//
// class dirichlet_bound : public bound
// {
// public:
//
//     // Constructors
//     dirichlet_bound(...);
//
//     virtual ~dirichlet_bound();
//
//     void adapt_mat(...) const override;
//     void get_boundaries(...) const override;
// };
//
// Then provide factory functions to instantiate the right boundary condition,
// and work with the base class only:
//
// enum class bound_type
// {
//     dirichlet,
//     neuman,
//     robin
// };
//
// bound* make_bound(bound_type b);
//
// You can also split the adapt_mat / get_boundaries functions
// to update lower and upper parts. This way you can choose different conditions
// on each boundary.

class bound
{
public:

	bound(const payoff& f, const mesh& grille, const std::vector<double>& conditions, std::string method = "Dirichlet");
	bound(const payoff& f, const mesh& grille, std::string method = "Dirichlet");
	
	void adapt_mat(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, double theta, double r, const std::vector<double>& sigma) const;
	void get_boundaries(std::vector<double>& sol, double T, double dt, int i, double r) const;
		
private:
	
	const payoff& b_f;
	const mesh& b_mesh;
	
	std::string b_method;
	double b_conditions_up;
	double b_conditions_down;
};
	

#endif
