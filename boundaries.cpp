#include "boundaries.hpp"
#include "tools.hpp"

#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

bound::bound(const payoff& f, const mesh& grille, const std::vector<double>& conditions, std::string method)
	: b_f(f), b_mesh(grille), b_method(method)
{
	if ((!(CaseSensitiveIsEqual(b_method,"Dirichlet"))) && (!(CaseSensitiveIsEqual(b_method,"Neumann"))))
	{
		std::cout<< "please enter Dirichelet or Neumann only" << std::endl;
	}
	b_conditions_down = conditions[0];
	b_conditions_up = conditions[1];
}

bound::bound(const payoff& f, const mesh& grille, std::string method)
	: b_f(f), b_mesh(grille), b_method(method)
{
	if ((!(CaseSensitiveIsEqual(b_method,"Dirichlet"))) && (!(CaseSensitiveIsEqual(b_method,"Neumann"))))
	{
		std::cout<< "please enter Dirichlet or Neumann only" << std::endl;
	}
	
	if (CaseSensitiveIsEqual(b_method,"Dirichlet"))
	{
		
	}
	else if (CaseSensitiveIsEqual(b_method,"Neumann"))
	{
		double h = 0.00001;
		b_conditions_down = std::exp(b_mesh.get_Smin())*(b_f.getpayoff(std::exp(b_mesh.get_Smin()) + h) - b_f.getpayoff(std::exp(b_mesh.get_Smin())))/h;
		b_conditions_up = std::exp(b_mesh.get_Smax())*(b_f.getpayoff(std::exp(b_mesh.get_Smax()) + h) - b_f.getpayoff(std::exp(b_mesh.get_Smax())))/h;;
	}
}

void bound::adapt_mat(std::vector<std::vector<double>>& mat, std::vector<std::vector<double>>& mat_inv, double theta, double r, const std::vector<double>& sigma) const
{
	double dx = b_mesh.get_dx();
	double dt = b_mesh.get_dt();
	
        // Design: This is WAY TOO COMPLICATED!
        // Compute the tridiag system (LHS matrix and RHS vector) first, and then apply the
        // boundary conditions on them.
	if (CaseSensitiveIsEqual(b_method,"Neumann"))
	{
	
		mat[1][0] = -1/dx;
		mat[2][0] = 1/dx;
		mat[0][b_mesh.get_nx()-1] = -1/dx;
		mat[1][b_mesh.get_nx()-1] = 1/dx;
			
		mat_inv[1][1] = 1+dt*theta*(sqr(sigma[1]/dx)*0.5 + r - (1./2.*sigma[1]*sigma[1] - r)/(2*dx));
		mat_inv[0][1] = -dt*theta/(2*dx)*(-sqr(sigma[1])/dx - sqr(sigma[1])/2. + r)*dx;
		
		mat_inv[1][b_mesh.get_nx()-2] = 1+dt*theta*(sqr(sigma[b_mesh.get_nx()-2]/dx)*0.5 + r + (1./2.*sigma[b_mesh.get_nx()-2]*sigma[b_mesh.get_nx()-2] - r)/(2*dx));
		mat_inv[2][b_mesh.get_nx()-2] = dx * dt*theta/(2*dx)*(-sqr(sigma[b_mesh.get_nx()-2])/dx + sqr(sigma[b_mesh.get_nx()-2])/2. - r);
	}
	else if (CaseSensitiveIsEqual(b_method,"Dirichlet"))
	{
		
		mat_inv[1][b_mesh.get_nx()-1] = std::exp(r*dt);
		mat_inv[1][0] = std::exp(r*dt);
	}
	
}

void bound::get_boundaries(std::vector<double>& sol, double T, double dt, int i, double r) const
{
	if (dt*i == T)
	{
		for(int j=1; j<sol.size()-1; ++j)
		{
			sol[j] = b_f.getpayoff(std::exp(b_mesh.get_Smin() + j*b_mesh.get_dx()));
		}
		
		if (CaseSensitiveIsEqual(b_method,"Dirichlet"))
		{
			sol[0] = b_f.getpayoff(std::exp(b_mesh.get_Smin())*std::exp(r*T));
			sol[sol.size()-1] = b_f.getpayoff(std::exp(b_mesh.get_Smax())*std::exp(r*T));
		}
		else if (CaseSensitiveIsEqual(b_method,"Neumann"))
		{
			sol[0] = b_f.getpayoff(std::exp(b_mesh.get_Smin()));
			sol[sol.size()-1] = b_f.getpayoff(std::exp(b_mesh.get_Smax()));
		}
	}
	else
	{	
		if (CaseSensitiveIsEqual(b_method,"Dirichlet"))
		{	
			// Nothing to transform when using Dirichlet
		}
		else if (CaseSensitiveIsEqual(b_method,"Neumann"))
		{
			sol[sol.size()-1] = sol[b_mesh.get_nx()-1] * b_mesh.get_dx() + sol[b_mesh.get_nx()-2];
			sol[0] = -sol[0] * b_mesh.get_dx() + sol[1];
		}
	}
}
