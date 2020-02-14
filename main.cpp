#include "closed_form.hpp"
#include "payoff.hpp"
#include "model.hpp"
#include "boundaries.hpp"
#include "solver.hpp"
#include "mesh.hpp"
#include "tools.hpp"

#include <iostream>
#include <algorithm>
#include <chrono>
#include <string>

int main(int argc, char* argv[])
{	

	//std::function<double(double)> f = [&](double d1) { return std::max(d1 - 100, 0.0); };

	payoff pp = payoff("Call", {100});
	
	double mat = 1;
	double S = 100;
	double theta = 1./2;
	
	int nx = 1000;
	int nt = 252;
	
	char delim = ';'; // Change delim to modify the csv separator of sigma.csv and rate.csv
	
	std::vector<double> r;
	std::vector<std::vector<double>> sigma;
	
	std::string user_input;
	
        // Design: sigma should not be hardcoded but asked to the user (at K, Mat)
        // if the user wants to input a matrix.
	mesh grille(S, mat, nx, nt, 0.2);
	
        // This is a good start for giving the user the ability to input her own
        // data. However, this still lacks flexibility. A more complete approach
        // would be to provide a DataReader base class, and different implementations:
        // one that reads data from a file, one that asks on the console, on that could
        // read from a DataBase (the idea is not to provide all possible implementations,
        // but to provide an interface easy to implement for who wants to add a data source).
	while (user_input != "y" && user_input != "n")
	{
		std::cout <<  "Do you want to input your own sigma matrix ? y/n ";
		std::cin >> user_input;
	}
	
	if (user_input == "y")
	{
		std::string tmp;
		grille.export_empty_sigma();
		std::cout << "Please input values in the sigma.csv just outputed and press ENTER";
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		
		grille.read_sigma(sigma, delim);
		std::cout << "Sigma read successfully" << std::endl;
	}
	else
	{
		std::cout << "Using sigma coded in the main.cpp" << std::endl;
		sigma.resize(nx, std::vector<double>(nt));
		for (int i= 0; i < nx; ++i)
		{
			for (int j=0; j < nt; ++j)
			{
				sigma[i][j] = 0.2;
			}
		}
	}
	
	user_input = "";
	
	while (user_input != "y" && user_input != "n")
	{
		std::cout <<  "Do you want to input your own rate vector ? y/n ";
		std::cin >> user_input;
	}
	
	if (user_input == "y")
	{
		std::string tmp;
		grille.export_empty_rate();
		std::cout << "Please input values in the rate.csv just outputed and press ENTER";
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		
		grille.read_rate(r, delim);
		std::cout << "Rate read successfully" << std::endl;
	}
	else
	{
		std::cout << "Using rate coded in the main.cpp" << std::endl;
		r.resize(nt);
		for (int i= 0; i < nt; ++i)
		{
			r[i]= 0.04;
		}
	}
	
	try
	{
	auto start = std::chrono::steady_clock::now();
	model model_pde(sigma, r, nt, nx);
		
	bound boundaries(pp, grille, "Neumann");
	solver_edp solver_model(model_pde, grille, boundaries, pp, theta);
	solver_model.solve_pde(1);

	auto end = std::chrono::steady_clock::now();
	
	std::cout << "Time taken solving : " << std::chrono::duration <double, std::milli> (end - start).count() << " ms" << std::endl;
	
	//solver_model.export_csv(); // Export the results to a csv file called output.csv
	solver_model.print_results(); // Uncomment this line to display the results in the terminal
	}
	catch(std::exception& e)
	{
		std::cout << e.what() << std::endl;
		std::cout << "Exiting programme without solving." << std::endl;
	}
	
	
	return 0;

}
