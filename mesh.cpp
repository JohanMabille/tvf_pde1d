#include "tools.hpp"
#include "mesh.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>


mesh::mesh(double S, double T, int& n_x, int n_t, double sigma)
	: m_nt(n_t), m_S(S), m_sigma(sigma)
{
	if (n_x % 2 == 0)
	{
		n_x += 1;
	}
	
	m_nx = n_x;
	
	m_dt = T/(m_nt - 1);
	
	m_Smax = std::log(S) + 5 * sigma * std::sqrt(T);
	m_Smin = std::log(S) - 5 * sigma * std::sqrt(T);
	
	m_dx = (m_Smax - m_Smin)/(m_nx - 1);
}

void mesh::print_mesh() const
{
	for(int j=0; j < m_nx; ++j)
	{
		std::cout << std::exp(m_Smin + j*m_dx) << ",";
	}
	
	std::cout << std::endl;
	
	for(int i=0; i < m_nt; ++i)
	{
		std::cout << i*m_dt << std::endl;
	}
}

void mesh::export_empty_sigma(std::string f_name) const
{
	std::ofstream f(f_name);
	
	f << ",";
	
	for (int i = 0; i < m_nt - 1; ++i)
	{
		f << i * m_dt << ",";
	}
	
	f << (m_nt - 1) * m_dt << "\n";
	
	for (int i = 0; i < m_nx;++i)
	{
		f << std::exp(m_Smin + i * m_dx)  << "\n";
	}
	
	f.close();
}

void mesh::export_empty_rate(std::string f_name) const
{
	std::ofstream f(f_name);
	
	for (int i = 0; i < m_nt - 1; ++i)
	{
		f << i * m_dt << ",";
	}
	
	f << (m_nt - 1) * m_dt << "\n";
	
	f.close();
}

void mesh::read_sigma(std::vector<std::vector<double>>& M, char sep) const
{
	
	// Got it from the internet and then adapted for the code
	
	std::ifstream in("sigma.csv");
	std::string line;
	int cpt_r = 0;

	while (std::getline(in, line))                        // read a whole line of the file
	{
		if (cpt_r > 0)
		{
			std::stringstream ss(line);                   // put it in a stringstream (internal stream)
			std::vector<double> row;
			std::string data;
			while (std::getline(ss, data, sep))           // read (string) items up to a comma
			{
				row.push_back(std::stod(data));         // use stod() to convert to double; put in row vector
			}
			if ( row.size() > 0)
			{
				row.erase(row.begin(), row.begin()+1);      // Delete first col which is the x
				M.push_back(row); 				        // add non-empty rows to matrix
			}	
		}

		cpt_r += 1;
	}
}

void mesh::read_rate(std::vector<double>& r, char sep) const
{
	std::ifstream in("rate.csv");
	std::string line;
	int cpt_r = 0;

	while (std::getline(in, line))                        // read a whole line of the file
	{
		if (cpt_r > 0)
		{
			std::stringstream ss(line);                   // put it in a stringstream (internal stream)
			std::vector<double> row;
			std::string data;
			while (std::getline(ss, data, sep))           // read (string) items up to a comma
			{
				r.push_back(std::stod(data));         // use stod() to convert to double; put in row vector
			}

		}

		cpt_r += 1;
	}
}

double mesh::get_dx() const
{
	return m_dx;
}

double mesh::get_dt() const
{
	return m_dt;
}

double mesh::get_Smax() const
{
	return m_Smax;
}

double mesh::get_Smin() const
{
	return m_Smin;
}

int mesh::get_nx() const
{
	return m_nx;
}

int mesh::get_nt() const
{
	return m_nt;
}

double mesh::get_S() const
{
	return m_S;
}

double mesh::get_sigma() const
{
	return m_sigma;
}

double mesh::get_mat() const
{
	return (m_nt-1)*m_dt;
}
