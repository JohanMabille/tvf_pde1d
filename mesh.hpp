#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <string>

class mesh
{
public:

	mesh(double S, double T, int& n_x, int n_t, double sigma);
	
	void print_mesh() const;
        // Design: the mesh should not be responsible for this.
	void export_empty_sigma(std::string f_name = "sigma.csv") const;
	void read_sigma(std::vector<std::vector<double>>& M, char sep) const;
	void export_empty_rate(std::string f_name = "rate.csv") const;
	void read_rate(std::vector<double>& r, char sep) const;
	
        // Design: the mesh should be responsible for returning
        // x(i). This would avoid computing exp(s_min+(i)*s_mesh.get_dx())
        // everywhere. Same for t(i).
	double get_dt() const;
	double get_dx() const;
        // Design: consistency between S (spot) and Smax and Smin (log(spot))
	double get_Smax() const;
	double get_Smin() const;
	double get_S() const;
	double get_mat() const;
	double get_sigma() const;
	
        // Design: consider using size_t instead of int since
        // these numbrs cannot be negative
	int get_nx() const;
	int get_nt() const;
	
private:
	
	int m_nx;
	int m_nt;
	
	double m_Smax;
	double m_Smin;
	double m_dt;
	double m_dx;
	double m_S;
	double m_sigma;
};
	

#endif
