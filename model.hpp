#ifndef MODEL_HPP
#define MODEL_HPP
#include "payoff.hpp"

// Design: this could get more flexibility. The model should be an abstract class
// with virtual methods. The current model becomes the default implementation,
// you can consider other models (constant vol and constant rate, mix of constant / variable
// rates and vol, etc ...).
class model
{
public:

	model(double sigma, double r, int n_t, int n_x);
	model(const std::vector<double>& sigma, double r, int n_t, int n_x);
	model(double sigma, const std::vector<double>& r, int n_t, int n_x);
	model(const std::vector<double>& sigma, const std::vector<double>& r, int n_t, int n_x);
	model(const std::vector<std::vector<double>>& sigma, double r, int n_t, int n_x);
	model(const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, int n_t, int n_x);
	
	void get_vol_col(std::vector<double>& mat, int i) const;
	double get_r(int i) const;
	double get_r_avg() const;
	std::vector<double> get_r() const;
	std::vector<std::vector<double>> getSigma() const;
	
private:
	
	std::vector<double> m_r;
	std::vector<std::vector<double>> m_sigma;
};

#endif
