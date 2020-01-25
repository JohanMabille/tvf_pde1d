#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <vector>
#include <string>
#include <functional>

class payoff
{
public:
	explicit payoff(std::string name = "", const std::vector<double>& parameters = { 0 }, const std::function<double(double)>& fct = [](double d1) { return d1; });
	std::string getname() const;
	std::function<double(double)> getpayoff() const;
	double getpayoff(double x) const;
	
	std::vector<double> getparameters() const;

private:
	std::string m_name;
	std::vector<double> param;
	
	std::function<double(double)> payoff_fct;
};


#endif
