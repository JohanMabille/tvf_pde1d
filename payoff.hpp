#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <vector>
#include <string>
#include <functional>

class payoff
{
public:
	explicit payoff(std::string name = "", const std::vector<double>& parameters = { 0 }, const std::function<double(double)>& fct = [](double d1) { return d1; });
        // Design: all these methods should return const reference
	std::string getname() const;
	std::function<double(double)> getpayoff() const;
	double getpayoff(double x) const;
	
	std::vector<double> getparameters() const;

private:
        // Design: name is convenient for User, but not for code. Consider transforming
        // the name into an enum value to avoid costly conversions each time you
        // need to check the type of payoff.
	std::string m_name;
	std::vector<double> param;
	
	std::function<double(double)> payoff_fct;
};


#endif
