#include "payoff.hpp"
#include "tools.hpp"

#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <cctype>
#include <iostream>

payoff::payoff(std::string name, const std::vector<double>& parameters, const std::function<double(double)>& fct)
		:m_name(name), param(parameters)
	{
		std::string s;

		if (CaseSensitiveIsEqual(m_name, "Call")) {
			payoff_fct = [&](double d1) { return std::max(d1 - param[0], 0.0); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Put")) {
			payoff_fct = [&](double d1) { return std::max(param[0] - d1, 0.0); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Straddle")) {
			payoff_fct = [&](double d1) { return std::max(param[0] - d1, d1 - param[0]); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Call Spread")) {
			payoff_fct = [&](double d1) { return std::max(0.0, std::min(d1 - param[0], param[1] - param[0])); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Put Spread")) {
			payoff_fct = [&](double d1) { return std::min(param[1] - param[0], std::max(param[1] - d1, 0.0)); };
		}
		else if (CaseSensitiveIsEqual(m_name, "Strangle")) {
			payoff_fct = [&](double d1) { return std::max(param[0] - d1, 0.0) + std::max(d1 - param[1], 0.0); };
		}
		else {
			payoff_fct = fct;
		}

	}

std::vector<double> payoff::getparameters() const
{
	return param;
}

std::string payoff::getname() const
{
	return m_name;
}

std::function<double(double)> payoff::getpayoff() const
{
	return payoff_fct;
}

double payoff::getpayoff(double x) const
{
	return payoff_fct(x);
}


