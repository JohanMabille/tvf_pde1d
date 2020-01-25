#include "tools.hpp"
#include <exception>
#include <iostream>
#include <exception>

double sqr(double x)
{
	return x*x;
}

void size_check(const std::vector<std::vector<double>>& m_sigma, const std::vector<double>& m_r, int nt, int nx)
{
	std::string msg = "";
	
	if (m_sigma.size() != nx)
	{
		msg = msg + "Sigma size mispecified in x. ";
	}
	if (m_sigma[0].size() != nt)
	{
		msg = msg + "Sigma size mispecified in t. ";
	}
	if (m_r.size() != nt)
	{
		msg = msg + "Rate size mispecified in t.";
	}
	
	if (msg != "")
	{
		throw std::runtime_error(msg);
	}
}

bool CaseSensitiveIsEqual(std::string str1, std::string str2)
{
	return ((str1.size() == str2.size()) && std::equal(str1.begin(), str1.end(), str2.begin(), [](char c1, char c2) {
		return (c1 == c2 || std::toupper(c1) == std::toupper(c2));
		}));
}