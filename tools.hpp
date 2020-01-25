#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <vector>
#include <string>

double sqr(double x);
void size_check(const std::vector<std::vector<double>>& sigma, const std::vector<double>& r, int nt, int nx);
bool CaseSensitiveIsEqual(std::string str1, std::string str2);
#endif