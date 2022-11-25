#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
#include <cstddef>

extern const double boltzmann_constant;

struct double3 {
  double x;
  double y;
  double z;
};

void splitString(const std::string& data, const std::string& delim, std::vector<std::string>& dest);

#endif // COMMON_H
