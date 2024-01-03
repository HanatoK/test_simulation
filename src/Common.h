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
  double3& operator+=(const double3& rhs);
  double3& operator-=(const double3& rhs);
  double3& operator*=(const double3& rhs);
  double3& operator*=(const double& rhs);
  double3& operator/=(const double& rhs);
  static double3 exp(const double3 lhs);
  static double3 sqrt(const double3 lhs);
};

double3 operator+(double3 lhs, const double3& rhs);
double3 operator-(double3 lhs, const double3& rhs);
double3 operator*(double3 lhs, const double3& rhs);
double3 operator*(const double& lhs, double3 rhs);
double3 operator*(double3 lhs, const double& rhs);
double3 operator/(double3 lhs, const double& rhs);
double3 operator/(const double& lhs, double3 rhs);


void splitString(const std::string& data, const std::string& delim, std::vector<std::string>& dest);

template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp = 2) {
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::abs(x - y) <
             std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
         // unless the result is subnormal
         || std::abs(x - y) < std::numeric_limits<T>::min();
}

double beta(double temperature);

#endif // COMMON_H
