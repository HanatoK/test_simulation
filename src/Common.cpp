#include "Common.h"
#include <cmath>

const double boltzmann_constant = 0.0019872041;

void splitString(const std::string& data, const std::string& delim, std::vector<std::string>& dest) {
  size_t index = 0, new_index = 0;
  std::string tmpstr;
  while (index != data.length()) {
    new_index = data.find(delim, index);
    if (new_index != std::string::npos) tmpstr = data.substr(index, new_index - index);
    else tmpstr = data.substr(index, data.length());
    if (!tmpstr.empty()) {
      dest.push_back(tmpstr);
    }
    if (new_index == std::string::npos) break;
    index = new_index + 1;
  }
}

double beta(double temperature) {
  return 1.0 / (temperature * boltzmann_constant);
}

// double3& double3::operator+=(const double3& rhs) {
//   this->x += rhs.x;
//   this->y += rhs.y;
//   this->z += rhs.z;
//   return *this;
// }
//
// double3& double3::operator-=(const double3& rhs) {
//   this->x -= rhs.x;
//   this->y -= rhs.y;
//   this->z -= rhs.z;
//   return *this;
// }
//
// double3& double3::operator*=(const double3& rhs) {
//   this->x *= rhs.x;
//   this->y *= rhs.y;
//   this->z *= rhs.z;
//   return *this;
// }
//
// double3& double3::operator*=(const double& rhs) {
//   this->x *= rhs;
//   this->y *= rhs;
//   this->z *= rhs;
//   return *this;
// }
//
// double3& double3::operator/=(const double& rhs) {
//   this->x /= rhs;
//   this->y /= rhs;
//   this->z /= rhs;
//   return *this;
// }
//
// double3 operator+(double3 lhs, const double3& rhs) {
//   lhs += rhs;
//   return lhs;
// }
//
// double3 operator-(double3 lhs, const double3& rhs) {
//   lhs -= rhs;
//   return lhs;
// }
//
// double3 operator*(double3 lhs, const double3& rhs) {
//   lhs *= rhs;
//   return lhs;
// }
//
// double3 operator*(const double& lhs, double3 rhs) {
//   rhs *= lhs;
//   return rhs;
// }
//
// double3 operator*(double3 lhs, const double& rhs) {
//   lhs *= rhs;
//   return lhs;
// }
//
// double3 operator/(double3 lhs, const double& rhs) {
//   lhs /= rhs;
//   return lhs;
// }
//
// double3 double3::exp(const double3 lhs) {
//   return double3{std::exp(lhs.x), std::exp(lhs.y), std::exp(lhs.z)};
// }
//
// double3 double3::sqrt(const double3 lhs) {
//   return double3{std::sqrt(lhs.x), std::sqrt(lhs.y), std::sqrt(lhs.z)};
// }
//
// double3 operator/(const double& lhs, double3 rhs) {
//   return double3{lhs / rhs.x, lhs / rhs.y, lhs / rhs.z};
// }
