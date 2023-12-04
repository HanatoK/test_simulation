#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
#include <cstddef>

extern const double boltzmann_constant;

/* unit conversion:
  delta_t: ps
  gamma: g/mol*(ps^-1)
  potential: kcal/mol
  force: kcal/(mol*angstrom)
  r: angstrom

  force * delta_t / gamma:
  =kcal/(mol*angstrom) * ps / (g/mol*(ps^-1))
  =kcal/(angstrom) * ps*ps / g
  =4.184 * kJ * ps * ps / (g * angstrom)
  =4.184 * 1e3 *1e3 g * m^2 * s^-2 * ps^2 / (g*angstrom)
  =4.184*1e6 * (1e10)^2 * angstrom^2 * s^-2 * ps^2 / (angstrom)
  =4.184*1e26 * (1e12)^-2 ps^-2 * ps^2 * angstrom
  =4.184*1e26 * 1e-24 angstrom
  =418.4 angstrom
*/

inline static constexpr double conversion_factor = 418.4;


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
