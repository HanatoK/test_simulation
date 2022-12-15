#include "Potential.h"
#include <fstream>
#include <fmt/format.h>
#include <iostream>

double rmsd(double3 x, double3 y) {
  double sum = 0.0;
  sum += (x.x - y.x) * (x.x - y.x);
  sum += (x.y - y.y) * (x.y - y.y);
  sum += (x.z - y.z) * (x.z - y.z);
  return std::sqrt(sum);
}

void test_grad(BSPotential& potential, const std::vector<double>& p) {
  const auto analytical_diff = potential.getGradients(double3{p[0], p[1], p[2]});
  const auto numerical_diff = potential.getNumericalGradient(double3{p[0], p[1], p[2]});
  std::cout << "Error = " << rmsd(analytical_diff, numerical_diff) << std::endl;
}

void dump_potential(BSPotential& potential) {
  std::ofstream ofs("potential.dat");
  double3 point{0, 0, 0};
  double x_lower = -6.0;
  double x_upper = 6.0;
  double x_width = 0.02;
  size_t nx = std::nearbyint((x_upper - x_lower) / x_width);
  double y_lower = -6.0;
  double y_upper = 6.0;
  double y_width = 0.02;
  size_t ny = std::nearbyint((y_upper - y_lower) / y_width);
  ofs << "# 2\n";
  ofs << "# " << x_lower << " " << x_width << " " << nx << " 0\n";
  ofs << "# " << y_lower << " " << y_width << " " << nx << " 0\n";
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      point.x = x_lower + (i + 0.5) * x_width;
      point.y = y_lower + (j + 0.5) * y_width;
      const double p = potential.getPotential(point);
      ofs << point.x << " " << point.y << " " << p << std::endl;
    }
  }
}

int main() {
  BSPotential potential(2.0, 2.2, 1.0 / (300.0 * 0.0019872041));
  std::vector<double> p{0.0, -0.05, 0};
  test_grad(potential, p);
  p[0] = -0.6;
  test_grad(potential, p);
  p[1] = 1.1;
  test_grad(potential, p);
  dump_potential(potential);
  return 0;
}