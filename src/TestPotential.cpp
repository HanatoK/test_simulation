#include "Potential.h"
#include <fstream>
#include <fmt/format.h>
#include <iostream>

struct axis_2d {
  double x_lower;
  double x_upper;
  double x_width;
  double y_lower;
  double y_upper;
  double y_width;
};

double rmsd(double3 x, double3 y) {
  double sum = 0.0;
  sum += (x.x - y.x) * (x.x - y.x);
  sum += (x.y - y.y) * (x.y - y.y);
  sum += (x.z - y.z) * (x.z - y.z);
  return std::sqrt(sum);
}

void test_grad(Potential& potential, const std::vector<double>& p) {
  const auto analytical_diff = potential.getGradients(double3{p[0], p[1], p[2]});
  double eps = 0.01;
  for (int i = 0; i < 5; ++i) {
    const auto numerical_diff = potential.getNumericalGradient(double3{p[0], p[1], p[2]}, eps);
    std::cout << "Epsilon = " << eps << ", error = " << rmsd(analytical_diff, numerical_diff) << std::endl;
    eps /= 10.0;
  }
}

void dump_potential(
  Potential& potential, const char* filename, const axis_2d& ax) {
  std::ofstream ofs(filename);
  double3 point{0, 0, 0};
  const size_t nx = std::nearbyint((ax.x_upper - ax.x_lower) / ax.x_width);
  const size_t ny = std::nearbyint((ax.y_upper - ax.y_lower) / ax.y_width);
  ofs << "# 2\n";
  ofs << "# " << ax.x_lower << " " << ax.x_width << " " << nx << " 0\n";
  ofs << "# " << ax.y_lower << " " << ax.y_width << " " << nx << " 0\n";
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      point.x = ax.x_lower + (i + 0.5) * ax.x_width;
      point.y = ax.y_lower + (j + 0.5) * ax.y_width;
      const double p = potential.getPotential(point);
      ofs << point.x << " " << point.y << " " << p << std::endl;
    }
  }
}

int main() {
  {
    // BS potential
    BSPotential potential(2.0, 2.2, 1.0 / (300.0 * 0.0019872041));
    std::vector<double> p{0.0, -0.05, 0};
    axis_2d ax{-6.0, 6.0, 0.02, -6.0, 6.0, 0.02};
    test_grad(potential, p);
    p[0] = -0.6;
    test_grad(potential, p);
    p[1] = 1.1;
    test_grad(potential, p);
    dump_potential(potential, "BSPotential.dat", ax);
  }
  {
    // MB potential
    MBPotential potential;
    std::vector<double> p{0.0, -0.05, 0};
    axis_2d ax{-2.0, 1.0, 0.01, -1.0, 2.0, 0.01};
    test_grad(potential, p);
    p[0] = -0.6;
    test_grad(potential, p);
    p[1] = 0.75;
    test_grad(potential, p);
    dump_potential(potential, "MBPotential.dat", ax);
  }
  return 0;
}