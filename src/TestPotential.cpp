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

double rmsd(const std::vector<double>& __restrict pos_x,
            const std::vector<double>& __restrict pos_y,
            const std::vector<double>& __restrict pos_z,
            const std::vector<double>& __restrict ref_x,
            const std::vector<double>& __restrict ref_y,
            const std::vector<double>& __restrict ref_z) {
  double sum = 0.0;
  const size_t num_atoms = pos_x.size();
  for (size_t i = 0; i < num_atoms; ++i) {
    sum += (pos_x[i] - ref_x[i]) * (pos_x[i] - ref_x[i]);
    sum += (pos_y[i] - ref_y[i]) * (pos_y[i] - ref_y[i]);
    sum += (pos_z[i] - ref_z[i]) * (pos_z[i] - ref_z[i]);
  }
  if (num_atoms > 0) sum /= num_atoms;
  return std::sqrt(sum);
}

void test_grad(Potential& potential, const std::vector<double>& p) {
  std::vector<double> analytical_diff_x(p.size());
  std::vector<double> analytical_diff_y(p.size());
  std::vector<double> analytical_diff_z(p.size());
  std::vector<double> numerical_diff_x(p.size());
  std::vector<double> numerical_diff_y(p.size());
  std::vector<double> numerical_diff_z(p.size());
  potential.getGradients(
    std::vector<double>{p[0]}, std::vector<double>{p[1]}, std::vector<double>{p[2]},
    analytical_diff_x, analytical_diff_y, analytical_diff_z);
  double eps = 0.01;
  for (int i = 0; i < 5; ++i) {
    potential.getNumericalGradient(
      std::vector<double>{p[0]}, std::vector<double>{p[1]}, std::vector<double>{p[2]},
      numerical_diff_x, numerical_diff_y, numerical_diff_z, eps);
    std::cout << "Epsilon = " << eps << ", error = "
              << rmsd(analytical_diff_x, analytical_diff_y, analytical_diff_z,
                      numerical_diff_x, numerical_diff_y, numerical_diff_z) << std::endl;
    eps /= 10.0;
  }
}

void dump_potential(
  Potential& potential, const char* filename, const axis_2d& ax) {
  std::ofstream ofs(filename);
  std::vector<double> pos_x(1);
  std::vector<double> pos_y(1);
  std::vector<double> pos_z(1);
  const size_t nx = std::nearbyint((ax.x_upper - ax.x_lower) / ax.x_width);
  const size_t ny = std::nearbyint((ax.y_upper - ax.y_lower) / ax.y_width);
  ofs << "# 2\n";
  ofs << "# " << ax.x_lower << " " << ax.x_width << " " << nx << " 0\n";
  ofs << "# " << ax.y_lower << " " << ax.y_width << " " << nx << " 0\n";
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      pos_x[0] = ax.x_lower + (i + 0.5) * ax.x_width;
      pos_y[0] = ax.y_lower + (j + 0.5) * ax.y_width;
      pos_z[0] = 0.0;
      const double p = potential.getPotential(pos_x, pos_y, pos_z);
      ofs << pos_x[0] << " " << pos_y[0] << " " << p << std::endl;
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
