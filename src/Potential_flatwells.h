#include "Potential.h"
#include <cmath>

class Potential_flatwells: public Potential {
public:
  Potential_flatwells(): Potential() {}
  double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z
  ) const override {
    return 5.0 * std::exp (-20.0 * pos_x[0] * pos_x[0]);
  }
  void getGradients(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad
  ) const override {
    // double3 grad{0, 0, 0};
    pos_x_grad[0] = 5.0 * std::exp (-20.0 * pos_x[0] * pos_x[0]) * (-20.0) * 2.0 * pos_x[0];
    // return grad;
  }
};
