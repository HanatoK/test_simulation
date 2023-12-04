#ifndef POTENTIAL_COSINE_H
#define POTENTIAL_COSINE_H

#include "Potential.h"
#include <cmath>

class Potential_cosine: public Potential {
public:
  Potential_cosine(): Potential() {}
  double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z
  ) const override {
    return 1.0 + 0.5 * std::cos(10.0 * pos_x[0]);
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
    pos_x_grad[0] = -1.0 * 0.5 * std::sin(10.0 * pos_x[0]) * 10.0;
    // return grad;
  }
};

#endif // POTENTIAL_COSINE_H
