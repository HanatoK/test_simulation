#ifndef POTENTIAL_COSINE_H
#define POTENTIAL_COSINE_H

#include "Potential.h"
#include <cmath>

class Potential_cosine: public Potential {
public:
  Potential_cosine(): Potential() {}
  double getPotential(double3 pos) const override {
    return 1.0 + 0.5 * std::cos(10.0 * pos.x);
  }
  double3 getGradients(double3 pos) const override {
    double3 grad{0, 0, 0};
    grad.x = -1.0 * 0.5 * std::sin(10.0 * pos.x) * 10.0;
    return grad;
  }
};

#endif // POTENTIAL_COSINE_H
