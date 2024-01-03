#include "Potential.h"
#include <cmath>

class Potential_flatwells: public Potential {
public:
  Potential_flatwells(): Potential() {}
  double getPotential(double3 pos) const override {
    return 5.0 * std::exp (-20.0 * pos.x * pos.x);
  }
  double3 getGradients(double3 pos) const override {
    double3 grad{0, 0, 0};
    grad.x = 5.0 * std::exp (-20.0 * pos.x * pos.x) * (-20.0) * 2.0 * pos.x;
    return grad;
  }
};
