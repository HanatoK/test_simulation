#include "Common.h"

class BSPotential {
private:
  double m_omega;
  double m_x0;
  double m_beta;
  double m_big_omega_square;
  double m_Delta;
  double Ux(double x) const;
  double dUx_dx(double x) const;
public:
  BSPotential(double omega, double x0, double beta):
    m_omega(omega), m_x0(x0), m_beta(beta),
    m_big_omega_square(1.01 * omega * omega),
    m_Delta(omega * omega * x0 * x0 / 4.0) {}
  double getPotential(double3 pos) const;
  double3 getForces(double3 pos) const;
  double3 getGradients(double3 pos) const;
};
