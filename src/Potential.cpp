#include "Potential.h"
#include <cmath>

double BSPotential::Ux(double x) const {
  if (x <= -0.5 * m_x0) {
    return (-m_Delta + m_omega * m_omega * (x + m_x0) * (x + m_x0) / 2.0) / m_beta;
  } else if ((x > -0.5 * m_x0) && (x < 0.5 * m_x0)) {
    return (-m_omega * m_omega * x * x / 2.0) / m_beta;
  } else {
    return (-m_Delta + m_omega * m_omega * (x - m_x0) * (x - m_x0) / 2.0) / m_beta;
  }
}

double BSPotential::dUx_dx(double x) const {
  if (x <= -0.5 * m_x0) {
    return (m_omega * m_omega * (x + m_x0)) / m_beta;
  } else if ((x > -0.5 * m_x0) && (x < 0.5 * m_x0)) {
    return (-m_omega * m_omega * x) / m_beta;
  } else {
    return (m_omega * m_omega * (x - m_x0)) / m_beta;
  }
}

double BSPotential::getPotential(double3 pos) const {
  const double ux = Ux(pos.x);
  return ux + (m_big_omega_square * (pos.x - pos.y) * (pos.x - pos.y) / 2.0) / m_beta;
}

double3 BSPotential::getGradients(double3 pos) const {
  double3 grad;
  grad.x = dUx_dx(pos.x) + (m_big_omega_square * (pos.x - pos.y)) / m_beta;
  grad.y = -(m_big_omega_square * (pos.x - pos.y)) / m_beta;
  grad.z = 0;
  return grad;
}

double3 BSPotential::getForces(double3 pos) const {
  double3 grad = getGradients(pos);
  grad.x *= -1.0;
  grad.y *= -1.0;
  grad.z *= -1.0;
  return grad;
}
