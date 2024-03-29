#include "Potential.h"
#include <cmath>

double3 Potential::getForces(double3 pos) const {
  double3 grad = getGradients(pos);
  grad.x *= -1.0;
  grad.y *= -1.0;
  grad.z *= -1.0;
  return grad;
}

double3 Potential::getNumericalGradient(double3 pos, double epsilon) const {
  const double dV_dx = (getPotential(double3{pos.x + epsilon, pos.y, pos.z}) - getPotential(double3{pos.x - epsilon, pos.y, pos.z})) / (2.0 * epsilon);
  const double dV_dy = (getPotential(double3{pos.x, pos.y + epsilon, pos.z}) - getPotential(double3{pos.x, pos.y - epsilon, pos.z})) / (2.0 * epsilon);
  const double dV_dz = (getPotential(double3{pos.x, pos.y, pos.z + epsilon}) - getPotential(double3{pos.x, pos.y, pos.z - epsilon})) / (2.0 * epsilon);
  return double3{dV_dx, dV_dy, dV_dz};
}

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

double MBPotential::subterm_i(double x, double y, size_t i) const {
  const double diff_x = x - param_x_0[i];
  const double diff_y = y - param_y_0[i];
  const double tmp = param_a[i] * diff_x * diff_x + param_b[i] * diff_x * diff_y + param_c[i] * diff_y * diff_y;
  return param_A[i] * std::exp(tmp);
}

void MBPotential::subterm_dxdy_i(double x, double y, size_t i, double& dx, double& dy) const {
  const double diff_x = x - param_x_0[i];
  const double diff_y = y - param_y_0[i];
  const double tmp = param_a[i] * diff_x * diff_x + param_b[i] * diff_x * diff_y + param_c[i] * diff_y * diff_y;
  const double dtmp_dx = 2.0 * param_a[i] * diff_x + param_b[i] * diff_y;
  const double dtmp_dy = param_b[i] * diff_x + 2.0 * param_c[i] * diff_y;
  dx = param_A[i] * std::exp(tmp) * dtmp_dx;
  dy = param_A[i] * std::exp(tmp) * dtmp_dy;
}

double MBPotential::getPotential(double3 pos) const {
  double p = 0;
  for (size_t i = 0; i < 4; ++i) {
    p += subterm_i(pos.x, pos.y, i);
  }
  return m_factor * p;
}

double3 MBPotential::getGradients(double3 pos) const {
  double3 grad{0, 0, 0};
  double dx, dy;
  for (size_t i = 0; i < 4; ++i) {
    subterm_dxdy_i(pos.x, pos.y, i, dx, dy);
    grad.x += dx;
    grad.y += dy;
  }
  grad.x *= m_factor;
  grad.y *= m_factor;
  return grad;
}

double TripleWellAlpha::getPotential(double3 pos) const {
  const auto& x = pos.x;
  const auto& y = pos.y;
  const double tmp1 = 3.0 * std::exp(-x * x) * (std::exp(-(y - 1.0 / 3.0) * (y - 1.0 / 3.0) / m_alpha) - std::exp(-(y - 5.0 / 3.0) * (y - 5.0 / 3.0) / m_alpha));
  const double tmp2 = -5.0 * std::exp(-y * y / m_alpha) * (std::exp(-(x - 1.0) * (x - 1.0)) + std::exp(-(x + 1.0) * (x + 1.0)));
  const double tmp3 = 0.2 * x * x * x * x;
  const double tmp4 = 0.2 * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) / (m_alpha * m_alpha);
  return tmp1 + tmp2 + tmp3 + tmp4;
}

double3 TripleWellAlpha::getGradients(double3 pos) const {
  double3 grad;
  const auto& x = pos.x;
  const auto& y = pos.y;
  const double dtmp1_dx = 3.0 * std::exp(-x * x) * (std::exp(-(y - 1.0 / 3.0) * (y - 1.0 / 3.0) / m_alpha) - std::exp(-(y - 5.0 / 3.0) * (y - 5.0 / 3.0) / m_alpha)) * (-2.0 * x);
  const double dtmp1_dy = 3.0 * std::exp(-x * x) * (std::exp(-(y - 1.0 / 3.0) * (y - 1.0 / 3.0) / m_alpha) * (-2.0 * (y - 1.0 / 3.0) / m_alpha) - std::exp(-(y - 5.0 / 3.0) * (y - 5.0 / 3.0) / m_alpha) * (-2.0 * (y - 5.0 / 3.0) / m_alpha));
  const double dtmp2_dx = -5.0 * std::exp(-y * y / m_alpha) * (std::exp(-(x - 1.0) * (x - 1.0)) * (-2.0 * (x - 1.0)) + std::exp(-(x + 1.0) * (x + 1.0)) * (-2.0 * (x + 1.0)));
  const double dtmp2_dy = -5.0 * std::exp(-y * y / m_alpha) * (-2.0 * y / m_alpha) * (std::exp(-(x - 1.0) * (x - 1.0)) + std::exp(-(x + 1.0) * (x + 1.0)));
  const double dtmp3_dx = 0.2 * 4 * x * x * x;
  const double dtmp3_dy = 0.0;
  const double dtmp4_dx = 0.0;
  const double dtmp4_dy = 0.2 * 4 * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) / (m_alpha * m_alpha);
  grad.x = dtmp1_dx + dtmp2_dx + dtmp3_dx + dtmp4_dx;
  grad.y = dtmp1_dy + dtmp2_dy + dtmp3_dy + dtmp4_dy;
  grad.z = 0;
  return grad;
}
