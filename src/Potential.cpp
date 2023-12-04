#include "Potential.h"
#include <cmath>

void Potential::getForces(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z,
  std::vector<double>& f_x,
  std::vector<double>& f_y,
  std::vector<double>& f_z) const {
  getGradients(pos_x, pos_y, pos_z, f_x, f_y, f_z);
  for (size_t i = 0; i < f_x.size(); ++i) {
    f_x[i] *= -1.0;
    f_y[i] *= -1.0;
    f_z[i] *= -1.0;
  }
}

void Potential::getNumericalGradient(
  std::vector<double> pos_x,
  std::vector<double> pos_y,
  std::vector<double> pos_z,
  std::vector<double>& __restrict pos_x_grad,
  std::vector<double>& __restrict pos_y_grad,
  std::vector<double>& __restrict pos_z_grad, double epsilon) const {
  for (size_t i = 0; i < pos_x.size(); ++i) {
    {
      const double saved_xi = pos_x[i];
      pos_x[i] = saved_xi + epsilon;
      const double V_next = getPotential(pos_x, pos_y, pos_z);
      pos_x[i] = saved_xi - epsilon;
      const double V_prev = getPotential(pos_x, pos_y, pos_z);
      pos_x_grad[i] = (V_next - V_prev) / (2.0 * epsilon);
    }
    {
      const double saved_yi = pos_y[i];
      pos_y[i] = saved_yi + epsilon;
      const double V_next = getPotential(pos_x, pos_y, pos_z);
      pos_y[i] = saved_yi - epsilon;
      const double V_prev = getPotential(pos_x, pos_y, pos_z);
      pos_y_grad[i] = (V_next - V_prev) / (2.0 * epsilon);
    }
    {
      const double saved_zi = pos_z[i];
      pos_z[i] = saved_zi + epsilon;
      const double V_next = getPotential(pos_x, pos_y, pos_z);
      pos_z[i] = saved_zi - epsilon;
      const double V_prev = getPotential(pos_x, pos_y, pos_z);
      pos_z_grad[i] = (V_next - V_prev) / (2.0 * epsilon);
    }
  }
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

double BSPotential::getPotential(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z) const {
  const double ux = Ux(pos_x[0]);
  return ux + (m_big_omega_square * (pos_x[0] - pos_y[0]) * (pos_x[0] - pos_y[0]) / 2.0) / m_beta;
}

void BSPotential::getGradients(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z,
  std::vector<double>& __restrict pos_x_grad,
  std::vector<double>& __restrict pos_y_grad,
  std::vector<double>& __restrict pos_z_grad) const {
  // double3 grad;
  pos_x_grad[0] = dUx_dx(pos_x[0]) + (m_big_omega_square * (pos_x[0] - pos_y[0])) / m_beta;
  pos_y_grad[0] = -(m_big_omega_square * (pos_x[0] - pos_y[0])) / m_beta;
  pos_z_grad[0] = 0;
  // return grad;
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

double MBPotential::getPotential(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z
) const {
  double p = 0;
  for (size_t i = 0; i < 4; ++i) {
    p += subterm_i(pos_x[0], pos_y[0], i);
  }
  return m_factor * p;
}

void MBPotential::getGradients(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z,
  std::vector<double>& __restrict pos_x_grad,
  std::vector<double>& __restrict pos_y_grad,
  std::vector<double>& __restrict pos_z_grad
) const {
  double dx, dy;
  for (size_t i = 0; i < 4; ++i) {
    subterm_dxdy_i(pos_x[0], pos_y[0], i, dx, dy);
    pos_x_grad[0] += dx;
    pos_y_grad[0] += dy;
  }
  pos_x_grad[0] *= m_factor;
  pos_y_grad[0] *= m_factor;
}

double TripleWellAlpha::getPotential(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z
) const {
  const auto& x = pos_x[0];
  const auto& y = pos_y[0];
  const double tmp1 = 3.0 * std::exp(-x * x) * (std::exp(-(y - 1.0 / 3.0) * (y - 1.0 / 3.0) / m_alpha) - std::exp(-(y - 5.0 / 3.0) * (y - 5.0 / 3.0) / m_alpha));
  const double tmp2 = -5.0 * std::exp(-y * y / m_alpha) * (std::exp(-(x - 1.0) * (x - 1.0)) + std::exp(-(x + 1.0) * (x + 1.0)));
  const double tmp3 = 0.2 * x * x * x * x;
  const double tmp4 = 0.2 * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) / (m_alpha * m_alpha);
  return tmp1 + tmp2 + tmp3 + tmp4;
}

void TripleWellAlpha::getGradients(
  const std::vector<double>& __restrict pos_x,
  const std::vector<double>& __restrict pos_y,
  const std::vector<double>& __restrict pos_z,
  std::vector<double>& __restrict pos_x_grad,
  std::vector<double>& __restrict pos_y_grad,
  std::vector<double>& __restrict pos_z_grad
) const {
  const auto& x = pos_x[0];
  const auto& y = pos_y[0];
  const double dtmp1_dx = 3.0 * std::exp(-x * x) * (std::exp(-(y - 1.0 / 3.0) * (y - 1.0 / 3.0) / m_alpha) - std::exp(-(y - 5.0 / 3.0) * (y - 5.0 / 3.0) / m_alpha)) * (-2.0 * x);
  const double dtmp1_dy = 3.0 * std::exp(-x * x) * (std::exp(-(y - 1.0 / 3.0) * (y - 1.0 / 3.0) / m_alpha) * (-2.0 * (y - 1.0 / 3.0) / m_alpha) - std::exp(-(y - 5.0 / 3.0) * (y - 5.0 / 3.0) / m_alpha) * (-2.0 * (y - 5.0 / 3.0) / m_alpha));
  const double dtmp2_dx = -5.0 * std::exp(-y * y / m_alpha) * (std::exp(-(x - 1.0) * (x - 1.0)) * (-2.0 * (x - 1.0)) + std::exp(-(x + 1.0) * (x + 1.0)) * (-2.0 * (x + 1.0)));
  const double dtmp2_dy = -5.0 * std::exp(-y * y / m_alpha) * (-2.0 * y / m_alpha) * (std::exp(-(x - 1.0) * (x - 1.0)) + std::exp(-(x + 1.0) * (x + 1.0)));
  const double dtmp3_dx = 0.2 * 4 * x * x * x;
  const double dtmp3_dy = 0.0;
  const double dtmp4_dx = 0.0;
  const double dtmp4_dy = 0.2 * 4 * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) * (y - 1.0 / 3.0) / (m_alpha * m_alpha);
  pos_x_grad[0] = dtmp1_dx + dtmp2_dx + dtmp3_dx + dtmp4_dx;
  pos_y_grad[0] = dtmp1_dy + dtmp2_dy + dtmp3_dy + dtmp4_dy;
  pos_z_grad[0] = 0;
}
