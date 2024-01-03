#ifndef POTENTIAL_LJ_H
#define POTENTIAL_LJ_H

#include "Potential.h"
#include <cmath>

// a simple LJ potential
class LJPotential: public Potential {
private:
  std::vector<std::vector<double>> m_epsilon;
  std::vector<std::vector<double>> m_sigma;
  std::vector<std::vector<double>> m_sigma_12;
  std::vector<std::vector<double>> m_sigma_6;
public:
  LJPotential(std::vector<std::vector<double>> epsilon, std::vector<std::vector<double>> sigma):m_epsilon(epsilon), m_sigma(sigma), m_sigma_12(m_sigma), m_sigma_6(m_sigma) {
    for (size_t i = 0; i < m_sigma_12.size(); ++i) {
      for (size_t j = 0; j < m_sigma_12.size(); ++j) {
        m_sigma_6[i][j] = m_sigma[i][j] * m_sigma[i][j] * m_sigma[i][j];
        m_sigma_6[i][j] *= m_sigma_6[i][j];
        m_sigma_12[i][j] = m_sigma_6[i][j] * m_sigma_6[i][j];
      }
    }
  }
  double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z) const override {
    double potential = 0;
    for (size_t i = 0; i < pos_x.size(); ++i) {
      for (size_t j = i + 1; j < pos_x.size(); ++j) {
        double dist2_ij = 0;
        dist2_ij += (pos_x[i] - pos_x[j]) * (pos_x[i] - pos_x[j]);
        dist2_ij += (pos_y[i] - pos_y[j]) * (pos_y[i] - pos_y[j]);
        dist2_ij += (pos_z[i] - pos_z[j]) * (pos_z[i] - pos_z[j]);
        const double dist6_ij = dist2_ij * dist2_ij * dist2_ij;
        const double dist12_ij = dist6_ij * dist6_ij;
        const double V_ij = 4.0 * m_epsilon[i][j] * (m_sigma_12[i][j] / dist12_ij - m_sigma_6[i][j] / dist6_ij);
        potential += V_ij;
      }
    }
    return potential;
  }
  void getGradients(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad) const override {
    for (size_t i = 0; i < pos_x.size(); ++i) {
      for (size_t j = i + 1; j < pos_x.size(); ++j) {
        double dist2_ij = 0;
        dist2_ij += (pos_x[i] - pos_x[j]) * (pos_x[i] - pos_x[j]);
        dist2_ij += (pos_y[i] - pos_y[j]) * (pos_y[i] - pos_y[j]);
        dist2_ij += (pos_z[i] - pos_z[j]) * (pos_z[i] - pos_z[j]);
        const double dist_ij = std::sqrt(dist2_ij);
        const double dist6_ij = dist2_ij * dist2_ij * dist2_ij;
        const double dist7_ij = dist6_ij * dist_ij;
        const double dist13_ij = dist7_ij * dist6_ij;
        const double factor = 4.0 * m_epsilon[i][j] * (12.0 * m_sigma_12[i][j] / dist13_ij - 6.0 * m_sigma_6[i][j] / dist7_ij);
        // TODO
      }
    }
  }
};

#endif // POTENTIAL_LJ_H
