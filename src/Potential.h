#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <vector>

class Potential {
public:
  Potential() = default;
  virtual ~Potential() = default;
  virtual double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z) const = 0;
  virtual void getGradients(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad) const = 0;
  virtual void getForces(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict f_x,
    std::vector<double>& __restrict f_y,
    std::vector<double>& __restrict f_z) const final;
  virtual void getNumericalGradient(
    std::vector<double> pos_x,
    std::vector<double> pos_y,
    std::vector<double> pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad, double epsilon = 0.00001) const final;
};

class BSPotential: public Potential {
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
    Potential(), m_omega(omega), m_x0(x0), m_beta(beta),
    m_big_omega_square(1.01 * omega * omega),
    m_Delta(omega * omega * x0 * x0 / 4.0) {}
  double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z
  ) const override;
  void getGradients(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad
  ) const override;
};

class MBPotential: public Potential {
private:
  double m_factor;
  static constexpr double param_A[] = {-200, -100, -170, 15};
  static constexpr double param_a[] = {-1, -1, -6.5, 0.7};
  static constexpr double param_b[] = {0, 0, 11, 0.6};
  static constexpr double param_c[] = {-10, -10, -6.5, 0.7};
  static constexpr double param_x_0[] = {1, 0, -0.5, -1};
  static constexpr double param_y_0[] = {0, 0.5, 1.5, 1};
  double subterm_i(double x, double y, size_t i) const;
  void subterm_dxdy_i(double x, double y, size_t i, double& dx, double& dy) const;
public:
  MBPotential(double factor = 1.0 / 20.0): Potential(), m_factor(factor) {}
  double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z
  ) const override;
  void getGradients(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad
  ) const override;
};

class TripleWellAlpha: public Potential {
private:
  double m_alpha;
public:
  TripleWellAlpha(double alpha = 1.0): m_alpha(alpha) {}
  double getPotential(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z
  ) const override;
  void getGradients(
    const std::vector<double>& __restrict pos_x,
    const std::vector<double>& __restrict pos_y,
    const std::vector<double>& __restrict pos_z,
    std::vector<double>& __restrict pos_x_grad,
    std::vector<double>& __restrict pos_y_grad,
    std::vector<double>& __restrict pos_z_grad
  ) const override;
};

#endif // POTENTIAL_H
