#include "Simulation.h"
#include "Common.h"

void Simulation::initializeVelocities() {
  for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
    const double alpha = std::sqrt(conversion_factor / (beta() * m_atom_group.m_mass[i]));
    m_atom_group.m_vel_x[i] = randGaussian();
    m_atom_group.m_vel_y[i] = randGaussian();
    m_atom_group.m_vel_z[i] = randGaussian();
    m_atom_group.m_vel_x[i] *= alpha;
    m_atom_group.m_vel_y[i] *= alpha;
    m_atom_group.m_vel_z[i] *= alpha;
  }
}

void Simulation::runLangevinDynamics(
    int64_t steps, const double timestep, const double friction,
    std::function<void(
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> forceFunction,
    std::function<double(
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict)> potentialFunction,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> forceCallback,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> velocityCallback,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> positionCallback,
    std::function<void(double&)> kineticEnergyCallback,
    std::function<void(double&)> potentialEnergyCallback,
    std::function<void(int64_t)> stepCallback,
    std::function<void()> runCallback) {
  // double3 frictions = {friction, friction, friction};
  const std::vector<double> friction_x(m_atom_group.m_num_atoms, friction);
  const std::vector<double> friction_y(m_atom_group.m_num_atoms, friction);
  const std::vector<double> friction_z(m_atom_group.m_num_atoms, friction);
  runLangevinDynamics(
    steps, timestep,
    friction_x, friction_y, friction_z,
    forceFunction,
    potentialFunction,
    forceCallback,
    velocityCallback,
    positionCallback,
    kineticEnergyCallback,
    potentialEnergyCallback,
    stepCallback,
    runCallback);
}

void Simulation::runLangevinDynamics(
  int64_t steps, const double timestep,
  // const double3& frictions,
  const std::vector<double>& __restrict friction_x,
  const std::vector<double>& __restrict friction_y,
  const std::vector<double>& __restrict friction_z,
  std::function<void(
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> forceFunction,
    std::function<double(
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict)> potentialFunction,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> forceCallback,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> velocityCallback,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> positionCallback,
  std::function<void(double&)> kineticEnergyCallback,
  std::function<void(double&)> potentialEnergyCallback,
  std::function<void(int64_t)> stepCallback,
  std::function<void()> runCallback) {
  std::vector<double> force_x(m_atom_group.m_num_atoms, 0);
  std::vector<double> force_y(m_atom_group.m_num_atoms, 0);
  std::vector<double> force_z(m_atom_group.m_num_atoms, 0);
  positionCallback(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
  double Ek = kineticEnergy();
  kineticEnergyCallback(Ek);
  double Ep = potentialFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
  potentialEnergyCallback(Ep);
  stepCallback(m_step);
  forceFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z, force_x, force_y, force_z);
  forceCallback(force_x, force_y, force_z);
  runCallback();
  velocityCallback(m_atom_group.m_vel_x, m_atom_group.m_vel_y, m_atom_group.m_vel_z);
  // BAOAB
  const auto half_timestep = 0.5 * timestep;
  // const auto half_timestep_divide_mass = half_timestep / m_mass * conversion_factor;
  for (m_step = 1; m_step <= steps; ++m_step) {
    for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
      const auto factor1_x = std::exp(-1.0 * friction_x[i] * timestep);
      const auto factor1_y = std::exp(-1.0 * friction_y[i] * timestep);
      const auto factor1_z = std::exp(-1.0 * friction_z[i] * timestep);
      const auto factor2_x = std::sqrt(conversion_factor / (beta() * m_atom_group.m_mass[i])) * std::sqrt(1.0 - std::exp(-2.0 * friction_x[i] * timestep));
      const auto factor2_y = std::sqrt(conversion_factor / (beta() * m_atom_group.m_mass[i])) * std::sqrt(1.0 - std::exp(-2.0 * friction_y[i] * timestep));
      const auto factor2_z = std::sqrt(conversion_factor / (beta() * m_atom_group.m_mass[i])) * std::sqrt(1.0 - std::exp(-2.0 * friction_z[i] * timestep));
      const auto half_timestep_divide_mass = half_timestep / m_atom_group.m_mass[i] * conversion_factor;
      // update v_{i+1/2}
      m_atom_group.m_vel_x[i] += half_timestep_divide_mass * force_x[i];
      m_atom_group.m_vel_y[i] += half_timestep_divide_mass * force_y[i];
      m_atom_group.m_vel_z[i] += half_timestep_divide_mass * force_z[i];
      // update x_{i+1/2}
      m_atom_group.m_pos_x[i] += half_timestep * m_atom_group.m_vel_x[i];
      m_atom_group.m_pos_y[i] += half_timestep * m_atom_group.m_vel_y[i];
      m_atom_group.m_pos_z[i] += half_timestep * m_atom_group.m_vel_z[i];
      // Langevin thermostat, full step
      m_atom_group.m_vel_x[i] = factor1_x * m_atom_group.m_vel_x[i] + factor2_x * randGaussian();
      m_atom_group.m_vel_y[i] = factor1_y * m_atom_group.m_vel_y[i] + factor2_y * randGaussian();
      m_atom_group.m_vel_z[i] = factor1_z * m_atom_group.m_vel_z[i] + factor2_z * randGaussian();
      // Langevin thermostat, full step
      m_atom_group.m_vel_x[i] = factor1_x * m_atom_group.m_vel_x[i] + factor2_x * randGaussian();
      m_atom_group.m_vel_y[i] = factor1_y * m_atom_group.m_vel_y[i] + factor2_y * randGaussian();
      m_atom_group.m_vel_z[i] = factor1_z * m_atom_group.m_vel_z[i] + factor2_z * randGaussian();
      // update x_{i+1}
      m_atom_group.m_pos_x[i] += half_timestep * m_atom_group.m_vel_x[i];
      m_atom_group.m_pos_y[i] += half_timestep * m_atom_group.m_vel_y[i];
      m_atom_group.m_pos_z[i] += half_timestep * m_atom_group.m_vel_z[i];
    }
    // collect CVs and write traj
    positionCallback(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
    Ep = potentialFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
    potentialEnergyCallback(Ep);
    stepCallback(m_step);
    // update f_{i+1}
    forceFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z, force_x, force_y, force_z);
    forceCallback(force_x, force_y, force_z);
    for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
      const auto half_timestep_divide_mass = half_timestep / m_atom_group.m_mass[i] * conversion_factor;
      // update v_{i+1}
      m_atom_group.m_vel_x[i] += half_timestep_divide_mass * force_x[i];
      m_atom_group.m_vel_y[i] += half_timestep_divide_mass * force_y[i];
      m_atom_group.m_vel_z[i] += half_timestep_divide_mass * force_z[i];
    }
    Ek = kineticEnergy();
    kineticEnergyCallback(Ek);
    runCallback();
    velocityCallback(m_atom_group.m_vel_x, m_atom_group.m_vel_y, m_atom_group.m_vel_z);
  }
}

void Simulation::runBrownianDynamics(
  int64_t steps, const double timestep,
  const std::vector<double>& friction_x,
  const std::vector<double>& friction_y,
  const std::vector<double>& friction_z,
  std::function<void(
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> forceFunction,
    std::function<double(
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict,
      const std::vector<double>& __restrict)> potentialFunction,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> forceCallback,
    std::function<void(
      std::vector<double>& __restrict,
      std::vector<double>& __restrict,
      std::vector<double>& __restrict)> positionCallback,
  std::function<void(double&)> potentialEnergyCallback,
  std::function<void(int64_t)> stepCallback,
  std::function<void()> runCallback) {
  std::vector<double> force_x(m_atom_group.m_num_atoms, 0);
  std::vector<double> force_y(m_atom_group.m_num_atoms, 0);
  std::vector<double> force_z(m_atom_group.m_num_atoms, 0);
  positionCallback(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
  double Ep = potentialFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
  potentialEnergyCallback(Ep);
  stepCallback(m_step);
  forceFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z, force_x, force_y, force_z);
  forceCallback(force_x, force_y, force_z);
  runCallback();
  // const auto inv_gamma = timestep / (m_mass * frictions) * conversion_factor;
  // const auto factor = std::sqrt(2.0 / beta()) * double3::sqrt(inv_gamma);
  for (m_step = 1; m_step <= steps; ++m_step) {
    for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
      const auto inv_gamma_x = timestep / (m_atom_group.m_mass[i] * friction_x[i]) * conversion_factor;
      const auto inv_gamma_y = timestep / (m_atom_group.m_mass[i] * friction_y[i]) * conversion_factor;
      const auto inv_gamma_z = timestep / (m_atom_group.m_mass[i] * friction_z[i]) * conversion_factor;
      const auto factor_x = std::sqrt(2.0 / beta()) * std::sqrt(inv_gamma_x);
      const auto factor_y = std::sqrt(2.0 / beta()) * std::sqrt(inv_gamma_y);
      const auto factor_z = std::sqrt(2.0 / beta()) * std::sqrt(inv_gamma_z);
      // update x_{i+1}
      m_atom_group.m_pos_x[i] += inv_gamma_x * force_x[i] + factor_x * randGaussian();
      m_atom_group.m_pos_y[i] += inv_gamma_y * force_y[i] + factor_y * randGaussian();
      m_atom_group.m_pos_z[i] += inv_gamma_z * force_z[i] + factor_z * randGaussian();
    }
    positionCallback(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
    Ep = potentialFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z);
    potentialEnergyCallback(Ep);
    stepCallback(m_step);
    // update f_{i+1}
    forceFunction(m_atom_group.m_pos_x, m_atom_group.m_pos_y, m_atom_group.m_pos_z, force_x, force_y, force_z);
    forceCallback(force_x, force_y, force_z);
    runCallback();
  }
}
