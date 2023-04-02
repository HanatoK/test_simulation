#include "Simulation.h"

void Simulation::initializeVelocities() {
  m_velocities.x = randGaussian();
  m_velocities.y = randGaussian();
  m_velocities.z = randGaussian();
  const double alpha = std::sqrt(conversion_factor / (beta() * m_mass));
  m_velocities.x *= alpha;
  m_velocities.y *= alpha;
  m_velocities.z *= alpha;
}

void Simulation::runLangevinDynamics(
  int64_t steps, const double timestep,
  const double friction,
  std::function<double3(double3)> forceFunction,
  std::function<double(double3)> potentialFunction,
  std::function<void(double3&)> forceCallback,
  std::function<void(double3&)> velocityCallback,
  std::function<void(double3&)> positionCallback,
  std::function<void(double&)> kineticEnergyCallback,
  std::function<void(double&)> potentialEnergyCallback,
  std::function<void(int64_t)> stepCallback,
  std::function<void()> runCallback) {
  double3 frictions = {friction, friction, friction};
  runLangevinDynamics(
    steps, timestep, frictions,
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
  const double3& frictions,
  std::function<double3(double3)> forceFunction,
  std::function<double(double3)> potentialFunction,
  std::function<void(double3&)> forceCallback,
  std::function<void(double3&)> velocityCallback,
  std::function<void(double3&)> positionCallback,
  std::function<void(double&)> kineticEnergyCallback,
  std::function<void(double&)> potentialEnergyCallback,
  std::function<void(int64_t)> stepCallback,
  std::function<void()> runCallback) {
  positionCallback(m_positions);
  double Ek = kineticEnergy();
  kineticEnergyCallback(Ek);
  double Ep = potentialFunction(m_positions);
  potentialEnergyCallback(Ep);
  stepCallback(m_step);
  runCallback();
  double3 force = forceFunction(m_positions);
  forceCallback(force);
  velocityCallback(m_velocities);
  const auto factor1 = double3::exp(-1.0 * frictions * timestep);
  const auto factor2 = std::sqrt(conversion_factor / (beta() * m_mass)) *
                       double3::sqrt(double3{1.0, 1.0, 1.0} - double3::exp(-2.0 * frictions * timestep));
  // BAOAB
  const auto half_timestep = 0.5 * timestep;
  const auto half_timestep_divide_mass = half_timestep / m_mass * conversion_factor;
  for (m_step = 1; m_step <= steps; ++m_step) {
    // update v_{i+1/2}
    m_velocities += half_timestep_divide_mass * force;
    // update x_{i+1/2}
    m_positions += half_timestep * m_velocities;
    // Langevin thermostat, full step
    m_velocities = factor1 * m_velocities + factor2 * randGaussian3();
    // update x_{i+1}
    m_positions += half_timestep * m_velocities;
    // collect CVs and write traj
    positionCallback(m_positions);
    Ep = potentialFunction(m_positions);
    potentialEnergyCallback(Ep);
    stepCallback(m_step);
    runCallback();
    // update f_{i+1}
    force = forceFunction(m_positions);
    forceCallback(force);
    // update v_{i+1}
    m_velocities += half_timestep_divide_mass * force;
    Ek = kineticEnergy();
    kineticEnergyCallback(Ek);
    velocityCallback(m_velocities);
  }
}
