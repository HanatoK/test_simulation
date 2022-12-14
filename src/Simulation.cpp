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
  const double factor1x = std::exp(-1.0 * frictions.x * timestep);
  const double factor1y = std::exp(-1.0 * frictions.y * timestep);
  const double factor1z = std::exp(-1.0 * frictions.z * timestep);
  const double factor2x = std::sqrt(conversion_factor / (beta() * m_mass)) *
                         std::sqrt(1.0 - std::exp(-2.0 * frictions.x * timestep));
  const double factor2y = std::sqrt(conversion_factor / (beta() * m_mass)) *
                         std::sqrt(1.0 - std::exp(-2.0 * frictions.y * timestep));
  const double factor2z = std::sqrt(conversion_factor / (beta() * m_mass)) *
                         std::sqrt(1.0 - std::exp(-2.0 * frictions.z * timestep));
  // BAOAB
  for (m_step = 1; m_step <= steps; ++m_step) {
    // update v_{i+1/2}
    m_velocities.x += conversion_factor * 0.5 * timestep * force.x / m_mass;
    m_velocities.y += conversion_factor * 0.5 * timestep * force.y / m_mass;
    m_velocities.z += conversion_factor * 0.5 * timestep * force.z / m_mass;
    // update x_{i+1/2}
    m_positions.x += 0.5 * m_velocities.x * timestep;
    m_positions.y += 0.5 * m_velocities.y * timestep;
    m_positions.z += 0.5 * m_velocities.z * timestep;
    // Langevin thermostat, full step
    m_velocities.x = factor1x * m_velocities.x + factor2x * randGaussian();
    m_velocities.y = factor1y * m_velocities.y + factor2y * randGaussian();
    m_velocities.z = factor1z * m_velocities.z + factor2z * randGaussian();
    // update x_{i+1}
    m_positions.x += 0.5 * m_velocities.x * timestep;
    m_positions.y += 0.5 * m_velocities.y * timestep;
    m_positions.z += 0.5 * m_velocities.z * timestep;
    // collect CVs and write traj
    positionCallback(m_positions);
    Ek = kineticEnergy();
    kineticEnergyCallback(Ek);
    Ep = potentialFunction(m_positions);
    potentialEnergyCallback(Ep);
    stepCallback(m_step);
    runCallback();
    // update f_{i+1}
    force = forceFunction(m_positions);
    forceCallback(force);
    // update v_{i+1}
    m_velocities.x += conversion_factor * 0.5 * timestep * force.x / m_mass;
    m_velocities.y += conversion_factor * 0.5 * timestep * force.y / m_mass;
    m_velocities.z += conversion_factor * 0.5 * timestep * force.z / m_mass;
    velocityCallback(m_velocities);
  }
}
