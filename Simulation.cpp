#include "Simulation.h"

const double boltzmann_constant = 0.001985875;

void Simulation::initializeVelocities() {
  m_velocities.x = randGaussian();
  m_velocities.y = randGaussian();
  m_velocities.z = randGaussian();
  const double alpha = std::sqrt(1 / (beta() * m_mass));
  m_velocities.x *= alpha;
  m_velocities.y *= alpha;
  m_velocities.z *= alpha;
}

void Simulation::runLangevinDynamics(
  int64_t steps, const double timestep,
  double friction,
  std::function<double3(double3)> forceFunction,
  std::function<void(double3&)> forceCallback,
  std::function<void(double3&)> velocityCallback,
  std::function<void(double3&)> positionCallback,
  std::function<void(double&)> kineticEnergyCallback,
  std::function<void()> runCallback) {
  double3 force = forceFunction(m_positions);
  positionCallback(m_positions);
  double Ek = kineticEnergy();
  kineticEnergyCallback(Ek);
  forceCallback(force);
  velocityCallback(m_velocities);
  runCallback();
  const double factor1 = std::exp(-1.0 * friction * timestep);
  const double factor2 = std::sqrt(1.0 / (beta() * m_mass)) *
                         std::sqrt(1.0 - std::exp(-2.0 * friction * timestep));
  // BAOAB
  for (int64_t i = 0; i < steps; ++i) {
    // update v_{i+1/2}
    m_velocities.x += 0.5 * timestep * force.x / m_mass;
    m_velocities.y += 0.5 * timestep * force.y / m_mass;
    m_velocities.z += 0.5 * timestep * force.z / m_mass;
    // update x_{i+1/2}
    m_positions.x += 0.5 * m_velocities.x * timestep; 
    m_positions.y += 0.5 * m_velocities.y * timestep;
    m_positions.z += 0.5 * m_velocities.z * timestep;
    // Langevin thermostat, full step
    m_velocities.x = factor1 * m_velocities.x + factor2 * randGaussian();
    m_velocities.y = factor1 * m_velocities.y + factor2 * randGaussian();
    m_velocities.z = factor1 * m_velocities.z + factor2 * randGaussian();
    // update x_{i+1}
    m_positions.x += 0.5 * m_velocities.x * timestep; 
    m_positions.y += 0.5 * m_velocities.y * timestep;
    m_positions.z += 0.5 * m_velocities.z * timestep;
    positionCallback(m_positions);
    Ek = kineticEnergy();
    kineticEnergyCallback(Ek);
    // update f_{i+1}
    force = forceFunction(m_positions);
    forceCallback(force);
    // update v_{i+1}
    m_velocities.x += 0.5 * timestep * force.x / m_mass;
    m_velocities.y += 0.5 * timestep * force.y / m_mass;
    m_velocities.z += 0.5 * timestep * force.z / m_mass;
    velocityCallback(m_velocities);
    runCallback();
  }
}

double3 forcesFromPositions(double3 pos) {
  double3 f{0, 0, 0};
  const auto grad = getGradients(pos.x, pos.y, pos.z);
  // boundary x
  const double spring_constant = 10.0;
  if (pos.x < -9.5) {
    f.x = -1.0 * spring_constant * (pos.x - (-9.5));
  } else if (pos.x > 9.5) {
    f.x = -1.0 * spring_constant * (pos.x - 9.5);
  } else {
    f.x = -1.0 * grad[0];
  }
  // boundary y
  if (pos.y < -9.5) {
    f.y = -1.0 * spring_constant * (pos.y - (-9.5));
  } else if (pos.y > 9.5) {
    f.y = -1.0 * spring_constant * (pos.y - 9.5);
  } else {
    f.y = -1.0 * grad[1];
  }
  return f;
}
