#include "Simulation.h"
#include "Common.h"

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

void Simulation::runBrownianDynamics(
    int64_t steps, const double timestep,
    const double3& frictions,
    std::function<double3(double3)> forceFunction,
    std::function<double(double3)> potentialFunction,
    std::function<void(double3&)> forceCallback,
    std::function<void(double3&)> positionCallback,
    std::function<void(double&)> potentialEnergyCallback,
    std::function<void(int64_t)> stepCallback,
    std::function<void()> runCallback
  ) {
  positionCallback(m_positions);
  double Ep = potentialFunction(m_positions);
  potentialEnergyCallback(Ep);
  stepCallback(m_step);
  runCallback();
  double3 force = forceFunction(m_positions);
  forceCallback(force);
  /* unit conversion:
    delta_t: ps
    gamma: g/mol*(ps^-1)
    potential: kcal/mol
    force: kcal/(mol*angstrom)
    r: angstrom

    force * delta_t / gamma:
    =kcal/(mol*angstrom) * ps / (g/mol*(ps^-1))
    =kcal/(angstrom) * ps*ps / g
    =4.184 * kJ * ps * ps / (g * angstrom)
    =4.184 * 1e3 *1e3 g * m^2 * s^-2 * ps^2 / (g*angstrom)
    =4.184*1e6 * (1e10)^2 * angstrom^2 * s^-2 * ps^2 / (angstrom)
    =4.184*1e26 * (1e12)^-2 ps^-2 * ps^2 * angstrom
    =4.184*1e26 * 1e-24 angstrom
    =418.4 angstrom
    */
  const auto inv_gamma = timestep / (m_mass * frictions) * conversion_factor;
  const auto factor = std::sqrt(2.0 / beta()) * double3::sqrt(inv_gamma);
  for (m_step = 1; m_step <= steps; ++m_step) {
    // update x_{i+1}
    m_positions += inv_gamma * force + factor * randGaussian3();
    positionCallback(m_positions);
    Ep = potentialFunction(m_positions);
    potentialEnergyCallback(Ep);
    stepCallback(m_step);
    runCallback();
    // update f_{i+1}
    force = forceFunction(m_positions);
    forceCallback(force);
  }
}
