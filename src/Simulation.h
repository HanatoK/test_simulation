#ifndef SIMULATION_H
#define SIMULATION_H

#include "Common.h"
#include <functional>
#include <random>
#include <iostream>
#include <fstream>
#include <fmt/format.h>

// simulation of a single atom
class Simulation {
private:
  const double conversion_factor = 418.4;
  int64_t m_step;
  double m_mass;
  double3 m_forces;
  double3 m_velocities;
  double3 m_positions;
  double m_temperatue;
  std::random_device m_random_device;
  std::mt19937 m_random_generator;
  std::normal_distribution<> m_normal_distribution;
public:
  Simulation(double mass, double temperature = 0,
             double3 positions = double3{0, 0, 0})
  : m_step(0), m_mass(mass), m_positions(positions),
    m_temperatue(temperature),
    m_random_generator(m_random_device()) {
    m_forces = double3{0, 0, 0};
    m_velocities = double3{0, 0, 0};
  }
  double randGaussian() {
    return m_normal_distribution(m_random_generator);
  }
  int64_t getStep() const {
    return m_step;
  }
  double kineticEnergy() const {
    double ek = 0;
    ek += (m_velocities.x * m_velocities.x +
           m_velocities.y * m_velocities.y +
           m_velocities.z * m_velocities.z) * m_mass * 0.5;
    // convert to kcal/mol
    return ek / conversion_factor;
  }
  double beta() const {
    return 1.0 / (m_temperatue * boltzmann_constant);
  }
  void initializeVelocities();
  void runLangevinDynamics(
    int64_t steps, const double timestep, const double friction,
    std::function<double3(double3)> forceFunction,
    std::function<double(double3)> potentialFunction,
    std::function<void(double3&)> forceCallback,
    std::function<void(double3&)> velocityCallback,
    std::function<void(double3&)> positionCallback,
    std::function<void(double&)> kineticEnergyCallback,
    std::function<void(double&)> potentialEnergyCallback,
    std::function<void(int64_t)> stepCallback,
    std::function<void()> runCallback);
  // allow anisotropic diffusivities
  void runLangevinDynamics(
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
    std::function<void()> runCallback);
};

class Reporter {
private:
  double3 m_forces;
  double3 m_velocities;
  double3 m_positions;
  double m_kineticEnergy;
  double m_potentialEnergy;
  int64_t m_stride;
  int64_t m_step;
  std::ofstream m_ofs_trajectory;
public:
  Reporter(int64_t stride, const std::string& outputname):
    m_forces{0, 0, 0},
    m_velocities{0, 0, 0},
    m_positions{0, 0, 0},
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_stride(stride),
    m_step(0) {
    m_ofs_trajectory.open(outputname.c_str());
    m_ofs_trajectory << "# step x y z v_x v_y v_z f_x f_y f_z Ek Ep\n";
  }
  void recordForces(const double3& f) {m_forces = f;}
  void recordVelocities(const double3& v) {m_velocities = v;}
  void recordPositions(const double3& r) {m_positions = r;}
  void recordKineticEnergy(const double& Ek) {m_kineticEnergy = Ek;}
  void recordPotentialEnergy(const double& Ep) {m_potentialEnergy = Ep;}
  void recordStep(const int64_t& step) {m_step = step;}
  void report() {
    if (m_step % m_stride == 0)
      m_ofs_trajectory << fmt::format(" {:>15d} {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                                      " {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                                      " {:15.10f} {:15.10f} {:15.10f}\n", m_step,
                                      m_positions.x, m_positions.y, m_positions.z,
                                      m_velocities.x, m_velocities.y, m_velocities.z,
                                      m_forces.x, m_forces.y, m_forces.z,
                                      m_kineticEnergy, m_potentialEnergy);
  }
};

#endif // SIMULATION_H
