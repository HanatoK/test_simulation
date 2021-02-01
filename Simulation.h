#ifndef SIMULATION_H
#define SIMULATION_H

#include "Potential.h"
#include <functional>
#include <random>
#include <iostream>
#include <fmt/format.h>

extern const double boltzmann_constant;

struct double3 {
  double x;
  double y;
  double z;
}; 

// simulation of a single atom
class Simulation {
private:
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
  : m_mass(mass), m_positions(positions), m_temperatue(temperature),
    m_random_generator(m_random_device()) {
    m_forces = double3{0, 0, 0};
    m_velocities = double3{0, 0, 0};
  }
  double randGaussian() {
    return m_normal_distribution(m_random_generator);
  }
  double kineticEnergy() const {
    double ek = 0;
    ek += (m_velocities.x * m_velocities.x +
           m_velocities.y * m_velocities.y +
           m_velocities.z * m_velocities.z) * m_mass * 0.5;
    return ek;
  }
  double beta() const {
    return 1.0 / (m_temperatue * boltzmann_constant);
  }
  void initializeVelocities();
  void runLangevinDynamics(
    int64_t steps, const double timestep, double friction,
    std::function<double3(double3)> forceFunction,
    std::function<void(double3&)> forceCallback,
    std::function<void(double3&)> velocityCallback,
    std::function<void(double3&)> positionCallback,
    std::function<void(double&)> kineticEnergyCallback,
    std::function<void()> runCallback);
};

class Reporter {
private:
  double3 m_forces;
  double3 m_velocities;
  double3 m_positions;
  double m_kineticEnergy;
  double m_potentialEnergy;
  std::ofstream m_ofs_trajectory;
public:
  Reporter(const std::string& outputname):
    m_forces{0, 0, 0},
    m_velocities{0, 0, 0},
    m_positions{0, 0, 0},
    m_kineticEnergy(0),
    m_potentialEnergy(0) {
    m_ofs_trajectory.open(outputname.c_str());
    m_ofs_trajectory << "# x y z v_x v_y v_z f_x f_y f_z Ek Ep\n";
  }
  void recordForces(const double3& f) {m_forces = f;}
  void recordVelocities(const double3& v) {m_velocities = v;}
  void recordPositions(const double3& r) {m_positions = r;}
  void recordKineticEnergy(const double& Ek) {m_kineticEnergy = Ek;}
  void recordPotentialEnergy() {
    m_potentialEnergy = getPotential(m_positions.x, m_positions.y, m_positions.z);
  }
  void report() {
    recordPotentialEnergy();
    m_ofs_trajectory << fmt::format(" {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                                    " {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                                    " {:15.10f} {:15.10f} {:15.10f}\n",
                                    m_positions.x, m_positions.y, m_positions.z,
                                    m_velocities.x, m_velocities.y, m_velocities.z,
                                    m_forces.x, m_forces.y, m_forces.z,
                                    m_kineticEnergy, m_potentialEnergy);
  }
};

double3 forcesFromPositions(double3 pos);

#endif // SIMULATION_H
