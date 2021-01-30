#include "cv.h"
#include <functional>
#include <random>
#include <iostream>
#include <fmt/format.h>

const double boltzmann_constant = 0.001985875;

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

void Simulation::initializeVelocities() {
  m_velocities.x = randGaussian();
  m_velocities.y = randGaussian();
  m_velocities.z = randGaussian();
  // use veloctity-rescaling
  const double Ek_tilde = 0.5 * 3 * 1 / beta();
  const double alpha = std::sqrt(Ek_tilde / kineticEnergy());
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
    m_ofs_trajectory << "# x y z v_x v_y v_z f_x f_y f_z Ek\n";
  }
  void recordForces(const double3& f) {m_forces = f;}
  void recordVelocities(const double3& v) {m_velocities = v;}
  void recordPositions(const double3& r) {m_positions = r;}
  void recordKineticEnergy(const double& Ek) {m_kineticEnergy = Ek;}
//   void recordPotentialEnergy(const double& Ep) {m_potentialEnergy = Ep;}
  void report() {
    m_ofs_trajectory << fmt::format(" {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                                    " {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                                    " {:15.10f} {:15.10f}\n",
                                    m_positions.x, m_positions.y, m_positions.z,
                                    m_velocities.x, m_velocities.y, m_velocities.z,
                                    m_forces.x, m_forces.y, m_forces.z,
                                    m_kineticEnergy);
  }
};

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
    // NOTE: gradients are from MTD hills Vb(x), and
    //       since V(X) = -Vb(X) and F(X) = -dV(X)/dX,
    //       F(X) = dVb(X)/dX, multiplying -1.0 is not required.
    f.x = grad[0];
  }
  // boundary y
  if (pos.y < -9.5) {
    f.y = -1.0 * spring_constant * (pos.y - (-9.5));
  } else if (pos.y > 9.5) {
    f.y = -1.0 * spring_constant * (pos.y - 9.5);
  } else {
    // NOTE: gradients are from MTD hills Vb(x), and
    //       since V(X) = -Vb(X) and F(X) = -dV(X)/dX,
    //       F(X) = dVb(X)/dX, multiplying -1.0 is not required.
    f.y = grad[1];
  }
  return f;
}

int main() {
  Reporter reporter("XY.traj");
  Simulation simulation(12.0, 300.0, double3{-5.0, -5.0, 0.0});
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    1000000, 0.5, 10.0, forcesFromPositions,
    [&](const double3& f){reporter.recordForces(f);},
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](){reporter.report();});
  return 0;
}
