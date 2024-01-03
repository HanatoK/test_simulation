#ifndef SIMULATION_H
#define SIMULATION_H

#include "Common.h"
#include "AtomGroup.h"
#include <functional>
#include <random>
#include <fstream>
#include <string>
#include <array>
#include <fmt/format.h>

// simulation of a single atom
class Simulation {
private:
  int64_t m_step;
  AtomGroup& m_atom_group;
  double m_temperatue;
  std::random_device m_random_device;
  std::mt19937 m_random_generator;
  std::normal_distribution<> m_normal_distribution;
public:
  Simulation(AtomGroup& atom_group, double temperature = 0)
  : m_step(0), m_atom_group(atom_group),
    m_temperatue(temperature),
    m_random_generator(m_random_device()) {
  }
  double randGaussian() {
    return m_normal_distribution(m_random_generator);
  }
  int64_t getStep() const {
    return m_step;
  }
  double kineticEnergy() const {
    double ek = 0;
    for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
      double ek_i = 0;
      ek_i += m_atom_group.m_vel_x[i] * m_atom_group.m_vel_x[i];
      ek_i += m_atom_group.m_vel_y[i] * m_atom_group.m_vel_y[i];
      ek_i += m_atom_group.m_vel_z[i] * m_atom_group.m_vel_z[i];
      ek_i *= m_atom_group.m_mass[i];
      ek += ek_i;
    }
    ek *= 0.5;
    // convert to kcal/mol
    return ek / conversion_factor;
  }
  double beta() const {
    return 1.0 / (m_temperatue * boltzmann_constant);
  }
  void initializeVelocities();
  void runLangevinDynamics(
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
    std::function<void()> runCallback);
  // allow anisotropic diffusivities
  void runLangevinDynamics(
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
    std::function<void()> runCallback);
  void runBrownianDynamics(
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
    std::function<void()> runCallback);
};

class Reporter {
private:
  // double3 m_forces;
  const std::vector<double>* m_f_x;
  const std::vector<double>* m_f_y;
  const std::vector<double>* m_f_z;
  const AtomGroup& m_atom_group;
  double m_kineticEnergy;
  double m_potentialEnergy;
  int64_t m_stride;
  int64_t m_step;
  std::ofstream m_ofs_trajectory;
public:
  Reporter(int64_t stride, const AtomGroup& atom_group, const std::string& outputname):
    m_atom_group(atom_group),
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_stride(stride),
    m_step(0) {
    m_ofs_trajectory.open(outputname.c_str());
    m_ofs_trajectory << fmt::format("#{:>10s} ", "step");
    for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
      m_ofs_trajectory << fmt::format("{:>12s} ", "x_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "y_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "z_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "v_x_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "v_y_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "v_z_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "f_x_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "f_y_" + std::to_string(i))
                       << fmt::format("{:>12s} ", "f_z_" + std::to_string(i));
    }
    m_ofs_trajectory << fmt::format("{:>12s} {:>12s}\n", "Ek", "Ep");
  }
  void recordForces(
    const std::vector<double>* f_x,
    const std::vector<double>* f_y,
    const std::vector<double>* f_z) {
    m_f_x = f_x;
    m_f_y = f_y;
    m_f_z = f_z;
  }
  // void recordAtomGroup(const AtomGroup* atom_group) {
  //   m_atom_group = atom_group;
  // }
  void recordKineticEnergy(const double& Ek) {m_kineticEnergy = Ek;}
  void recordPotentialEnergy(const double& Ep) {m_potentialEnergy = Ep;}
  void recordStep(const int64_t& step) {m_step = step;}
  void report() {
    if (m_step % m_stride == 0) {
      m_ofs_trajectory << fmt::format(" {:>10d} ", m_step);
      if (m_f_x && m_f_y && m_f_z) {
        for (size_t i = 0; i < m_atom_group.m_num_atoms; ++i) {
          m_ofs_trajectory << fmt::format(
            "{:12.5e} {:12.5e} {:12.5e} "
            "{:12.5e} {:12.5e} {:12.5e} "
            "{:12.5e} {:12.5e} {:12.5e} ",
            m_atom_group.m_pos_x[i],
            m_atom_group.m_pos_y[i],
            m_atom_group.m_pos_z[i],
            m_atom_group.m_vel_x[i],
            m_atom_group.m_vel_y[i],
            m_atom_group.m_vel_z[i],
            (*m_f_x)[i], (*m_f_y)[i], (*m_f_z)[i]);
        }
      }
      m_ofs_trajectory << fmt::format("{:12.5e} {:12.5e}\n", m_kineticEnergy, m_potentialEnergy);
    }
  }
};

#endif // SIMULATION_H
