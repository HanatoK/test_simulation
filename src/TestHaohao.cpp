#include "Potential_cosine.h"
#include "Potential_flatwells.h"
#include "Simulation.h"
#include "Bias.h"

const double timestep = 0.001;
const double mass = 12.0;
const int64_t total_steps = 100000000;

void Flatwells() {
  const size_t num_atoms = 1;
  AtomGroup atoms(num_atoms, std::vector<double>(num_atoms, mass));
  atoms.m_pos_x[0] = -2.0;
  atoms.m_pos_y[0] = -2.0;
  atoms.m_pos_z[0] = 0.0;
  Reporter reporter(100, atoms, "XYZ_flatwells.traj");
  Potential_flatwells potential;
  Simulation simulation(atoms, 300.0);
  HarmonicWalls restraint(
    std::vector{-5.0}, std::vector{5.0}, std::vector{1000.0});
  simulation.runLangevinDynamics(
    total_steps, timestep, 10.0,
    [&potential](const std::vector<double>& __restrict pos_x,
                 const std::vector<double>& __restrict pos_y,
                 const std::vector<double>& __restrict pos_z,
                 std::vector<double>& __restrict f_x,
                 std::vector<double>& __restrict f_y,
                 std::vector<double>& __restrict f_z) {
                   potential.getForces(pos_x, pos_y, pos_z, f_x, f_y, f_z);
                 },
    [&](const std::vector<double>& __restrict pos_x,
        const std::vector<double>& __restrict pos_y,
        const std::vector<double>& __restrict pos_z){
      restraint.update_value({pos_x[0]});
      return potential.getPotential(pos_x, pos_y, pos_z);
    },
    [&](std::vector<double>& __restrict f_x,
        std::vector<double>& __restrict f_y,
        std::vector<double>& __restrict f_z){
      const auto& restraint_force = restraint.force();
      f_x[0] += restraint_force[0];
      reporter.recordForces(&f_x, &f_y, &f_z);
    },
    [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
    [](std::vector<double>& __restrict,
        std::vector<double>& __restrict,
        std::vector<double>& __restrict) {},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

void Cosine() {
  const size_t num_atoms = 1;
  AtomGroup atoms(num_atoms, std::vector<double>(num_atoms, mass));
  atoms.m_pos_x[0] = -2.0;
  atoms.m_pos_y[0] = -2.0;
  atoms.m_pos_z[0] = 0.0;
  Reporter reporter(100, atoms, "XYZ_cosine.traj");
  Potential_cosine potential;
  Simulation simulation(atoms, 300.0);
  HarmonicWalls restraint(
    std::vector{-5.0}, std::vector{5.0}, std::vector{1000.0});
  simulation.runLangevinDynamics(
    total_steps, timestep, 10.0,
    [&potential](const std::vector<double>& __restrict pos_x,
                 const std::vector<double>& __restrict pos_y,
                 const std::vector<double>& __restrict pos_z,
                 std::vector<double>& __restrict f_x,
                 std::vector<double>& __restrict f_y,
                 std::vector<double>& __restrict f_z) {
                   potential.getForces(pos_x, pos_y, pos_z, f_x, f_y, f_z);
                 },
    [&](const std::vector<double>& __restrict pos_x,
        const std::vector<double>& __restrict pos_y,
        const std::vector<double>& __restrict pos_z){
      restraint.update_value({pos_x[0]});
      return potential.getPotential(pos_x, pos_y, pos_z);
    },
    [&](std::vector<double>& __restrict f_x,
        std::vector<double>& __restrict f_y,
        std::vector<double>& __restrict f_z){
      const auto& restraint_force = restraint.force();
      f_x[0] += restraint_force[0];
      reporter.recordForces(&f_x, &f_y, &f_z);
    },
    [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
    [](std::vector<double>& __restrict,
        std::vector<double>& __restrict,
        std::vector<double>& __restrict) {},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

int main() {
  Flatwells();
  Cosine();
  return 0;
};
