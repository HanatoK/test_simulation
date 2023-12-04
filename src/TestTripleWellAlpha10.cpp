#include "Simulation.h"
#include "Potential.h"
#include "Bias.h"
#include <string>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 100000000;

void TripleWellAlpha10Eq() {
  const size_t num_atoms = 1;
  AtomGroup atoms(num_atoms, std::vector<double>(num_atoms, mass));
  atoms.m_pos_x[0] = -1.0;
  atoms.m_pos_y[0] = 0.0;
  atoms.m_pos_z[0] = 0.0;
  Reporter reporter(100, atoms, "TripleWellAlpha10_XYZ_100_10.traj");
  TripleWellAlpha potential(10.0);
  Simulation simulation(atoms, 300.0);
  HarmonicWalls restraint({-2.0, -1.0}, {2.0, 2.0}, {100.0, 100.0});
  // double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
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
      [&potential](const std::vector<double>& __restrict pos_x,
                   const std::vector<double>& __restrict pos_y,
                   const std::vector<double>& __restrict pos_z) {
                     return potential.getPotential(pos_x, pos_y, pos_z);
                   },
      [&](std::vector<double>& __restrict f_x,
          std::vector<double>& __restrict f_y,
          std::vector<double>& __restrict f_z){
        const auto& restraint_force = restraint.force();
        f_x[0] += restraint_force[0];
        f_y[0] += restraint_force[1];
        // report the total force
        reporter.recordForces(&f_x, &f_y, &f_z);
      },
      [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
      [&](const std::vector<double>& __restrict pos_x,
          const std::vector<double>& __restrict pos_y,
          const std::vector<double>& __restrict pos_z){
        restraint.update_value({pos_x[0], pos_y[0]});
      },
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){reporter.recordStep(step);},
      [&](){reporter.report();});
}

int main() {
  TripleWellAlpha10Eq();
  return 0;
}
