#include "Simulation.h"
#include "Potential.h"
#include <string>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 1000000000;

void UnbiasedSimulationsMB() {
  const size_t num_atoms = 1;
  AtomGroup atoms(num_atoms, std::vector<double>(num_atoms, mass));
  atoms.m_pos_x[0] = -0.5;
  atoms.m_pos_y[0] = 1.5;
  atoms.m_pos_z[0] = 0.0;
  Reporter reporter(100, atoms, "MB_XYZ_10_10.traj");
  MBPotential potential;
  Simulation simulation(atoms, 300.0);
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
      total_steps, timestep, 100.0,
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
        reporter.recordForces(&f_x, &f_y, &f_z);
      },
      [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
      [&](const std::vector<double>& __restrict pos_x,
          const std::vector<double>& __restrict pos_y,
          const std::vector<double>& __restrict pos_z) {},
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){reporter.recordStep(step);},
      [&](){reporter.report();});
}

int main() {
  UnbiasedSimulationsMB();
  return 0;
}
