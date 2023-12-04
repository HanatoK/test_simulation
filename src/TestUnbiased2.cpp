#include "Simulation.h"
#include "Potential.h"
#include <string>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 300000000;

// gamma_x == 10*gamma_y
void UnbiasedSimulations2() {
  const size_t num_atoms = 1;
  AtomGroup atoms(num_atoms, std::vector<double>(num_atoms, mass));
  atoms.m_pos_x[0] = -2.0;
  atoms.m_pos_y[0] = -2.0;
  atoms.m_pos_z[0] = 0.0;
  const std::vector<double> friction_x(num_atoms, 100.0);
  const std::vector<double> friction_y(num_atoms, 10.0);
  const std::vector<double> friction_z(num_atoms, 10.0);
  Reporter reporter(100, atoms, "XYZ_100_10.traj");
  BSPotential potential(2.0, 2.2, 1.0 / (300.0 * 0.0019872041));
  // double3 frictions{100.0, 10.0, 10.0};
  Simulation simulation(atoms, 300.0);
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
      total_steps, timestep, friction_x, friction_y, friction_z,
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
      [&reporter](std::vector<double>& __restrict f_x,
                  std::vector<double>& __restrict f_y,
                  std::vector<double>& __restrict f_z) {
                    reporter.recordForces(&f_x, &f_y, &f_z);
                  },
      [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
      [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
      [&reporter](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&reporter](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&reporter](int64_t step){reporter.recordStep(step);},
      [&reporter](){reporter.report();});
}

int main() {
  UnbiasedSimulations2();
  return 0;
}
