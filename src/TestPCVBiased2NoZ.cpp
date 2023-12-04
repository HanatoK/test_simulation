#include "Simulation.h"
#include "Potential.h"
#include "Histogram.h"
#include "Bias.h"
#include "PCV.h"
#include <fstream>
#include <string>
#include <vector>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 300000000;

// gamma_x == 10*gamma_y
void PCVBiasedSimulations2(const std::string& path_filename = "../data/path_new3.txt") {
  const size_t num_atoms = 1;
  AtomGroup atoms(num_atoms, std::vector<double>(num_atoms, mass));
  atoms.m_pos_x[0] = -2.0;
  atoms.m_pos_y[0] = -2.0;
  atoms.m_pos_z[0] = 0.0;
  // setup bias
  std::vector<Axis> ax{Axis(0.01, 0.99, 98)};
  std::vector<Axis> mtd_ax{Axis(0.01, 0.99, 98)};
  BiasWTMeABF bias(ax, mtd_ax, 0.1, 300.0 * 0.0019872041 / (0.01 * 0.01), 300.0, 8.0, timestep);
  // HarmonicWalls restraint({0.01, -0.1}, {0.99, 1.0}, {1000000.0, 100.0});
  Reporter reporter(100, atoms, "PCV_100_10_b.traj");
  // std::ofstream ofs_restraint_traj("PCV_restraint_100_10.dat");
  // ofs_restraint_traj << fmt::format("#{:>10s} {:>12s} {:>12s} {:>12s}\n", "step", "x", "y", "restraint_energy");
  std::ofstream ofs_bias_traj("PCV_bias_100_10_pcv.dat");
  ofs_bias_traj << fmt::format("#{:>10s} {:>12s} {:>12s} {:>12s}\n", "step", "s", "r_s", "fb_s");
  std::ofstream ofs_hill_traj("PCV_bias_100_10_pcv.hills");
  ofs_hill_traj << fmt::format("#{:>10s} ", "step")
                << fmt::format("{:>12s} {:>12s} ", "x", "y")
                << fmt::format("{:>12s} {:>12s} ", "sigma_x", "sigma_y")
                << fmt::format("{:>12s}\n", "hill_h");
  BSPotential potential(2.0, 2.2, 1.0 / (300.0 * 0.0019872041));
  PathCV pcv(path_filename);
  std::ofstream ofs_pcv_traj("PCV_bias_100_10_pcv.traj");
  ofs_pcv_traj << fmt::format("#{:>10s} ", "step")
               << fmt::format("{:>12s} {:>12s} ", "pcv_s", "pcv_z")
               << fmt::format("{:>12s} {:>12s} ", "x", "y")
               << fmt::format("{:>12s} {:>12s} ", "dsdx", "dsdy")
               << fmt::format("{:>12s} {:>12s}\n", "dzdx", "dzdy");
  Simulation simulation(atoms, 300.0);
  const std::vector<double> friction_x(num_atoms, 100.0);
  const std::vector<double> friction_y(num_atoms, 10.0);
  const std::vector<double> friction_z(num_atoms, 100.0);
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
      [&](std::vector<double>& __restrict f_x,
          std::vector<double>& __restrict f_y,
          std::vector<double>& __restrict f_z) {
        // apply the biasing force (backward)
        std::vector<double> bias_force(ax.size(), 0);
        // the MTD hill and other biasing forces to the extended system is updated here
        bias.updateExtendedLagrangian();
        // write the hill immediately
        bias.writeHills(ofs_hill_traj);
        bias.applyBiasForce(bias_force);
        const auto& grad_s = pcv.get_dsdx();
        const auto& grad_z = pcv.get_dzdx();
        // chain rule
        f_x[0] += bias_force[0] * grad_s[0];
        f_y[0] += bias_force[0] * grad_s[1];
        // // restraint force along s and z
        // const auto& restraint_force_sz = restraint.force();
        // // apply the restraint force along s
        // f_x[0] += restraint_force_sz[0] * grad_s[0];
        // f_y[0] += restraint_force_sz[0] * grad_s[1];
        // // apply the restraint force along z
        // f_x[0] += restraint_force_sz[1] * grad_z[0];
        // f_y[0] += restraint_force_sz[1] * grad_z[1];
        // report the total force
        reporter.recordForces(&f_x, &f_y, &f_z);
      },
      [](std::vector<double>& __restrict,
         std::vector<double>& __restrict,
         std::vector<double>& __restrict) {},
      [&](const std::vector<double>& __restrict pos_x,
          const std::vector<double>& __restrict pos_y,
          const std::vector<double>& __restrict pos_z){
        // reporter.recordPositions(r);
        // update CVs (forward) and collect CZAR
        pcv.update_value(pos_x[0], pos_y[0]);
        bias.updateCV({pcv.get_s()});
        // restraint.update_value({pcv.get_s(), pcv.get_z()});
      },
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){
        bias.recordStep(step);
        reporter.recordStep(step);
      },
      [&](){
        bias.writeTrajectory(ofs_bias_traj);
        bias.writeOutput("PCV_bias_100_10_step");
        const auto step = simulation.getStep();
        if (step % 100 == 0) {
          // ofs_restraint_traj << fmt::format("  {:>10d} {:12.5e} {:12.5e}\n",
          //                                   step, fmt::join(restraint.position(), " "), restraint.energy());
          ofs_pcv_traj << fmt::format(" {:>10d} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}\n",
                                      step, pcv.get_s(), pcv.get_z(), fmt::join(pcv.get_point(), " "),
                                      fmt::join(pcv.get_dsdx(), " "), fmt::join(pcv.get_dzdx(), " "));
        }
        reporter.report();
      });
  bias.writeOutput("PCV_bias_100_10");
}

int main(int argc, char* argv[]) {
  if (argc < 2) PCVBiasedSimulations2();
  else PCVBiasedSimulations2(argv[1]);
  return 0;
}
