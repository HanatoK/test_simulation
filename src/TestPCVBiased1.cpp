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

// gamma_x == gamma_y
void PCVBiasedSimulations1(const std::string& path_filename = "../data/path_new3.txt") {
  // setup bias
  std::vector<Axis> ax{Axis(0.01, 0.99, 98)};
  std::vector<Axis> mtd_ax{Axis(0.01, 0.99, 98)};
  BiasWTMeABF bias(ax, mtd_ax, 0.1, 300.0 * 0.0019872041 / (0.01 * 0.01), 300.0, 8.0, timestep);
  HarmonicWalls restraint({0.01, -0.1}, {0.99, 1.0}, {1000000.0, 100.0});
  Reporter reporter(100, "PCV_10_10_b.traj");
  std::ofstream ofs_restraint_traj("PCV_restraint_10_10.dat");
  ofs_restraint_traj << "# step s z restraint_energy\n";
  std::ofstream ofs_bias_traj("PCV_bias_10_10_pcv.dat");
  ofs_bias_traj << "# step s r_s fb_s\n";
  std::ofstream ofs_hill_traj("PCV_bias_10_10_pcv.hills");
  ofs_hill_traj << "# step s sigma_s height\n";
  BSPotential potential(2.0, 2.2, 1.0 / (300.0 * 0.0019872041));
  PathCV pcv(path_filename);
  std::ofstream ofs_pcv_traj("PCV_bias_10_10_pcv.traj");
  ofs_pcv_traj << "# step pcv_s pcv_z x y dsdx dsdy dzdx dzdy\n";
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 100.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
      total_steps, timestep, frictions,
      [&](const double3& r){
        auto force =  potential.getForces(r);
        return force;
      },
      [&](const double3& r){return potential.getPotential(r);},
      [&](double3& f){
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
        f.x += bias_force[0] * grad_s[0];
        f.y += bias_force[0] * grad_s[1];
        // restraint force along s and z
        const auto& restraint_force_sz = restraint.force();
        // apply the restraint force along s
        f.x += restraint_force_sz[0] * grad_s[0];
        f.y += restraint_force_sz[0] * grad_s[1];
        // apply the restraint force along z
        f.x += restraint_force_sz[1] * grad_z[0];
        f.y += restraint_force_sz[1] * grad_z[1];
        // report the total force
        reporter.recordForces(f);
      },
      [&](const double3& v){reporter.recordVelocities(v);},
      [&](const double3& r){
        reporter.recordPositions(r);
        // update CVs (forward) and collect CZAR
        pcv.update_value(r.x, r.y);
        bias.updateCV({pcv.get_s()});
        restraint.update_value({pcv.get_s(), pcv.get_z()});
      },
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){
        bias.recordStep(step);
        reporter.recordStep(step);
      },
      [&](){
        bias.writeTrajectory(ofs_bias_traj);
        bias.writeOutput("PCV_bias_10_10_step");
        const auto step = simulation.getStep();
        if (step % 100 == 0) {
          ofs_restraint_traj << fmt::format("  {:>15d} {:15.10f} {:15.10f}\n",
                                            step, fmt::join(restraint.position(), " "), restraint.energy());
          ofs_pcv_traj << fmt::format(" {:>15d} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n",
                                      step, pcv.get_s(), pcv.get_z(), fmt::join(pcv.get_point(), " "),
                                      fmt::join(pcv.get_dsdx(), " "), fmt::join(pcv.get_dzdx(), " "));
        }
        reporter.report();
      });
  bias.writeOutput("PCV_bias_10_10");
}

int main(int argc, char* argv[]) {
  if (argc < 2) PCVBiasedSimulations1();
  else PCVBiasedSimulations1(argv[1]);
  return 0;
}