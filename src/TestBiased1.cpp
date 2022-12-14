#include "Simulation.h"
#include "Potential.h"
#include "Histogram.h"
#include "Bias.h"
#include <fstream>
#include <string>
#include <vector>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 300000000;

// gamma_x == gamma_y
void BiasedSimulations1() {
  // setup bias
  std::vector<Axis> ax{Axis(-7, 7, 280), Axis(-7, 7, 280)};
  std::vector<Axis> mtd_ax{Axis(-7, 7, 280), Axis(-7, 7, 280)};
  BiasWTMeABF2D bias(ax, mtd_ax, 0.1, 300.0*0.0019872041/(0.05*0.05), 300.0, 8.0, timestep);
  HarmonicWalls restraint({-7.0, -7.0}, {7.0, 7.0}, 8000.0);
  Reporter reporter(100, "XYZ_10_10_b.traj");
  std::ofstream ofs_restraint_traj("restraint_10_10.dat");
  ofs_restraint_traj << "# step x y restraint_energy\n";
  std::ofstream ofs_bias_traj("bias_10_10.dat");
  ofs_bias_traj << "# step x y r_x r_y fb_rx fb_ry\n";
  std::ofstream ofs_hill_traj("bias_10_10.hills");
  ofs_hill_traj << "# step x y sigma_x sigma_y height\n";
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
      total_steps, timestep, frictions,
      [&](const double3& r){return potential.getForces(r);},
      [&](const double3& r){return potential.getPotential(r);},
      [&](double3& f){
        // apply the biasing force (backward)
        std::vector<double> bias_force(ax.size(), 0);
        // the MTD hill and other biasing forces to the extended system is updated here
        bias.updateExtendedLagrangian();
        // write the hill immediately
        bias.writeHills(ofs_hill_traj);
        bias.applyBiasForce(bias_force);
        // add restraint force
        const auto& restraint_force = restraint.force();
        f.x += bias_force[0] + restraint_force[0];
        f.y += bias_force[1] + restraint_force[1];
        reporter.recordForces(f);
      },
      [&](const double3& v){reporter.recordVelocities(v);},
      [&](const double3& r){
        reporter.recordPositions(r);
        // update CVs (forward) and collect CZAR
        bias.updateCV({r.x, r.y});
        restraint.update_value({r.x, r.y});
      },
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){
        bias.recordStep(step);
        reporter.recordStep(step);
      },
      [&](){
        bias.writeTrajectory(ofs_bias_traj);
        bias.writeOutput("bias_10_10_step");
        const auto step = simulation.getStep();
        if (step % 100 == 0) {
          ofs_restraint_traj << fmt::format("  {:>15d} {:15.10f} {:15.10f}\n",
                                            step, fmt::join(restraint.position(), " "), restraint.energy());
        }
        reporter.report();
      });
  bias.writeOutput("bias_10_10");
}

int main() {
  BiasedSimulations1();
  return 0;
}