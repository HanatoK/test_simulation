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

// 10*gamma_x == gamma_y
void BiasedSimulations3() {
  // setup bias
  std::vector<Axis> ax{Axis(-7, 7, 280), Axis(-7, 7, 280)};
  std::vector<Axis> mtd_ax{Axis(-7, 7, 280), Axis(-7, 7, 280)};
  BiasWTMeABF2D bias(ax, mtd_ax, 0.1, 300.0*0.0019872041/(0.05*0.05), 300.0, 8.0, timestep);
  Reporter reporter(100, "XYZ_10_100_b.traj");
  std::ofstream ofs_bias_traj("bias_10_100.dat");
  ofs_bias_traj << "# step x y r_x r_y fb_rx fb_ry\n";
  std::ofstream ofs_hill_traj("bias_10_100.hills");
  ofs_hill_traj << "# step x y sigma_x sigma_y height\n";
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 100.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
      total_steps, timestep, frictions,
      [&](const double3& r){
        auto force =  potential.getForces(r);
        const auto restraint_force = restraintForce(r);
        force.x += restraint_force.x;
        force.y += restraint_force.y;
        force.z += restraint_force.z;
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
        f.x += bias_force[0];
        f.y += bias_force[1];
        reporter.recordForces(f);
      },
      [&](const double3& v){reporter.recordVelocities(v);},
      [&](const double3& r){
        reporter.recordPositions(r);
        // compute and update CVs
        bias.updateCV(std::vector{r.x, r.y});
      },
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){
        bias.recordStep(step);
        reporter.recordStep(step);
      },
      [&](){
        bias.writeTrajectory(ofs_bias_traj);
        bias.writeOutput("bias_10_100_step_");
        reporter.report();
      });
  bias.writeOutput("bias_10_100");
}

int main() {
  BiasedSimulations3();
  return 0;
}