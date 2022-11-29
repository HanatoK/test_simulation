#include "Simulation.h"
#include "Potential.h"
#include "Histogram.h"
#include "Bias.h"
#include <string>
#include <vector>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 400000000;

void UnbiasedSimulations1() {
  Reporter reporter(100, "XYZ_10_10.traj");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){return potential.getForces(r);},
    [&](const double3& r){return potential.getPotential(r);},
    [&](const double3& f){reporter.recordForces(f);},
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

void UnbiasedSimulations2() {
  Reporter reporter(100, "XYZ_100_10.traj");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{100.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){return potential.getForces(r);},
    [&](const double3& r){return potential.getPotential(r);},
    [&](const double3& f){reporter.recordForces(f);},
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

void UnbiasedSimulations3() {
  Reporter reporter(100, "XYZ_10_100.traj");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 100.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){return potential.getForces(r);},
    [&](const double3& r){return potential.getPotential(r);},
    [&](const double3& f){reporter.recordForces(f);},
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

void BiasedSimulations1() {
  // setup bias
  std::vector<Axis> ax{Axis(-6, 6, 240), Axis(-6, 6, 240)};
  std::vector<Axis> mtd_ax{Axis(-6, 6, 240), Axis(-6, 6, 240)};
  BiasWTMeABF2D bias(ax, mtd_ax, 0.1, 300.0*0.0019872041/(0.05*0.05), 300.0, 8.0, timestep, "bias_10_10.hills");
  Reporter reporter(100, "XYZ_10_10_b.traj");
  std::ofstream ofs_bias_traj("bias_10_10.dat");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){
      return potential.getForces(r);
    },
    [&](const double3& r){return potential.getPotential(r);},
    [&](double3& f){
      // apply the biasing force
      bias.applyBiasForce(f);
      reporter.recordForces(f);
    },
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){
      reporter.recordPositions(r);
      // compute and update CVs
      bias.updateCV(r);
      bias.updateExtendedLagrangian();
    },
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){
      bias.recordStep(step);
      reporter.recordStep(step);
      if (step % 100 == 0) {
        bias.writeTrajectory(ofs_bias_traj);
      }
      if (step % 100000000 == 0) {
        const std::string filename = "bias_10_10_step_" + std::to_string(step);
        bias.writeOutput(filename);
      }
    },
    [&](){reporter.report();});
  bias.writeOutput("bias_10_10");
}

void BiasedSimulations2() {
  // setup bias
  std::vector<Axis> ax{Axis(-6, 6, 240), Axis(-6, 6, 240)};
  std::vector<Axis> mtd_ax{Axis(-6, 6, 240), Axis(-6, 6, 240)};
  BiasWTMeABF2D bias(ax, mtd_ax, 0.1, 300.0*0.0019872041/(0.05*0.05), 300.0, 8.0, timestep, "bias_10_10.hills");
  Reporter reporter(100, "XYZ_100_10_b.traj");
  std::ofstream ofs_bias_traj("bias_100_10.dat");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{100.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){
      return potential.getForces(r);
    },
    [&](const double3& r){return potential.getPotential(r);},
    [&](double3& f){
      // apply the biasing force
      bias.applyBiasForce(f);
      reporter.recordForces(f);
    },
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){
      reporter.recordPositions(r);
      // compute and update CVs
      bias.updateCV(r);
      bias.updateExtendedLagrangian();
    },
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){
      bias.recordStep(step);
      reporter.recordStep(step);
      if (step % 100 == 0) {
        bias.writeTrajectory(ofs_bias_traj);
      }
      // write history
      if (step % 100000000 == 0) {
        const std::string filename = "bias_100_10_step_" + std::to_string(step);
        bias.writeOutput(filename);
      }
    },
    [&](){reporter.report();});
  bias.writeOutput("bias_100_10");
}

void BiasedSimulations3() {
  // setup bias
  std::vector<Axis> ax{Axis(-6, 6, 240), Axis(-6, 6, 240)};
  std::vector<Axis> mtd_ax{Axis(-6, 6, 240), Axis(-6, 6, 240)};
  BiasWTMeABF2D bias(ax, mtd_ax, 0.1, 300.0*0.0019872041/(0.05*0.05), 300.0, 8.0, timestep, "bias_10_10.hills");
  Reporter reporter(100, "XYZ_10_100_b.traj");
  std::ofstream ofs_bias_traj("bias_10_100.dat");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 100.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){
      return potential.getForces(r);
    },
    [&](const double3& r){return potential.getPotential(r);},
    [&](double3& f){
      // apply the biasing force
      bias.applyBiasForce(f);
      reporter.recordForces(f);
    },
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){
      reporter.recordPositions(r);
      // compute and update CVs
      bias.updateCV(r);
      bias.updateExtendedLagrangian();
    },
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){
      bias.recordStep(step);
      reporter.recordStep(step);
      if (step % 100 == 0) {
        bias.writeTrajectory(ofs_bias_traj);
      }
      // write history
      if (step % 100000000 == 0) {
        const std::string filename = "bias_10_100_step_" + std::to_string(step);
        bias.writeOutput(filename);
      }
    },
    [&](){reporter.report();});
  bias.writeOutput("bias_10_100");
}

int main() {
  // UnbiasedSimulations1();
  // UnbiasedSimulations2();
  // UnbiasedSimulations3();
  BiasedSimulations1();
  BiasedSimulations2();
  BiasedSimulations3();
  return 0;
}
