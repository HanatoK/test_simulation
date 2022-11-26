#include "Simulation.h"
#include "Potential.h"
#include "Axis.h"
#include "Bias.h"
#include <vector>

void UnbiasedSimulations1() {
  Reporter reporter(100, "XYZ_10_10.traj");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(1.0, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    100000000, 0.5, frictions,
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
  Simulation simulation(1.0, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{100.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    100000000, 0.5, frictions,
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
  Simulation simulation(1.0, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 100.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    100000000, 0.5, frictions,
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
  std::vector<Axis> mtd_ax{Axis(-7, 7, 280), Axis(-7, 7, 280)};
  BiasWTMeABF2D bias(ax, mtd_ax, 100.0, 300.0*0.0019872041/(0.05*0.05), 300.0, 10.0, 0.5, "bias_10_10.hills");
  Reporter reporter(100, "XYZ_10_10_b.traj");
  std::ofstream ofs_bias_traj("bias_10_10.dat");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(1.0, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    100000000, 0.5, frictions,
    [&](const double3& r){
      // compute and update CVs
      bias.positionCallback(r);
      return potential.getForces(r);
    },
    [&](const double3& r){return potential.getPotential(r);},
    [&](double3& f){
      reporter.recordForces(f);
      // apply the biasing force
      // bias.applyBiasForce(f);
    },
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){
      bias.recordStep(step);
      reporter.recordStep(step);
      if (step % 100 == 0) {
        bias.writeTrajectory(ofs_bias_traj);
      }
    },
    [&](){reporter.report();});
  bias.writeOutput("bias_10_10");
}

int main() {
  // UnbiasedSimulations1();
  // UnbiasedSimulations2();
  // UnbiasedSimulations3();
  BiasedSimulations1();
  return 0;
}
