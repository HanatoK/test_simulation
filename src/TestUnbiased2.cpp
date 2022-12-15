#include "Simulation.h"
#include "Potential.h"
#include <string>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 300000000;

// gamma_x == 10*gamma_y
void UnbiasedSimulations2() {
  Reporter reporter(100, "XYZ_100_10.traj");
  BSPotential potential(2.0, 2.2, 1.0 / (300.0 * 0.0019872041));
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

int main() {
  UnbiasedSimulations2();
  return 0;
}
