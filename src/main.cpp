#include "Simulation.h"
#include "Potential.h"

int main() {
  Reporter reporter(100, "XYZ.traj");
  BSPotential potential(2.0, 2.0, 1.0 / (300.0 * 0.0019872041));
  Simulation simulation(1.0, 300.0, double3{-2.0, -2.0, 0.0});
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    100000000, 0.5, 10.0,
    [&](const double3& r){return potential.getForces(r);},
    [&](const double3& r){return potential.getPotential(r);},
    [&](const double3& f){reporter.recordForces(f);},
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
  return 0;
}
