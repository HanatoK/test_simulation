#include "Simulation.h"

int main() {
  Reporter reporter("XY.traj");
  Simulation simulation(12.0, 300.0, double3{-5.0, -5.0, 0.0});
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
    1000000, 0.5, 10.0, forcesFromPositions,
    [&](const double3& f){reporter.recordForces(f);},
    [&](const double3& v){reporter.recordVelocities(v);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](){reporter.report();});
  return 0;
} 
