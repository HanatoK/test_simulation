#include "Simulation.h"
#include "Potential.h"
#include "Bias.h"
#include <string>

const double timestep = 0.0005;
const double mass = 12.0;
const int64_t total_steps = 100000000;

void TripleWellAlpha10Eq() {
  Reporter reporter(100, "TripleWellAlpha10_XYZ_100_10.traj");
  TripleWellAlpha potential(10.0);
  Simulation simulation(mass, 300.0, double3{-1.0, 0.0, 0.0});
  HarmonicWalls restraint({-2.0, -1.0}, {2.0, 2.0}, {100.0, 100.0});
  double3 frictions{10.0, 10.0, 10.0};
  simulation.initializeVelocities();
  simulation.runLangevinDynamics(
      total_steps, timestep, frictions,
      [&](const double3& r){return potential.getForces(r);},
      [&](const double3& r){return potential.getPotential(r);},
      [&](double3& f){
        reporter.recordForces(f);
        const auto& restraint_force = restraint.force();
        f.x += restraint_force[0];
        f.y += restraint_force[1];
      },
      [&](const double3& v){reporter.recordVelocities(v);},
      [&](const double3& r){
        reporter.recordPositions(r);
        restraint.update_value({r.x, r.y});
      },
      [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
      [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
      [&](int64_t step){reporter.recordStep(step);},
      [&](){reporter.report();});
}

int main() {
  TripleWellAlpha10Eq();
  return 0;
}
