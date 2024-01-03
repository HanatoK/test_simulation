#include "Potential_cosine.h"
#include "Potential_flatwells.h"
#include "Simulation.h"
#include "Bias.h"

const double timestep = 0.001;
const double mass = 12.0;
const int64_t total_steps = 100000000;

void Flatwells() {
  Reporter reporter(100, "XYZ_flatwells.traj");
  Potential_flatwells potential;
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  HarmonicWalls restraint(
    std::vector{-5.0}, std::vector{5.0}, std::vector{1000.0});
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){return potential.getForces(r);},
    [&](const double3& r){
      restraint.update_value({r.x});
      return potential.getPotential(r);
    },
    [&](double3& f){
      const auto& restraint_force = restraint.force();
      f.x += restraint_force[0];
      reporter.recordForces(f);
    },
    [&](const double3& vel){reporter.recordVelocities(vel);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

void Cosine() {
  Reporter reporter(100, "XYZ_cosine.traj");
  Potential_cosine potential;
  Simulation simulation(mass, 300.0, double3{-2.0, -2.0, 0.0});
  double3 frictions{10.0, 10.0, 10.0};
  HarmonicWalls restraint(
    std::vector{-5.0}, std::vector{5.0}, std::vector{1000.0});
  simulation.runLangevinDynamics(
    total_steps, timestep, frictions,
    [&](const double3& r){return potential.getForces(r);},
    [&](const double3& r){
      restraint.update_value({r.x});
      return potential.getPotential(r);
    },
    [&](double3& f){
      const auto& restraint_force = restraint.force();
      f.x += restraint_force[0];
      reporter.recordForces(f);
    },
    [&](const double3& vel){reporter.recordVelocities(vel);},
    [&](const double3& r){reporter.recordPositions(r);},
    [&](const double& Ep){reporter.recordPotentialEnergy(Ep);},
    [&](const double& Ek){reporter.recordKineticEnergy(Ek);},
    [&](int64_t step){reporter.recordStep(step);},
    [&](){reporter.report();});
}

int main() {
  Flatwells();
  Cosine();
  return 0;
};
