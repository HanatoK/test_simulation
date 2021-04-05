#include "../Grid.h"
#include "../Axis.h"
#include "../Tools.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main() {
  std::ifstream ifs_traj("../XY.traj");
  std::string line;
  std::vector<std::string> fields;
  const size_t skip_lines = 100000;
  const std::vector<Axis> ax{Axis(0, 10.0, 50)};
  HistogramScalar hist(ax);
  size_t current_line = 0;
  double count = 0;
  while (std::getline(ifs_traj, line)) {
    fields.clear();
    splitString(line, " ", fields);
    if (fields.size() > 0 && fields[0][0] != '#') {
      if (current_line > skip_lines) {
        const double potentialEnergy = std::stod(fields[10]);
        const double kineticEnergy = std::stod(fields[9]);
        const double totalEnergy = kineticEnergy + potentialEnergy;
        const std::vector<double> pos{totalEnergy};
        double value = 0;
        if (hist.get(pos, value)) {
          hist.set(pos, value + 1.0);
          count = count + 1.0;
        }
      }
      ++current_line;
    }
  }
  hist.applyFunction([=](double x){return x / count;});
  hist.writeToFile("dV_dist.dat");
  return 0;
}
