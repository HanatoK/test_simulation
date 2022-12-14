#include "Common.h"

const double boltzmann_constant = 0.0019872041;

void splitString(const std::string& data, const std::string& delim, std::vector<std::string>& dest) {
    size_t index = 0, new_index = 0;
    std::string tmpstr;
    while (index != data.length()) {
        new_index = data.find(delim, index);
        if (new_index != std::string::npos) tmpstr = data.substr(index, new_index - index);
        else tmpstr = data.substr(index, data.length());
        if (!tmpstr.empty()) {
            dest.push_back(tmpstr);
        }
        if (new_index == std::string::npos) break;
        index = new_index + 1;
    }
}

double beta(double temperature) {
  return 1.0 / (temperature * boltzmann_constant);
}
