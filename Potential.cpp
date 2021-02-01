#include "Potential.h"

potential p;

double getPotential(const double x, const double y, const double z) {
    double res = 0;
    p.meta.get(vector<double>{x, y}, res);
    return res;
}

vector<double> getGradients(const double x, const double y, const double z) {
    double potential;
    vector<double> grad(2);
    p.meta.computeEnergyGradients(vector<double>{x, y}, potential, grad);
    grad.push_back(0.0);
    return grad;
}

vector<double> getGaMDForces(const double x, const double y, const double z) {
    return p.getGaMDForces(x, y, z);
}
