#ifndef METADYNAMICS_H
#define METADYNAMICS_H

#include "Grid.h"
#include "Axis.h"

#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <iostream>

using std::vector;
using std::string;
// using std::function;
// using std::cout;
using std::ostream;

class Hill
{
public:
    Hill() {};
    Hill(size_t ndims);
    Hill(const vector<double>& centers, const vector<double>& sigmas, double height);
    virtual ~Hill() {}
    void init(size_t ndims);
    virtual void setParameters(const vector<double>& centers, const vector<double>& sigmas, double height);
    virtual double hillEnergy(const vector<double>& pos, const vector<Axis>& axes) const;
    virtual vector<double> hillGradients(const vector<double>& pos, const vector<Axis>& axes) const;
    virtual void hillEnergyGradients(const vector<double>& pos, const vector<Axis>& axes, vector<double>& gradients, double& potential) const;
    virtual void hillGradients(const vector<double>& pos, const vector<Axis>& axes, vector<double>& gradients, bool clear_init_gradients = false) const;
    void debugOutput(ostream& os) const;
    vector<double> mCenters;
    vector<double> mSigmas;
    double mHeight;
};

// sum hills
class MetaDynamics : public HistogramScalar
{
public:
    MetaDynamics() {}
    MetaDynamics(const vector<Axis>& ax);
    void addHill(const Hill& hill);
    void compute();
    void computeEnergyGradients(const vector<double>& position, double& potential, vector<double>& gradients) const;
    void readHillsTrajectory(const string& hill_filename);
    void writeGradients(const string& filename) const;
    void writeDummyCounts(const string& filename) const;
private:
    vector<Hill> mHills;
};

#endif
