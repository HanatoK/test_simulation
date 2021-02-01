#ifndef AXIS_H
#define AXIS_H
#include <cstddef>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <iostream>

using std::string;
using std::vector;

class Axis
{
public:
    Axis();
    Axis(double lowerBound, double upperBound, size_t bin, bool periodic = false);
    bool isInitialized() const;
    double width() const;
    size_t bin() const;
    bool isInBoundary(double value) const;
    bool isPeriodic() const;
    bool wrap(double& x) const;
    bool index(double x, size_t& idx) const;
    size_t index(double x) const;
    double distance(double start, double end) const;
    string infoHeader() const;
    vector<double> middlePoint() const;
    friend bool operator==(const Axis& lhs, const Axis& rhs);
private:
    double          mLowerBound;
    double          mUpperBound;
    size_t          mBin;
    double          mWidth;
    bool            mPeriodic;
};

#endif
