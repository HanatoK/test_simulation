#ifndef GRID_H
#define GRID_H

#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <iostream>

#include "Axis.h"
#include "Tools.h"

using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::reference_wrapper;

class HistogramBase
{
public:
    HistogramBase();
    HistogramBase(const vector<Axis>& ax);
    bool isInGrid(const vector<double>& pos) const;
    bool index(const vector<double>& pos, vector<size_t>& idx) const;
    bool address(const vector<double>& pos, size_t& addr) const;
    vector<Axis> getAxes() const;
    vector<vector<double>> getTable() const;
    double getGridSize() const;
    size_t getDimension() const;
protected:
    size_t                  mNDim;
    vector<Axis>            mAxes;
    vector<size_t>          mAccu;
    size_t                  mGridSize;
    vector<vector<double>>  mPointTable;

    void fillTable();
};

class HistogramValue: public HistogramBase
{
public:
    HistogramValue() {}
    HistogramValue(const vector<Axis>& ax);
    virtual ~HistogramValue() {};
    virtual bool set(const vector<double>& pos, double value = 1.0);
    virtual void fill(double value);
    virtual bool get(const vector<double>& pos, double& value) const;
    friend HistogramValue multiply(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue add(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue minus(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue divide(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator+(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator-(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator*(const HistogramValue& h1, const HistogramValue& h2);
    friend HistogramValue operator*(double x, const HistogramValue& h2);
    friend HistogramValue operator*(const HistogramValue& h1, double x);
    friend HistogramValue operator/(const HistogramValue& h1, const HistogramValue& h2);
    void applyFunction(std::function<double(double)> f);
    vector<double>& getRawData();
    const vector<double>& getRawData() const;
    virtual void writeToFile(const string& filename) const;
    virtual void readFromFile(const string& filename);
    void dump() const;
    void normalize();
    virtual HistogramValue reduceDimension(const vector<size_t> dims) const;
protected:
    vector<double>          mValue;
};

#endif
