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
// #include "Common.h"

using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
// using std::reference_wrapper;

class HistogramBase
{
public:
    HistogramBase();
    HistogramBase(const vector<Axis>& ax);
    bool isInGrid(const vector<double>& pos) const;
    bool index(const vector<double>& pos, vector<size_t>& idx) const;
    bool address(const vector<double>& pos, size_t& addr) const;
    vector<Axis> getAxes() const;
    const vector<vector<double>>& getTable() const;
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

class HistogramScalar: public HistogramBase {
public:
    HistogramScalar() {}
    HistogramScalar(const vector<Axis>& ax);
    virtual ~HistogramScalar() {}
    virtual bool set(const vector<double>& pos, double value = 1.0);
    virtual void fill(double value);
    virtual bool get(const vector<double>& pos, double& value) const;
    virtual bool add(const vector<double>& pos, double value);
    friend HistogramScalar multiply(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar add(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar minus(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar divide(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar operator+(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar operator-(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar operator*(const HistogramScalar& h1, const HistogramScalar& h2);
    friend HistogramScalar operator*(double x, const HistogramScalar& h2);
    friend HistogramScalar operator*(const HistogramScalar& h1, double x);
    friend HistogramScalar operator/(const HistogramScalar& h1, const HistogramScalar& h2);
    void applyFunction(std::function<double(double)> f);
    vector<double>& getRawData();
    const vector<double>& getRawData() const;
    virtual void writeToFile(const string& filename) const;
    virtual void readFromFile(const string& filename);
    void dump() const;
    void normalize();
    virtual HistogramScalar reduceDimension(const vector<size_t> dims) const;
protected:
    vector<double>          mValue;
};

class HistogramVector: public HistogramBase {
public:
  HistogramVector() {}
  HistogramVector(const vector<Axis>& ax, const size_t multiplicity);
  virtual ~HistogramVector() {}
  virtual bool set(const vector<double>& pos, const vector<double>& value);
  virtual bool get(const vector<double>& pos, vector<double>& value) const;
  virtual bool add(const vector<double>& pos, const vector<double>& value);
  virtual void writeToFile(const string& filename) const;
  virtual void readFromFile(const string& filename);
  vector<double>& getRawData();
  const vector<double>& getRawData() const;
protected:
  size_t                  mMultiplicity;
  vector<double>          mValue;
};

#endif
