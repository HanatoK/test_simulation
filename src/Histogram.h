/*
  PMFToolBox: A toolbox to analyze and post-process the output of
  potential of mean force calculations.
  Copyright (C) 2020  Haochuan Chen

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef HISTOGRAMBASE_H
#define HISTOGRAMBASE_H

// #include "base/graph.h"
// #include "base/helper.h"
#include "Common.h"

// #include <QDebug>
// #include <QFile>
// #include <QObject>
// #include <QString>
// #include <QTextStream>
// #include <QThread>

#include <cctype>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <fmt/format.h>

using std::size_t;

class HistogramBase;

class Axis {
public:
  Axis();
  Axis(double lowerBound, double upperBound, size_t bins,
       bool periodic = false);
  void setPeriodicity(bool periodic, double periodicLower,
                      double periodicUpper);
  double width() const;
  size_t bin() const;
  bool inBoundary(double x) const;
  size_t index(double x, bool *inBoundary = nullptr) const;
  double wrap(double x) const;
  std::string infoHeader() const;
  std::vector<double> getMiddlePoints() const;
  std::vector<double> getBoundaryPoints() const;
  double lowerBound() const;
  double upperBound() const;
  double setLowerBound(double newLowerBound);
  double setUpperBound(double newUpperBound);
  double setWidth(double new_width);
  double dist(double x, double reference) const;
  bool realPeriodic() const;
  bool periodic() const;
  double period() const;
  friend class HistogramBase;

private:
  double mLowerBound;
  double mUpperBound;
  size_t mBins;
  double mWidth;
  bool mPeriodic;
  double mPeriodicLowerBound;
  double mPeriodicUpperBound;
};

// QDebug operator<<(QDebug dbg, const Axis &ax);

class HistogramBase {
public:
  HistogramBase();
  virtual ~HistogramBase();
  explicit HistogramBase(const std::vector<Axis> &ax);
  virtual bool readFromStream(std::ifstream &ifs);
  virtual bool writeToStream(std::ofstream &ofs) const;
  bool isInGrid(const std::vector<double> &position) const;
  virtual std::vector<size_t> index(const std::vector<double> &position,
                                    bool *inBoundary = nullptr) const;
  virtual size_t address(const std::vector<double> &position,
                         bool *inBoundary = nullptr) const;
  virtual size_t address(const std::vector<size_t> &idx) const;
  std::vector<double> reverseAddress(size_t address,
                                     bool *inBoundary = nullptr) const;
  size_t histogramSize() const;
  size_t dimension() const;
  const std::vector<Axis> &axes() const;
  const std::vector<std::vector<double>>& pointTable() const;
  const std::vector<size_t>& pointTableAddr() const {
    return mPointTableAddr;
  }

protected:
  size_t mNdim;
  size_t mHistogramSize;
  std::vector<Axis> mAxes;
  std::vector<std::vector<double>> mPointTable;
  std::vector<size_t> mPointTableAddr;
  std::vector<size_t> mAccu;
  void fillTable();
};

// 1D histogram
template <typename T> class HistogramScalar : public virtual HistogramBase {
public:
  static_assert(std::is_arithmetic<T>::value,
                "HistogramScalar requires a scalar type!");
  HistogramScalar();
  explicit HistogramScalar(const std::vector<Axis> &ax);
  virtual ~HistogramScalar();
  virtual bool readFromStream(std::ifstream &ifs) override;
  virtual bool readFromFile(const std::string &filename);
  virtual bool writeToStream(std::ofstream &ofs) const override;
  virtual bool writeToFile(const std::string &filename) const;
  virtual T operator()(const std::vector<double> &position);
  virtual const T operator()(const std::vector<double> &position) const;
  virtual T &operator[](size_t addr);
  virtual const T &operator[](size_t addr) const;
  T sum() const;
  T minimum() const;
  const std::vector<T> &data() const;
  std::vector<T> &data();
  virtual std::vector<T> getDerivative(const std::vector<double> &pos,
                                       bool *inBoundary = nullptr) const;
  virtual bool set(const std::vector<double> &pos, const T &value);
  virtual void merge(const HistogramScalar<T> &source);

protected:
  std::vector<T> mData;
};

template <typename T> HistogramScalar<T>::HistogramScalar() : mData(0) {
  // qDebug() << "Calling" << Q_FUNC_INFO;
}

template <typename T>
HistogramScalar<T>::HistogramScalar(const std::vector<Axis> &ax)
    : HistogramBase(ax) {
  // qDebug() << "Calling" << Q_FUNC_INFO;
  mData.assign(mHistogramSize, T());
}

template <typename T> HistogramScalar<T>::~HistogramScalar() {
  // qDebug() << "Calling" << Q_FUNC_INFO;
}

template <typename T>
bool HistogramScalar<T>::readFromStream(std::ifstream &ifs) {
  // qDebug() << "Calling" << Q_FUNC_INFO;
  bool file_opened = HistogramBase::readFromStream(ifs);
  if (!file_opened)
    return file_opened;
  // read data into m_data
  std::string line;
  std::vector<double> pos(mNdim, 0);
  std::vector<std::string> tmpFields;
  mData.resize(mHistogramSize);
  while (std::getline(ifs, line)) {
    tmpFields.clear();
    splitString(line, " ", tmpFields);
    // skip blank lines
    if (tmpFields.size() == int(mNdim) + 1) {
      // skip unnecessary comment lines starting with #
      if (tmpFields.size() == 0) {
        continue;
      } else if (tmpFields[0][0] != '#') {
        for (size_t i = 0; i < mNdim; ++i) {
          pos[i] = std::stod(tmpFields[i]);
        }
        // find the position
        const size_t addr = address(pos);
        mData[addr] = std::stod(tmpFields[mNdim]);
      }
    } else {
      return false;
    }
  }
  return true;
}

template <typename T>
bool HistogramScalar<T>::readFromFile(const std::string &filename) {
  std::ifstream ifs(filename);
  if (ifs.is_open()) {
    return readFromStream(ifs);
  } else {
    std::cerr << "Failed to open file:" << filename;
    return false;
  }
}

template <typename T>
bool HistogramScalar<T>::writeToStream(std::ofstream &ofs) const {
  // qDebug() << "Calling" << Q_FUNC_INFO;
  bool file_opened = HistogramBase::writeToStream(ofs);
  if (!file_opened)
    return file_opened;
  // std::vector<double> pos(mNdim, 0);
  for (size_t i = 0; i < mHistogramSize; ++i) {
    ofs << fmt::format(" {:15.10f}", fmt::join(mPointTable[i], " "));
    // find the position
    const size_t& addr = mPointTableAddr[i];
    ofs << fmt::format(" {}", mData[addr]) << '\n';
  }
  // restore flags
  return true;
}

template <typename T>
bool HistogramScalar<T>::writeToFile(const std::string &filename) const {
  std::ofstream ofs(filename);
  if (ofs.is_open()) {
    return writeToStream(ofs);
  } else {
    return false;
  }
}

template <typename T>
T HistogramScalar<T>::operator()(const std::vector<double> &position) {
  bool inBoundary = true;
  const size_t addr = address(position, &inBoundary);
  if (inBoundary == false) {
    return T();
  } else {
    return (*this)[addr];
  }
}

template <typename T>
const T
HistogramScalar<T>::operator()(const std::vector<double> &position) const {
  bool inBoundary = true;
  const size_t addr = address(position, &inBoundary);
  if (inBoundary == false) {
    return T();
  } else {
    return (*this)[addr];
  }
}

template <typename T> T &HistogramScalar<T>::operator[](size_t addr) {
  return mData[addr];
}

template <typename T>
const T &HistogramScalar<T>::operator[](size_t addr) const {
  return mData[addr];
}

template <typename T> T HistogramScalar<T>::sum() const {
  return std::accumulate(mData.begin(), mData.end(), T(0));
}

template <typename T> T HistogramScalar<T>::minimum() const {
  const T result = *std::min_element(mData.begin(), mData.end());
  return result;
}

template <typename T> const std::vector<T> &HistogramScalar<T>::data() const {
  return mData;
}

template <typename T> std::vector<T> &HistogramScalar<T>::data() {
  return mData;
}

template <typename T>
std::vector<T> HistogramScalar<T>::getDerivative(const std::vector<double> &pos,
                                                 bool *inBoundary) const {
  size_t addr;
  addr = address(pos, inBoundary);
  if (inBoundary != nullptr) {
    if (*inBoundary == false)
      return std::vector<T>(mNdim);
  }
  const T data_this = mData[addr];
  std::vector<T> result(mNdim, T());
  for (size_t i = 0; i < mNdim; ++i) {
    const double bin_width = mAxes[i].width();
    const size_t addr_first = addr - mAccu[i] * mAxes[i].index(pos[i]) + 0;
    const size_t addr_last = addr_first + mAccu[i] * (mAxes[i].bin() - 1);
    if (addr == addr_first) {
      if (mAxes[i].realPeriodic()) {
        const T &data_next = mData[addr + mAccu[i]];
        const T &data_prev = mData[addr_last];
        result[i] = (data_next - data_prev) / (2.0 * bin_width);
      } else {
        const T &data_next = mData[addr + mAccu[i]];
        const T &data_next2 = mData[addr + mAccu[i] * 2];
        result[i] = (data_next2 * -1.0 + data_next * 4.0 - data_this * 3.0) /
                    (2.0 * bin_width);
      }
    } else if (addr == addr_last) {
      if (mAxes[i].realPeriodic()) {
        const T &data_prev = mData[addr - mAccu[i]];
        const T &data_next = mData[addr_first];
        result[i] = (data_next - data_prev) / (2.0 * bin_width);
      } else {
        const T &data_prev = mData[addr - mAccu[i]];
        const T &data_prev2 = mData[addr - mAccu[i] * 2];
        result[i] = (data_this * 3.0 - data_prev * 4.0 + data_prev2) /
                    (2.0 * bin_width);
      }
    } else {
      const T &data_prev = mData[addr - mAccu[i]];
      const T &data_next = mData[addr + mAccu[i]];
      result[i] = (data_next - data_prev) / (2.0 * bin_width);
    }
  }
  return result;
}

template <typename T>
bool HistogramScalar<T>::set(const std::vector<double> &pos, const T &value) {
  bool inBoundary = true;
  const size_t addr = address(pos, &inBoundary);
  if (inBoundary) {
    mData[addr] = value;
  }
  return inBoundary;
}

template <typename T>
void HistogramScalar<T>::merge(const HistogramScalar<T> &source) {
  for (size_t i = 0; i < source.histogramSize(); ++i) {
    bool inSourceBoundary = true;
    bool inThisBoundary = true;
    const std::vector<double> pos = source.reverseAddress(i, &inSourceBoundary);
    const size_t this_addr = this->address(pos, &inThisBoundary);
    if (inSourceBoundary && inThisBoundary) {
      this->mData[this_addr] += source[i];
    }
  }
}

// nD histogram
template <typename T> class HistogramVector : public virtual HistogramBase {
public:
  static_assert(std::is_arithmetic<T>::value,
                "HistogramVector requires a scalar type!");
  HistogramVector();
  HistogramVector(const std::vector<Axis> &, const size_t);
  virtual ~HistogramVector();
  virtual bool readFromStream(std::ifstream &ifs) override;
  bool readFromStream(std::ifstream &ifs, size_t multiplicity);
  virtual bool readFromFile(const std::string &filename);
  virtual bool writeToStream(std::ofstream &ofs) const override;
  virtual bool writeToFile(const std::string &filename) const;
  virtual std::vector<T> operator()(const std::vector<T> &) const;
  T &operator[](size_t);
  const T &operator[](size_t) const;
  size_t multiplicity() const;

protected:
  size_t mMultiplicity;
  std::vector<T> mData;
};

template <typename T>
HistogramVector<T>::HistogramVector()
    : HistogramBase(), mMultiplicity(0), mData(0) {
}

template <typename T>
HistogramVector<T>::HistogramVector(const std::vector<Axis> &ax,
                                    const size_t multiplicity)
    : HistogramBase(ax), mMultiplicity(multiplicity) {
  mData.resize(mHistogramSize * mMultiplicity, T());
}

template <typename T> HistogramVector<T>::~HistogramVector() {
}

template <typename T>
bool HistogramVector<T>::readFromStream(std::ifstream &ifs) {
  return readFromStream(ifs, 0);
}

template <typename T>
bool HistogramVector<T>::readFromStream(std::ifstream &ifs, size_t multiplicity) {
  bool file_opened = HistogramBase::readFromStream(ifs);
  if (!file_opened)
    return file_opened;
  // try to use the dimensionality as multiplicity if it is not specified
  mMultiplicity = multiplicity > 0 ? multiplicity : mNdim;
  std::string line;
  std::vector<double> pos(mNdim, 0);
  std::vector<std::string> tmpFields;
  mData.resize(mHistogramSize * mMultiplicity);
  size_t dataLines = 0;
  while (std::getline(ifs, line)) {
    tmpFields.clear();
    splitString(line, " ", tmpFields);
    // skip blank lines
    if (tmpFields.size() == static_cast<int>(mNdim + mMultiplicity)) {
      // skip unnecessary comment lines starting with #
      if (tmpFields[0][0] != '#') {
        bool ok = true;
        for (size_t i = 0; i < mNdim; ++i) {
          try {
            pos[i] = std::stod(tmpFields[i]);
          } catch (std::invalid_argument) {
            std::cerr << fmt::format("invalid conversion at line {}: {}", dataLines, line);
            ok = false;
          } catch (std::out_of_range) {
            std::cerr << fmt::format("value out of range: {}", tmpFields[i]);
            ok = false;
          }
        }
        // find the position
        const size_t addr = address(pos);
        for (size_t j = 0; j < mMultiplicity; ++j) {
          try {
            mData[addr * mMultiplicity + j] = std::stod(tmpFields[mNdim + j]);
          } catch (std::invalid_argument) {
            std::cerr << fmt::format("invalid conversion at line {}: {}", dataLines, line);
            ok = false;
          } catch (std::out_of_range) {
            std::cerr << fmt::format("value out of range: {}", tmpFields[mNdim + j]);
            ok = false;
          }
        }
        if (ok) {
          ++dataLines;
        }
      }
    }
  }
  // std::cout << ": expect " << mHistogramSize << " lines, read "
  //           << dataLines << "lines";
  return true;
}

template <typename T>
bool HistogramVector<T>::readFromFile(const std::string &filename) {
  std::ifstream ifs(filename);
  if (ifs.is_open()) {
    return readFromStream(ifs);
  } else {
    std::cerr << "Failed to open file:" << filename;
    return false;
  }
}

template <typename T>
bool HistogramVector<T>::writeToStream(std::ofstream &ofs) const {
  // qDebug() << "Calling" << Q_FUNC_INFO;
  bool file_opened = HistogramBase::writeToStream(ofs);
  if (!file_opened)
    return file_opened;
  // std::vector<double> pos(mNdim, 0);
  for (size_t i = 0; i < mHistogramSize; ++i) {
    ofs << fmt::format(" {:15.10f}", fmt::join(mPointTable[i], " "));
    const size_t& addr = mPointTableAddr[i];
    for (size_t k = 0; k < mMultiplicity; ++k) {
      ofs << fmt::format(" {:15.10f}", mData[addr * mMultiplicity + k]);
    }
    ofs << '\n';
  }
  return true;
}

template <typename T>
bool HistogramVector<T>::writeToFile(const std::string &filename) const {
  std::ofstream ofs(filename);
  if (ofs.is_open()) {
    return writeToStream(ofs);
  } else {
    return false;
  }
}

template <typename T>
std::vector<T> HistogramVector<T>::operator()(const std::vector<T> &pos) const {
  bool inBoundary = true;
  const size_t addr = address(pos, &inBoundary);
  if (inBoundary) {
    std::vector<T> result(mMultiplicity);
    std::copy_n(mData.begin() + addr * mMultiplicity, mMultiplicity,
                result.begin());
    return result;
  } else {
    return std::vector<T>(mMultiplicity, T());
  }
}

template <typename T> T &HistogramVector<T>::operator[](size_t addr_mult) {
  return mData[addr_mult];
}

template <typename T>
const T &HistogramVector<T>::operator[](size_t addr_mult) const {
  return mData[addr_mult];
}

template <typename T> size_t HistogramVector<T>::multiplicity() const {
  return mMultiplicity;
}

// class HistogramPMF : public HistogramScalar<double> {
// public:
//   HistogramPMF();
//   explicit HistogramPMF(const std::vector<Axis> &ax);
//   void toProbability(HistogramScalar<double> &probability, double kbt) const;
//   void fromProbability(const HistogramScalar<double> &probability, double kbt);
// };

// class HistogramProbability : public HistogramScalar<double> {
// public:
//   HistogramProbability();
//   explicit HistogramProbability(const std::vector<Axis> &ax);
//   virtual ~HistogramProbability();
//   void convertToFreeEnergy(double kbt);
//   HistogramProbability
//   reduceDimension(const std::vector<size_t> &new_dims) const;
// };

#endif // HISTOGRAMBASE_H
