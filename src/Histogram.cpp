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

#include "Histogram.h"

#include <fmt/ranges.h>
#include <algorithm>
#include <cmath>
#include <iterator>

HistogramBase::HistogramBase()
    : mNdim(0), mHistogramSize(0), mAxes(0), mPointTable(0), mAccu(0) {}

HistogramBase::~HistogramBase() {}

HistogramBase::HistogramBase(const std::vector<Axis> &ax)
    : mNdim(ax.size()), mAxes(ax), mAccu(mNdim) {
  if (mNdim == 0)
    return;
  mHistogramSize = 1;
  for (size_t i = 0; i < mNdim; ++i) {
    mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
    mHistogramSize *= mAxes[i].bin();
  }
  // qDebug() << "Dimensionality is " << mNdim;
  // qDebug() << "Histogram size is " << mHistogramSize;
  fillTable();
}

bool HistogramBase::readFromStream(std::ifstream &ifs) {
  std::string line;
  std::getline(ifs, line);
  std::vector<std::string> tmp;
  splitString(line, " ", tmp);
  if (tmp.size() < 2)
    return false;
  mNdim = std::stoull(tmp[1]);
  mAxes.resize(mNdim);
  mAccu.resize(mNdim);
  // now we know how many axes should be read
  for (size_t i = 0; i < mNdim; ++i) {
    std::getline(ifs, line);
    tmp.clear();
    splitString(line, " ", tmp);
    if (tmp.size() < 5)
      return false;
    // initialize each axis
    // format: # lower_bound bin_width num_bins is_periodic
    mAxes[i].mLowerBound = std::stod(tmp[1]);
    mAxes[i].mWidth = std::stod(tmp[2]);
    mAxes[i].mBins = std::stoull(tmp[3]);
    mAxes[i].mUpperBound =
        mAxes[i].mLowerBound + mAxes[i].mWidth * double(mAxes[i].mBins);
    mAxes[i].mPeriodic = (std::stoi(tmp[4]) == 0) ? false : true;
    if (mAxes[i].mPeriodic) {
      mAxes[i].mPeriodicLowerBound = mAxes[i].mLowerBound;
      mAxes[i].mPeriodicUpperBound = mAxes[i].mUpperBound;
    }
  }
  // initialize other variables
  mHistogramSize = 1;
  for (size_t i = 0; i < mNdim; ++i) {
    mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
    mHistogramSize *= mAxes[i].bin();
  }
  // initialize the table
  fillTable();
  return true;
}

bool HistogramBase::writeToStream(std::ofstream &ofs) const {
  ofs << "# " << mNdim << '\n';
  for (const auto &ax : mAxes) {
    ofs << ax.infoHeader() << '\n';
  }
  return true;
}

bool HistogramBase::isInGrid(const std::vector<double> &position) const {
  auto it_val = position.cbegin();
  auto it_ax = mAxes.cbegin();
  while (it_ax != mAxes.cend()) {
    if (!(it_ax->inBoundary(*it_val)))
      return false;
    ++it_ax;
    ++it_val;
  }
  return true;
}

std::vector<size_t> HistogramBase::index(const std::vector<double> &position,
                                         bool *inBoundary) const {
  std::vector<size_t> idx(mNdim, 0);
  for (size_t i = 0; i < mNdim; ++i) {
    if (inBoundary != nullptr) {
      idx[i] = mAxes[i].index(position[i], inBoundary);
      if (*inBoundary == false) {
        std::cerr << "Warning: position " << fmt::format("({})", fmt::join(position, ", ")) << " is not in boundary!";
        break;
      }
    } else {
      idx[i] = mAxes[i].index(position[i]);
    }
  }
  return idx;
}

size_t HistogramBase::address(const std::vector<double> &position,
                              bool *inBoundary) const {
  size_t addr = 0;
  for (size_t i = 0; i < mNdim; ++i) {
    addr += mAccu[i] * mAxes[i].index(position[i], inBoundary);
    if (inBoundary != nullptr) {
      if (*inBoundary == false) {
        break;
      }
    }
  }
  return addr;
}

size_t HistogramBase::address(const std::vector<size_t>& idx) const
{
  size_t addr = 0;
  for (size_t i = 0; i < mNdim; ++i) {
    addr += mAccu[i] * idx[i];
  }
  return addr;
}

std::vector<double> HistogramBase::reverseAddress(size_t address,
                                                  bool *inBoundary) const {
  std::vector<double> pos(mNdim, 0);
  if (address >= mHistogramSize) {
    if (inBoundary != nullptr) {
      *inBoundary = false;
    }
  } else {
    for (int i = mNdim - 1; i >= 0; --i) {
      const size_t index_i = static_cast<size_t>(
          std::floor(static_cast<double>(address) / mAccu[i]));
      pos[i] = mAxes[i].mLowerBound + (0.5 + index_i) * mAxes[i].width();
      address -= index_i * mAccu[i];
    }
    if (inBoundary != nullptr) {
      *inBoundary = true;
    }
  }
  return pos;
}

size_t HistogramBase::histogramSize() const { return mHistogramSize; }

size_t HistogramBase::dimension() const { return mNdim; }

const std::vector<Axis> &HistogramBase::axes() const { return mAxes; }

const std::vector<std::vector<double>> &HistogramBase::pointTable() const {
  return mPointTable;
}

void HistogramBase::fillTable() {
  std::vector<std::vector<double>> middlePoints(mNdim);
  for (size_t i = 0; i < mNdim; ++i) {
    middlePoints[i] = mAxes[i].getMiddlePoints();
  }
  auto tmpPointTable = std::vector<std::vector<double>>(mNdim, std::vector<double>(mHistogramSize, 0.0));
  for (size_t i = 0; i < mNdim; ++i) {
    size_t repeatAll = 1, repeatOne = 1;
    for (size_t j = i + 1; j < mNdim; ++j) {
      repeatOne *= middlePoints[j].size();
    }
    for (size_t j = 0; j < i; ++j) {
      repeatAll *= middlePoints[j].size();
    }
    const size_t in_i_sz = middlePoints[i].size();
    for (size_t l = 0; l < in_i_sz; ++l) {
      std::fill_n(tmpPointTable[i].begin() + l * repeatOne, repeatOne,
                  middlePoints[i][l]);
    }
    for (size_t k = 0; k < repeatAll - 1; ++k) {
      std::copy_n(tmpPointTable[i].begin(), repeatOne * in_i_sz,
                  tmpPointTable[i].begin() + repeatOne * in_i_sz * (k + 1));
    }
  }
  mPointTable.reserve(mHistogramSize);
  mPointTableAddr.reserve(mHistogramSize);
  for (size_t i = 0; i < mHistogramSize; ++i) {
    mPointTable.push_back(std::vector<double>(mNdim, 0));
    for (size_t j = 0; j < mNdim; ++j) {
      mPointTable[i][j] = tmpPointTable[j][i];
    }
    mPointTableAddr.push_back(address(mPointTable[i]));
  }
}

Axis::Axis()
    : mLowerBound(0.0), mUpperBound(0.0), mBins(0), mWidth(0.0),
      mPeriodic(false), mPeriodicLowerBound(0.0), mPeriodicUpperBound(0.0) {
}

Axis::Axis(double lowerBound, double upperBound, size_t bins, bool periodic)
    : mLowerBound(lowerBound), mUpperBound(upperBound), mBins(bins),
      mPeriodic(periodic) {
  mWidth = (mUpperBound - mLowerBound) / static_cast<double>(mBins);
  mPeriodicLowerBound = mLowerBound;
  mPeriodicUpperBound = mUpperBound;
}

void Axis::setPeriodicity(bool periodic, double periodicLower,
                          double periodicUpper) {
  mPeriodic = periodic;
  mPeriodicLowerBound = periodicLower;
  mPeriodicUpperBound = periodicUpper;
}

double Axis::width() const { return mWidth; }

size_t Axis::bin() const { return mBins; }

bool Axis::inBoundary(double x) const {
  x = wrap(x);
  if (x < mLowerBound || x > mUpperBound) {
    return false;
  } else {
    return true;
  }
}

size_t Axis::index(double x, bool *inBoundary) const {
  x = wrap(x);
  bool checkResult = true;
  if (inBoundary != nullptr) {
    checkResult = this->inBoundary(x);
    *inBoundary = checkResult;
  }
  if (checkResult == false) {
    return 0;
  }
  size_t idx = std::floor((x - mLowerBound) / mWidth);
  if (idx == mBins)
    --idx;
  return idx;
}

double Axis::wrap(double x) const {
  if (!mPeriodic)
    return x;
  if (x >= mPeriodicLowerBound && x <= mPeriodicUpperBound)
    return x;
  const double periodicity = mPeriodicUpperBound - mPeriodicLowerBound;
  if (x < mPeriodicLowerBound) {
    const double dist_to_lower = mPeriodicLowerBound - x;
    const int num_period_add = int(dist_to_lower / periodicity);
    const double tmp = std::abs(dist_to_lower / periodicity -
                                (std::nearbyint(dist_to_lower / periodicity)));
    if (almost_equal(tmp, 0.0)) {
      x += num_period_add * periodicity;
    } else {
      x += (num_period_add + 1) * periodicity;
    }
  }
  if (x > mPeriodicUpperBound) {
    const double dist_to_upper = x - mPeriodicUpperBound;
    const int num_period_subtract = int(dist_to_upper / periodicity);
    const double tmp = std::abs(dist_to_upper / periodicity -
                                (std::nearbyint(dist_to_upper / periodicity)));
    if (almost_equal(tmp, 0.0)) {
      x -= num_period_subtract * periodicity;
    } else {
      x -= (num_period_subtract + 1) * periodicity;
    }
  }
  return x;
}

std::string Axis::infoHeader() const {
  const int pbc = mPeriodic ? 1 : 0;
  const std::string str = fmt::format("# {:17.10e} {:17.10e} {} {}", mLowerBound, mWidth, mBins, pbc);
  return str;
}

std::vector<double> Axis::getMiddlePoints() const {
  double tmp = mLowerBound - 0.5 * mWidth;
  std::vector<double> result(mBins, 0.0);
  for (auto &i : result) {
    tmp += mWidth;
    i = tmp;
  }
  return result;
}

std::vector<double> Axis::getBoundaryPoints() const
{
  std::vector<double> result(mBins + 1);
  result[0] = mLowerBound;
  for (size_t i = 1; i < mBins + 1; ++i) {
    result[i] = result[i-1] + mWidth;
  }
  return result;
}

double Axis::lowerBound() const { return mLowerBound; }

double Axis::upperBound() const { return mUpperBound; }

double Axis::setLowerBound(double newLowerBound) {
  // keep bin width and reset lower bound
  mBins =
      mWidth > 0 ? std::nearbyintl((mUpperBound - newLowerBound) / mWidth) : 0;
  mLowerBound =
      mBins == 0 ? newLowerBound : mUpperBound - double(mBins) * mWidth;
  return mLowerBound;
}

double Axis::setUpperBound(double newUpperBound) {
  // keep bin width and reset upper bound
  mBins =
      mWidth > 0 ? std::nearbyintl((newUpperBound - mLowerBound) / mWidth) : 0;
  mUpperBound =
      mBins == 0 ? newUpperBound : mLowerBound + double(mBins) * mWidth;
  return mUpperBound;
}

double Axis::setWidth(double new_width) {
  if (new_width <= 0)
    return -1.0;
  mBins = std::nearbyintl((mUpperBound - mLowerBound) / new_width);
  mWidth = mBins == 0 ? new_width : (mUpperBound - mLowerBound) / double(mBins);
  return mWidth;
}

double Axis::dist(double x, double reference) const {
  if (!periodic()) {
    return x - reference;
  } else {
    x = wrap(x);
    reference = wrap(reference);
    const double dist = x - reference;
    if (std::abs(dist) > (period() * 0.5)) {
      if (reference > x) {
        return (dist + period());
      } else if (reference < x) {
        return (dist - period());
      } else {
        return dist;
      }
    } else {
      return dist;
    }
  }
}

bool Axis::realPeriodic() const {
  if (mPeriodic) {
    if (std::abs(mLowerBound - mPeriodicLowerBound) < mWidth &&
        std::abs(mUpperBound - mPeriodicUpperBound) < mWidth) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool Axis::periodic() const { return mPeriodic; }

double Axis::period() const {
  return mPeriodicUpperBound - mPeriodicLowerBound;
}
