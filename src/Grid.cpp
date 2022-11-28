#include "Grid.h"
#include "Common.h"

HistogramBase::HistogramBase():
mNDim(0), mAxes(0), mAccu(0), mGridSize(0), mPointTable(0) {}

HistogramBase::HistogramBase(const vector<Axis>& ax):
mNDim(ax.size()), mAxes(ax), mAccu(mNDim) {
    mAccu[0] = 1;
    mGridSize = 1;
    for (size_t i = 0; i < mNDim; ++i) {
        mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
        mGridSize *= mAxes[i].bin();
    }
    fillTable();
}

void HistogramBase::fillTable() {
    vector<vector<double>> middlePoint(mNDim);
    for (size_t i = 0; i < mNDim; ++i) {
        middlePoint[i] = mAxes[i].middlePoint();
    }
    mPointTable.assign(mNDim, vector<double>(mGridSize, 0.0));
    for (size_t i = 0; i < mNDim; ++i) {
        size_t repeatAll = 1, repeatOne = 1;
        for (size_t j = i + 1; j < mNDim; ++j) {
            repeatOne *= middlePoint[j].size();
        }
        for (size_t j = 0; j < i; ++j) {
            repeatAll *= middlePoint[j].size();
        }
        const size_t in_i_sz = middlePoint[i].size();
        for (size_t l = 0; l < in_i_sz; ++l) {
            std::fill_n(begin(mPointTable[i]) + l * repeatOne, repeatOne, middlePoint[i][l]);
        }
        for (size_t k = 0; k < repeatAll - 1; ++k) {
        std::copy_n(begin(mPointTable[i]), repeatOne * in_i_sz,
                    begin(mPointTable[i]) + repeatOne * in_i_sz * (k + 1));
        }
    }
}

bool HistogramBase::isInGrid(const vector<double>& pos) const {
    auto it_val = pos.cbegin();
    auto it_ax = mAxes.cbegin();
    while (it_ax != mAxes.cend()) {
        if (!(it_ax->isInBoundary((*it_val)))) return false;
        ++it_ax;
        ++it_val;
    }
    return true;
}

bool HistogramBase::index(const vector<double>& pos, vector<size_t>& idx) const {
    idx.assign(mNDim, 0);
    auto it_val = pos.cbegin();
    auto it_ax = mAxes.cbegin();
    auto it_idx = idx.begin();
    while (it_idx != idx.end()) {
        if (!(it_ax->index((*it_val), (*it_idx)))) return false;
    }
    return true;
}

bool HistogramBase::address(const vector<double>& pos, size_t& addr) const {
    addr = 0;
    for (size_t i = 0; i < mNDim; ++i) {
        size_t idx_i = 0;
        if (!(mAxes[i].index(pos[i], idx_i))) return false;
        addr += idx_i * mAccu[i];
//         std::cout << "pos[" << i << "] = " << pos[i] << " idx["<<i<<"] = " << idx_i << " ";
    }
//     std::cout << std::endl;
    return true;
}

vector<Axis> HistogramBase::getAxes() const {
    return mAxes;
}

const vector<vector<double>>& HistogramBase::getTable() const {
    return mPointTable;
}

double HistogramBase::getGridSize() const {
    return mGridSize;
}

size_t HistogramBase::getDimension() const {
    return mNDim;
}

HistogramScalar::HistogramScalar(const vector<Axis>& ax): HistogramBase(ax) {
    mValue.assign(mGridSize, 0.0);
}

bool HistogramScalar::set(const vector<double>& pos, double value) {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        mValue[addr] = value;
    }
    return inGrid;
}

bool HistogramScalar::add(const vector<double>& pos, double value) {
  size_t addr = 0;
  bool inGrid = address(pos, addr);
  if (inGrid) {
    mValue[addr] += value;
  }
  return inGrid;
}

void HistogramScalar::fill(double value) {
    mValue.assign(mGridSize, value);
}

bool HistogramScalar::get(const vector<double>& pos, double& value) const {
    size_t addr = 0;
    bool inGrid = address(pos, addr);
    if (inGrid) {
        value = mValue[addr];
    }
    return inGrid;
}

void HistogramScalar::writeToFile(const string& filename) const {
    ofstream ofs_histo(filename.c_str());
    vector<double> pos(mNDim, 0.0);
    double val = 0;
    ofs_histo << "# " << mNDim << '\n';
    for (size_t j = 0; j < mNDim; ++j) {
        ofs_histo << mAxes[j].infoHeader() << '\n';
    }
    ofs_histo.setf(std::ios::fixed);
    ofs_histo << std::setprecision(7);
    for (size_t i = 0; i < mGridSize; ++i) {
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
            ofs_histo << pos[j] << ' ';
        }
        get(pos, val);
        ofs_histo << val << ' ';
        ofs_histo << '\n';
    }
}

void HistogramScalar::readFromFile(const string& filename) {
    ifstream ifs_histo(filename.c_str());
    string line;
    string token{" "};
    vector<string> fields;
    // Parse first line
    std::getline(ifs_histo, line);
    splitString(line, token, fields);
    if (fields[0].compare("#") != 0) {
        std::cerr << "Histogram file reads error!" << std::endl;
        std::abort();
    } else {
        mNDim = std::stoul(fields[1]);
    }
    fields.clear();
    // Parse axes
    mAxes.clear();
    mAxes.resize(mNDim);
    for (size_t i = 0; i < mNDim; ++i) {
        std::getline(ifs_histo, line);
        splitString(line, token, fields);
        double lower, width, upper;
        size_t bins;
        bool periodic = false;
        if (fields[0].compare("#") != 0) {
            std::cerr << "Histogram file reads error!" << std::endl;
            std::abort();
        } else {
            lower = std::stod(fields[1]);
            width = std::stod(fields[2]);
            bins = std::stoul(fields[3]);
            int p = std::stoi(fields[4]);
            upper = lower + double(bins) * width;
            periodic = (p != 0) ? true : false;
            mAxes[i] = Axis(lower, upper, bins, periodic);
        }
        fields.clear();
    }
    // Initialize mAccu
    mAccu.resize(mNDim);
    mAccu[0] = 1;
    mGridSize = 1;
    for (size_t i = 0; i < mNDim; ++i) {
        mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
        mGridSize *= mAxes[i].bin();
    }
    mValue.assign(mGridSize, 0.0);
    // Initialize table
    fillTable();
    vector<double> pos(mNDim, 0);
    size_t data_count = 0;
    while(std::getline(ifs_histo, line)) {
        splitString(line, token, fields);
        if (fields.empty()) {
            continue;
        }
        if (fields[0].compare("#") != 0) {
            if (fields.size() != (mNDim + 1)) {
                std::cerr << "Histogram file reads error!" << std::endl;
                std::abort();
            }
            double value;
            for (size_t j = 0; j < mNDim; ++j) {
                pos[j] =std::stod(fields[j]);
            }
            value = std::stod(fields[mNDim]);
            set(pos, value);
            ++data_count;
            fields.clear();
        }
    }
    if (data_count != mGridSize) {
        std::cerr << "Histogram file reads error!" << std::endl;
    }
}

void HistogramScalar::dump() const {
    using std::cout;
    std::ios_base::fmtflags f(cout.flags());
    vector<double> pos(mNDim, 0.0);
    double val = 0;
    cout.setf(std::ios::fixed);
    cout << std::setprecision(7);
    for (size_t i = 0; i < mGridSize; ++i) {
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
            cout << pos[j] << ' ';
        }
        get(pos, val);
        for (size_t j = 0; j < mNDim; ++j) {
            cout << val << ' ';
        }
        cout << '\n';
    }
    cout << std::flush;
    cout.flags(f);
}

HistogramScalar minus(const HistogramScalar& h1, const HistogramScalar& h2) {
    HistogramScalar h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = h1.mValue[i] - h2.mValue[i];
    }
    return h3;
}

HistogramScalar multiply(const HistogramScalar& h1, const HistogramScalar& h2) {
    HistogramScalar h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = h1.mValue[i] * h2.mValue[i];
    }
    return h3;
}

HistogramScalar add(const HistogramScalar& h1, const HistogramScalar& h2) {
    HistogramScalar h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = h1.mValue[i] + h2.mValue[i];
    }
    return h3;
}

HistogramScalar divide(const HistogramScalar& h1, const HistogramScalar& h2) {
    HistogramScalar h3(h1.mAxes);
    // Assume h1 and h2 have the same axes.
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        if (h2.mValue[i] != 0) {
            h3.mValue[i] = h1.mValue[i] / h2.mValue[i];
        }
    }
    return h3;
}

HistogramScalar operator+(const HistogramScalar& h1, const HistogramScalar& h2) {
    return add(h1, h2);
}

HistogramScalar operator-(const HistogramScalar& h1, const HistogramScalar& h2) {
    return minus(h1, h2);
}

HistogramScalar operator*(const HistogramScalar& h1, const HistogramScalar& h2) {
    return multiply(h1, h2);
}

HistogramScalar operator*(double x, const HistogramScalar& h2) {
    HistogramScalar h3(h2.mAxes);
    for (size_t i = 0; i < h3.mGridSize; ++i) {
        h3.mValue[i] = x * h2.mValue[i];
    }
    return h3;
}

HistogramScalar operator*(const HistogramScalar& h1, double x) {
    return x * h1;
}

HistogramScalar operator/(const HistogramScalar& h1, const HistogramScalar& h2) {
    return divide(h1, h2);
}

void HistogramScalar::applyFunction(std::function<double(double)> f) {
    for (auto it = mValue.begin(); it != mValue.end(); ++it) {
        (*it) = f(*it);
    }
}

vector<double>& HistogramScalar::getRawData() {
    return mValue;
}

const vector<double>& HistogramScalar::getRawData() const {
    const vector<double>& ret = mValue;
    return ret;
}

void HistogramScalar::normalize() {
    double factor = std::accumulate(mValue.begin(), mValue.end(), 0.0);
    if (factor > 0) {
        applyFunction([factor](double x){return x / factor;});
    }
}

HistogramScalar HistogramScalar::reduceDimension(const vector<size_t> new_dims) const {
    vector<Axis> new_ax;
    for (size_t i = 0; i < new_dims.size(); ++i) {
        new_ax.push_back(mAxes.at(new_dims[i]));
    }
    HistogramScalar new_hist(new_ax);
    vector<double> pos(mNDim, 0.0);
    vector<double> new_pos(new_hist.getDimension(), 0.0);
    for (size_t i = 0; i < mGridSize; ++i) {
        double val = 0;
        double new_val = 0;
        for (size_t j = 0; j < mNDim; ++j) {
            pos[j] = mPointTable[j][i];
        }
        for (size_t k = 0; k < new_hist.getDimension(); ++k) {
            new_pos[k] = pos[new_dims[k]];
        }
        get(pos, val);
        new_hist.get(new_pos, new_val);
        new_hist.set(new_pos, new_val + val);
    }
    return new_hist;
}

HistogramVector::HistogramVector(const vector<Axis>& ax,
                                 const size_t multiplicity):
  HistogramBase(ax), mMultiplicity(multiplicity) {
  mValue.assign(mGridSize * mMultiplicity, 0.0);
}

bool HistogramVector::set(const vector<double>& pos, const vector<double>& value) {
  size_t addr = 0;
  bool inGrid = address(pos, addr);
  if (inGrid) {
    for (size_t i = 0; i < mMultiplicity; ++i) {
      mValue[addr*mMultiplicity+i] = value[i];
    }
  }
  return inGrid;
}

bool HistogramVector::get(const vector<double>& pos, vector<double>& value) const {
  size_t addr = 0;
  bool inGrid = address(pos, addr);
  if (inGrid) {
    for (size_t i = 0; i < mMultiplicity; ++i) {
      value[i] = mValue[addr*mMultiplicity+i];
    }
  }
  return inGrid;
}

bool HistogramVector::add(const vector<double>& pos, const vector<double>& value) {
  size_t addr = 0;
  bool inGrid = address(pos, addr);
  if (inGrid) {
    for (size_t i = 0; i < mMultiplicity; ++i) {
      mValue[addr*mMultiplicity+i] += value[i];
    }
  }
  return inGrid;
}

void HistogramVector::writeToFile(const string& filename) const {
  ofstream ofs_histo(filename.c_str());
  vector<double> pos(mNDim, 0.0);
  vector<double> val(mMultiplicity);
  ofs_histo << "# " << mNDim << '\n';
  for (size_t j = 0; j < mNDim; ++j) {
    ofs_histo << mAxes[j].infoHeader() << '\n';
  }
  ofs_histo.setf(std::ios::fixed);
  ofs_histo << std::setprecision(7);
  for (size_t i = 0; i < mGridSize; ++i) {
    for (size_t j = 0; j < mNDim; ++j) {
      pos[j] = mPointTable[j][i];
      ofs_histo << pos[j] << ' ';
    }
    get(pos, val);
    for (size_t k = 0; k < mMultiplicity; ++k) {
      ofs_histo << val[k] << ' ';
    }
    ofs_histo << '\n';
  }
}

void HistogramVector::readFromFile(const string& filename) {
  ifstream ifs_histo(filename.c_str());
  string line;
  string token{" "};
  vector<string> fields;
  // Parse first line
  std::getline(ifs_histo, line);
  splitString(line, token, fields);
  if (fields[0].compare("#") != 0) {
    std::cerr << "Histogram file reads error!" << std::endl;
    std::abort();
  } else {
    mNDim = std::stoul(fields[1]);
  }
  fields.clear();
  // Parse axes
  mAxes.clear();
  mAxes.resize(mNDim);
  for (size_t i = 0; i < mNDim; ++i) {
    std::getline(ifs_histo, line);
    splitString(line, token, fields);
    double lower, width, upper;
    size_t bins;
    bool periodic = false;
    if (fields[0].compare("#") != 0) {
      std::cerr << "Histogram file reads error!" << std::endl;
      std::abort();
    } else {
      lower = std::stod(fields[1]);
      width = std::stod(fields[2]);
      bins = std::stoul(fields[3]);
      int p = std::stoi(fields[4]);
      upper = lower + double(bins) * width;
      periodic = (p != 0) ? true : false;
      mAxes[i] = Axis(lower, upper, bins, periodic);
    }
    fields.clear();
  }
  // Initialize mAccu
  mAccu.resize(mNDim);
  mAccu[0] = 1;
  mGridSize = 1;
  for (size_t i = 0; i < mNDim; ++i) {
    mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
    mGridSize *= mAxes[i].bin();
  }
  // Initialize table
  fillTable();
  vector<double> pos(mNDim, 0);
  size_t data_count = 0;
  bool first_time = true;
  while(std::getline(ifs_histo, line)) {
    splitString(line, token, fields);
    if (fields.empty()) {
      continue;
    }
    if (fields[0].compare("#") != 0) {
      if (first_time) {
        mMultiplicity = fields.size() - mNDim;
        std::cout << "Read a histogram with multiplicity " << mMultiplicity << std::endl;
        mValue.assign(mGridSize * mMultiplicity, 0.0);
        first_time = false;
      }
      vector<double> value;
      for (size_t j = 0; j < mNDim; ++j) {
        pos[j] = std::stod(fields[j]);
      }
      for (size_t j = mNDim; j < fields.size(); ++j) {
        value.push_back(std::stod(fields[j]));
      }
      set(pos, value);
      ++data_count;
      fields.clear();
    }
  }
  if (data_count != mGridSize) {
    std::cerr << "Histogram file reads error!" << std::endl;
  }
}

vector<double>& HistogramVector::getRawData() {
    return mValue;
}

const vector<double>& HistogramVector::getRawData() const {
    const vector<double>& ret = mValue;
    return ret;
}
