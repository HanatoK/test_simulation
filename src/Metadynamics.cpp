#include "Metadynamics.h"
#include "Common.h"

#include <iomanip>
#include <fstream>

Hill::Hill(size_t ndims) {
    init(ndims);
}

void Hill::init(size_t ndims) {
    mCenters.resize(ndims);
    mSigmas.resize(ndims);
    mHeight = 0.0;
}

Hill::Hill(const vector<double>& centers, const vector<double>& sigmas, double height) {
    setParameters(centers, sigmas, height);
}

void Hill::setParameters(const vector<double>& centers, const vector<double>& sigmas, double height) {
    mCenters = centers;
    mSigmas = sigmas;
    mHeight = height;
}

double Hill::hillEnergy(const vector<double>& pos, const vector<Axis>& axes) const {
    double result = 0;
    for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
        const double dist = axes[i_cv].dist(pos[i_cv], mCenters[i_cv]);
        result += dist * dist / (2.0 * mSigmas[i_cv] * mSigmas[i_cv]);
    }
    result = mHeight * std::exp(-result);
    return result;
}

vector<double> Hill::hillGradients(const vector<double>& pos, const vector<Axis>& axes) const {
  vector<double> gradients(pos.size());
  for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
    const double dist = axes[i_cv].dist(mCenters[i_cv], pos[i_cv]);
    const double factor = -1.0 * hillEnergy(pos, axes) / (mSigmas[i_cv] * mSigmas[i_cv]);
    gradients[i_cv] = dist * factor;
  }
  return gradients;
}

void Hill::hillEnergyGradients(
    const vector<double>& pos, const vector<Axis>& axes,
    vector<double>& gradients, double& potential) const {
    double result = 0;
    const size_t N = mCenters.size();
    for (size_t i_cv = 0; i_cv < N; ++i_cv) {
        const double dist = axes[i_cv].dist(pos[i_cv], mCenters[i_cv]);
        result += dist * dist / (2.0 * mSigmas[i_cv] * mSigmas[i_cv]);
        gradients[i_cv] = -dist / (mSigmas[i_cv] * mSigmas[i_cv]);
    }
    if (result > 20.0) {
        potential = 0.0;
        for (size_t i_cv = 0; i_cv < N; ++i_cv) {
            gradients[i_cv] = 0;
        }
    } else {
        potential = mHeight * std::exp(-result);
    }
    for (size_t i_cv = 0; i_cv < N; ++i_cv) {
        gradients[i_cv] *= potential;
    }
}

void Hill::hillGradients(const vector<double>& pos, const vector<Axis>& axes, vector<double>& gradients, bool clear_init_gradients) const {
    if (clear_init_gradients) {
        std::cout << "clear gradients vector!\n";
        gradients.assign(pos.size(), 0);
    }
    for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
        const double dist = axes[i_cv].dist(mCenters[i_cv], pos[i_cv]);
        const double factor = -1.0 * hillEnergy(pos, axes) / (mSigmas[i_cv] * mSigmas[i_cv]);
        gradients[i_cv] += dist * factor;
        if (std::isnan(gradients[i_cv])) {
            std::cout << dist << ' ' << factor << std::endl;
        }
    }
}

void Hill::debugOutput(ostream& os) const {
    std::ios_base::fmtflags f(os.flags());
    os.setf(std::ios::fixed);
    os << std::setprecision(7);
    for (size_t i_cv = 0; i_cv < mCenters.size(); ++i_cv) {
        os << mCenters[i_cv] << "\t ";
    }
    for (size_t i_cv = 0; i_cv < mSigmas.size(); ++i_cv) {
        os << mSigmas[i_cv] << "\t ";
    }
    os << mHeight;
    os << "\n";
    os.flags(f);
}

MetaDynamics::MetaDynamics(const vector<Axis>& ax) : HistogramScalar(ax), mHills(0) {}

void MetaDynamics::addHill(const Hill& hill) {
    mHills.push_back(hill);
}

void MetaDynamics::compute() {
    for (size_t i_hill = 0; i_hill < mHills.size(); ++i_hill) {
        // TODO: Check whether hill is inside the grid
        const Hill& h = mHills[i_hill];
        vector<double> pos(dimension(), 0.0);
        for (size_t i_pos = 0; i_pos < histogramSize(); ++i_pos) {
            for (size_t j_dim = 0; j_dim < dimension(); ++j_dim) {
                pos[j_dim] = mPointTable[j_dim][i_pos];
            }
            const double old_energy = (*this)(pos);
            const double add_energy = h.hillEnergy(pos, mAxes);
            set(pos, old_energy + add_energy);
        }
    }
}

void MetaDynamics::computeEnergyGradients(const vector<double>& pos, double& potential, vector<double>& gradients) const {
    potential = 0;
    gradients.assign(dimension(), 0);
    for (size_t i_hill = 0; i_hill < mHills.size(); ++i_hill) {
        const Hill& h = mHills[i_hill];
        potential += h.hillEnergy(pos, mAxes);
        h.hillGradients(pos, mAxes, gradients, false);
    }
}

void MetaDynamics::readHillsTrajectory(const string& hill_filename) {
    std::ifstream ifs_hilltraj(hill_filename.c_str());
    string line;
    string token{" "};
    vector<string> fields;
    // Parse first line
    std::getline(ifs_hilltraj, line);
    splitString(line, token, fields);
    if (fields[0].compare("#") != 0) {
        std::cerr << "Hills trajectory file reads error!" << std::endl;
        std::abort();
    } else {
        mNdim = std::stoul(fields[1]);
    }
    fields.clear();
    // Parse axes
    mAxes.clear();
    mAxes.resize(mNdim);
    for (size_t i = 0; i < mNdim; ++i) {
        std::getline(ifs_hilltraj, line);
        splitString(line, token, fields);
        double lower, width, upper;
        size_t bins;
        bool periodic = false;
        if (fields[0].compare("#") != 0) {
            std::cerr << "Hills trajectory file reads error!" << std::endl;
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
    mAccu.resize(mNdim);
    mAccu[0] = 1;
    mHistogramSize = 1;
    for (size_t i = 0; i < mNdim; ++i) {
        mAccu[i] = (i == 0) ? 1 : (mAccu[i - 1] * mAxes[i - 1].bin());
        mHistogramSize *= mAxes[i].bin();
    }
    mData.assign(mHistogramSize, 0.0);
    // Initialize table
    fillTable();
    vector<double> pos(mNdim, 0);
    while(std::getline(ifs_hilltraj, line)) {
        splitString(line, token, fields);
        if (fields.empty()) {
            continue;
        }
        if (fields.size() < (2*mNdim + 2)) {
            std::cerr << "Hills trajectory file reads error!" << std::endl;
            std::abort();
        }
        if (fields[0].compare("#") != 0) {
            Hill h(mNdim);
            for (size_t i_cv = 0; i_cv < mNdim; ++i_cv) {
                h.mCenters[i_cv] = std::stod(fields[1+i_cv]);
                h.mSigmas[i_cv] = std::stod(fields[1+1*mNdim+i_cv])/2.0;
            }
            h.mHeight = std::stod(fields[1+2*mNdim]);
            addHill(h);
        }
        fields.clear();
    }
}

void MetaDynamics::writeGradients(const string& filename) const {
    std::ofstream ofs_histo(filename.c_str());
    vector<double> pos(dimension(), 0.0);
    vector<double> grads(dimension(), 0.0);
    double potential = 0;
    ofs_histo << "# " << dimension() << '\n';
    for (size_t j = 0; j < dimension(); ++j) {
        ofs_histo << mAxes[j].infoHeader() << '\n';
    }
    ofs_histo.setf(std::ios::fixed);
    ofs_histo << std::setprecision(7);
    for (size_t i = 0; i < histogramSize(); ++i) {
        for (size_t j = 0; j < dimension(); ++j) {
            pos[j] = mPointTable[j][i];
            ofs_histo << pos[j] << ' ';
        }
        computeEnergyGradients(pos, potential, grads);
        for (size_t j = 0; j < dimension(); ++j) {
            ofs_histo << grads[j] << ' ';
        }
        ofs_histo << '\n';
    }
}

void MetaDynamics::writeDummyCounts(const string& filename) const {
    std::ofstream ofs_histo(filename.c_str());
    vector<double> pos(dimension(), 0.0);
    ofs_histo << "# " << dimension() << '\n';
    for (size_t j = 0; j < dimension(); ++j) {
        ofs_histo << mAxes[j].infoHeader() << '\n';
    }
    ofs_histo.setf(std::ios::fixed);
    ofs_histo << std::setprecision(7);
    for (size_t i = 0; i < histogramSize(); ++i) {
        for (size_t j = 0; j < dimension(); ++j) {
            pos[j] = mPointTable[j][i];
            ofs_histo << pos[j] << ' ';
        }
        ofs_histo << 1 << ' ';
        ofs_histo << '\n';
    }
}
