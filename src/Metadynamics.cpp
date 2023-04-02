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