#ifndef BIAS_H
#define BIAS_H

#include "Histogram.h"
#include "Common.h"
#include "Metadynamics.h"
#include <random>

using std::random_device;
using std::mt19937;
using std::normal_distribution;

class CZARCount: public virtual HistogramScalar<size_t> {
public:
  CZARCount() {}
  CZARCount(const vector<Axis>& ax): HistogramBase(ax), HistogramScalar(ax) {
  }
  virtual ~CZARCount() {}
  virtual vector<double> getLogDerivative(const vector<double> &pos) const;
};

// TODO: new hill frequency
class BiasWTMeABF2D {
public:
  BiasWTMeABF2D();
  BiasWTMeABF2D(const vector<Axis>& ax,
                const vector<Axis>& mtd_ax,
                double tau, double kappa,
                double temperatue, double friction, double timestep,
                const std::string& hill_traj_filename);
  ~BiasWTMeABF2D();
  void updateExtendedLagrangian();
  void applyBiasForce(double3& force);
  double randGaussian();
  double beta() const;
  void writeOutput(const string& filename) const;
  void recordStep(const int64_t& step) {m_step = step;}
  void updateCV(const double3& position);
  void writeTrajectory(std::ostream& os) const;
  double sumHistoryHillsAtPosition(const std::vector<double>& pos) const;
private:
  bool                    m_first_time;
  int64_t                 m_step;
  // extended variables
  double                  m_mass;
  double3                 m_forces;
  double3                 m_applied_forces;
  double3                 m_velocities;
  double3                 m_positions;
  double                  m_kappa;
  double                  m_temperatue;
  double                  m_friction;
  double                  m_timestep;
  double                  m_factor1;
  double                  m_factor2;
  // real variables
  double3                 m_real_positions;
  // random generator
  random_device           m_random_device;
  mt19937                 m_random_generator;
  normal_distribution<>   m_normal_distribution;
  // ABF + MTD
  HistogramVector<double> m_bias_abf;
  HistogramVector<double> m_bias_mtd;
  HistogramScalar<double> m_mtd_sum_hills;
  HistogramScalar<size_t> m_count;
  double3                 m_bias_force;
  Hill                    m_tmp_current_hill;
  int64_t                 m_hill_freq;
  double                  m_hill_initial_height;
  std::vector<double>     m_hill_sigma;
  std::vector<double>     m_tmp_grid_pos;
  std::vector<double>     m_tmp_hill_gradient;
  std::vector<double>     m_tmp_pos;
  std::vector<double>     m_tmp_system_f;
  // CZAR estimator
  CZARCount               m_zcount;
  HistogramVector<double> m_zgrad;
  void updateForce();
  double3 biasForce(const double3& position);
  // MTD hill traj
  std::ofstream           m_hill_traj;
  std::vector<Hill>       m_history_hills;
  // restraint
  double3 restraintForce(const double3& position);
};

#endif // BIAS_H
