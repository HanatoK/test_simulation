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

class BiasExtendedLagrangianBase {
public:
  BiasExtendedLagrangianBase(
    size_t num_extended_cvs,
    const std::vector<double>& tau,
    const std::vector<double>& kappa,
    const std::vector<double>& temperature,
    const std::vector<double>& friction,
    double timestep);
  virtual ~BiasExtendedLagrangianBase() = default;
  void updateExtendedLagrangian();
  double randGaussian();
  void updateCV(const std::vector<double>& position);
  virtual void recordStep(const int64_t& step) {m_step = step;}
  // apply biasing force to the real variable
  void applyBiasForce(std::vector<double>& force);
protected:
  bool                    m_first_time;
  int64_t                 m_step;
  size_t                  m_dof;
  // extended variables
  std::vector<double>     m_mass;
  std::vector<double>     m_forces;
  std::vector<double>     m_bias_force;
  std::vector<double>     m_applied_forces;
  std::vector<double>     m_velocities;
  std::vector<double>     m_positions;
  std::vector<double>     m_kappa;
  std::vector<double>     m_temperature;
  std::vector<double>     m_friction;
  double                  m_timestep;
  std::vector<double>     m_factor1;
  std::vector<double>     m_factor2;
  // real variables
  std::vector<double>     m_real_positions;
  // random generator
  random_device           m_random_device;
  mt19937                 m_random_generator;
  normal_distribution<>   m_normal_distribution;
  // biasing force for the extended variable
  virtual std::vector<double> biasForce(const std::vector<double>& position) = 0;
  virtual void updateForce();
};

class BiasWTMeABF: public BiasExtendedLagrangianBase {
public:
  BiasWTMeABF(const vector<Axis>& ax,
              const vector<Axis>& mtd_ax,
              double tau, double kappa,
              double temperature, double friction, double timestep);
  ~BiasWTMeABF() override;
  // void updateExtendedLagrangian();
  void writeOutput(string filename, size_t freq = 100000000) const;
  void writeTrajectory(std::ostream& os, size_t freq = 100) const;
  void writeHills(std::ostream& os) const;
  double sumHistoryHillsAtPosition(const std::vector<double>& pos) const;
  void recordStep(const int64_t& step) override;
private:
  // ABF + MTD
  HistogramVector<double> m_bias_abf;
  HistogramVector<double> m_bias_mtd;
  HistogramScalar<double> m_mtd_sum_hills;
  HistogramScalar<size_t> m_count;
  Hill                    m_tmp_current_hill;
  int64_t                 m_hill_freq;
  double                  m_hill_initial_height;
  double                  m_bias_temperature;
  std::vector<double>     m_hill_sigma;
  std::vector<double>     m_tmp_grid_pos;
  std::vector<double>     m_tmp_hill_gradient;
  std::vector<double>     m_abf_bias_force;
  std::vector<double>     m_mtd_bias_force;
  double                  m_abf_force_factor;
  double                  m_fullsample;
  // CZAR estimator
  CZARCount               m_zcount;
  HistogramVector<double> m_zgrad;
  // MTD hill traj
  std::vector<Hill>       m_history_hills;
  // override methods
  void updateForce() override;
  std::vector<double> biasForce(const std::vector<double>& position) override;
};

class HarmonicWalls {
public:
  HarmonicWalls(std::vector<double> lower, std::vector<double> upper, double force_constant);
  HarmonicWalls(std::vector<double> lower, std::vector<double> upper, std::vector<double> force_constants);
  void update_value(const std::vector<double>& position);
  const std::vector<double>& position() const {return m_position;}
  double energy() const {return m_energy;}
  const std::vector<double>& force() const {return m_force;}
private:
  std::vector<double> m_lower;
  std::vector<double> m_upper;
  std::vector<double> m_constants;
  std::vector<double> m_position;
  std::vector<double> m_force;
  double m_energy;
};

#endif // BIAS_H
