#ifndef BIAS_H
#define BIAS_H

#include "Grid.h"
#include "Simulation.h"
#include "MetaDynamics.h"

using std::random_device;
using std::mt19937;
using std::normal_distribution;

class BiasWTMeABF2D {
public:
  BiasWTMeABF2D();
  BiasWTMeABF2D(const vector<Axis>& ax);
  ~BiasWTMeABF2D();
  void positionCallback(double3& position);
  void applyBiasForce(double3& force);
  double randGaussian();
  double beta() const;
private:
  bool                    m_first_time;
  // extended variables
  double                  m_mass;
  double3                 m_forces;
  double3                 m_velocities;
  double3                 m_positions;
  double                  m_temperatue;
  double                  m_kappa;
  double                  m_friction;
  double                  m_tau;
  double                  m_factor1;
  double                  m_factor2;
  // real variables
  double3                 m_real_positions;
  // random generator
  random_device           m_random_device;
  mt19937                 m_random_generator;
  normal_distribution<>   m_normal_distribution;
  // ABF + MTD
  HistogramVector         m_bias_abf;
  HistogramVector         m_bias_mtd;
  HistogramScalar         m_mtd_sum_hills;
  HistogramScalar         m_count;
  double3 updateForce(const double3& position);
  double3 biasForce(const double3& position);
};

#endif // BIAS_H
