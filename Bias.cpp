#include "Bias.h"

BiasWTMeABF2D::BiasWTMeABF2D() {}

BiasWTMeABF2D::BiasWTMeABF2D(
  const vector<Axis>& ax, double fict_mass, double kappa, double temperatue,
  double friction, double timestep):
  m_first_time(true), m_mass{fict_mass},
  m_forces{0, 0, 0}, m_velocities{0, 0, 0}, m_positions{0, 0, 0},
  m_kappa{kappa}, m_temperatue{temperatue}, m_friction{friction},
  m_timestep{timestep}, m_real_positions{0, 0, 0},
  m_random_generator(m_random_device),
  m_bias_abf(ax, ax.size()), m_bias_mtd(ax, ax.size()),
  m_mtd_sum_hills(ax), m_count(ax) {
  m_factor1 = std::exp(-1.0 * m_friction * m_timestep);
  m_factor2 = std::sqrt(1.0 / (beta() * m_mass)) *
              std::sqrt(1.0 - std::exp(-2.0 * m_friction * m_timestep));
}

BiasWTMeABF2D::~BiasWTMeABF2D() {}

double BiasWTMeABF2D::randGaussian() {
  return m_normal_distribution(m_random_generator);
}

double BiasWTMeABF2D::beta() const {
  return 1.0 / (m_temperatue * boltzmann_constant);
}

void BiasWTMeABF2D::applyBiasForce(double3& force) {
  force.x = m_kappa * (position.x - m_real_positions.x);
  force.y = m_kappa * (position.y - m_real_positions.y);
  force.z = m_kappa * (position.z - m_real_positions.z);
}

void BiasWTMeABF2D::positionCallback(double3& position) {
  m_real_positions = position;
  if (m_first_time) {
    m_positions = m_real_positions;
    m_velocities.x = randGaussian();
    m_velocities.y = randGaussian();
    m_velocities.z = randGaussian();
    const double alpha = std::sqrt(1 / (beta() * m_mass));
    m_velocities.x *= alpha;
    m_velocities.y *= alpha;
    m_velocities.z *= alpha;
    m_first_time = false;
  }
  // update v_{i+1/2}
  m_velocities.x += 0.5 * m_timestep * m_forces.x / m_mass;
  m_velocities.y += 0.5 * m_timestep * m_forces.y / m_mass;
  m_velocities.z += 0.5 * m_timestep * m_forces.z / m_mass;
  // update x_{i+1/2}
  m_positions.x += 0.5 * m_velocities.x * m_timestep; 
  m_positions.y += 0.5 * m_velocities.y * m_timestep;
  m_positions.z += 0.5 * m_velocities.z * m_timestep;
  // Langevin thermostat, full step
  m_velocities.x = m_factor1 * m_velocities.x + m_factor2 * randGaussian();
  m_velocities.y = m_factor1 * m_velocities.y + m_factor2 * randGaussian();
  m_velocities.z = m_factor1 * m_velocities.z + m_factor2 * randGaussian();
  // update x_{i+1}
  m_positions.x += 0.5 * m_velocities.x * m_timestep;
  m_positions.y += 0.5 * m_velocities.y * m_timestep;
  m_positions.z += 0.5 * m_velocities.z * m_timestep;
  // update f_{i+1}
  m_forces = updateForce(m_positions);
  m_forces += biasForce(m_positions);
  // update v_{i+1}
  m_velocities.x += 0.5 * m_timestep * m_forces.x / m_mass;
  m_velocities.y += 0.5 * m_timestep * m_forces.y / m_mass;
  m_velocities.z += 0.5 * m_timestep * m_forces.z / m_mass;
}

double3 BiasWTMeABF2D::updateForce(const double3& position) {
  double3 force{0, 0, 0};
  force.x = -1.0 * m_kappa * (position.x - m_real_positions.x);
  force.y = -1.0 * m_kappa * (position.y - m_real_positions.y);
  force.z = -1.0 * m_kappa * (position.z - m_real_positions.z);
  // store the instantaneous force to ABF
  const vector<double> tmp_pos{position.x, position.y};
  const vector<double> tmp_f{force.x, force.y};
  m_bias_abf.add(tmp_pos, tmp_f);
  m_count.add(tmp_pos, 1.0);
  // TODO: update MTD bias force grid
  // FIXME: Fixed MTD parameters are used!
  const Hill h(vector<double>{position.x, position.y},
               vector<double>{0.1, 0.1}, 0.1);
  const double bias_temperature = 3000.0;
  const vector<vector<double>>& point_table = m_bias_mtd.getTable();
  vector<double> grid_pos(2);
  vector<double> hill_gradient(2);
  for (size_t i = 0; i < point_table[0].size(); ++i) {
    for (size_t j = 0; j < point_table.size(); ++j) {
      grid_pos[j] = point_table[j][i];
    }
    double hill_energy = h.hillEnergy(grid_pos, m_bias_mtd.getAxes());
    hill_gradient = h.hillGradients(grid_pos, m_bias_mtd.getAxes());
    // well-tempered
    const double previous_Vbias = 0;
    m_mtd_sum_hills.get(grid_pos, previous_Vbias);
    const double well_tempered_factor = std::exp(-1.0 * previous_Vbias / (boltzmann_constant * bias_temperature));
    hill_energy *= well_tempered_factor;
    for (auto& grad: hill_gradient) {
      grad *= well_tempered_factor;
    }
    // Work-in-progress: store sum V and its derivatives
  }
  return force;
}

double3 BiasWTMeABF2D::biasForce(const double3& position) {
  double3 force{0, 0, 0};
  const vector<double> tmp_pos{position.x, position.y};
  vector<double> abf_bias_force(2, 0);
  vector<double> mtd_bias_force(2, 0);
  m_bias_abf.get(tmp_pos, abf_bias_force);
  m_bias_mtd.get(tmp_pos, mtd_bias_force);
  force.x = -1.0 * (abf_bias_force[0] + mtd_bias_force[0]);
  force.y = -1.0 * (abf_bias_force[1] + mtd_bias_force[1]);
  return force;
}
