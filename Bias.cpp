#include "Bias.h"

vector<double> CZARCount::getLogDerivative(const vector<double> &pos) const {
  vector<double> result(mNDim, 0.0);
  size_t addr = 0;
  bool inGrid = address(pos, addr);
  if (!inGrid) return result;
  const double count_this = mValue[addr];
  for (size_t i = 0; i < mNDim; ++i) {
    const double bin_width = mAxes[i].width();
    const size_t addr_first = addr - mAccu[i] * mAxes[i].index(pos[i]) + 0;
    const size_t addr_last = addr_first + mAccu[i] * (mAxes[i].bin() - 1);
    if (addr == addr_first) {
      const double count_next = mValue[addr + mAccu[i]];
      if (mAxes[i].isPeriodic()) {
        const double count_prev = mValue[addr_last];
        if (count_next > 0 && count_prev > 0) {
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2.0 * bin_width);
        }
      } else {
        const double count_next2 = mValue[addr + mAccu[i] * 2];
        if (count_next > 0 && count_next2 > 0 && count_this > 0) {
          result[i] =
            (std::log(count_next2) * (-1.0) + std::log(count_next) * 4.0 -
             std::log(count_this) * 3.0) /
            (2.0 * bin_width);
        }
      }
    } else if (addr == addr_last) {
      const double count_prev = mValue[addr - mAccu[i]];
      if (mAxes[i].isPeriodic()) {
        const double count_next = mValue[addr_first];
        if (count_next > 0 && count_prev > 0) {
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2.0 * bin_width);
        }
      } else {
        const double count_prev2 = mValue[addr - mAccu[i] * 2];
        if (count_prev > 0 && count_this > 0 && count_prev2 > 0) {
          result[i] = (std::log(count_this) * 3.0 - std::log(count_prev) * 4.0 +
                      std::log(count_prev2)) /
                      (2.0 * bin_width);
        }
      }
    } else {
      const double count_prev = mValue[addr - mAccu[i]];
      const double count_next = mValue[addr + mAccu[i]];
      if (count_next > 0 && count_prev > 0)
        result[i] =
          (std::log(count_next) - std::log(count_prev)) / (2 * bin_width);
    }
  }
  return result;
}

BiasWTMeABF2D::BiasWTMeABF2D() {}

BiasWTMeABF2D::BiasWTMeABF2D(
  const vector<Axis>& ax, double fict_mass, double kappa, double temperatue,
  double friction, double timestep):
  m_first_time(true), m_mass{fict_mass},
  m_forces{0, 0, 0}, m_velocities{0, 0, 0}, m_positions{0, 0, 0},
  m_kappa{kappa}, m_temperatue{temperatue}, m_friction{friction},
  m_timestep{timestep}, m_real_positions{0, 0, 0},
  m_random_generator(m_random_device()),
  m_bias_abf(ax, ax.size()), m_bias_mtd(ax, ax.size()),
  m_mtd_sum_hills(ax), m_count(ax), m_zcount(ax), m_zgrad(ax, ax.size()) {
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
  force.x = m_kappa * (m_positions.x - m_real_positions.x);
  force.y = m_kappa * (m_positions.y - m_real_positions.y);
  force.z = m_kappa * (m_positions.z - m_real_positions.z);
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
  // apply bias force
  double3 bias_force = biasForce(m_positions);
  m_forces.x += bias_force.x;
  m_forces.y += bias_force.y;
  m_forces.z += bias_force.z;
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
  const vector<double> tmp_system_f{force.x, force.y};
  m_bias_abf.add(tmp_pos, tmp_system_f);
  m_count.add(tmp_pos, 1.0);
  const Hill h(vector<double>{position.x, position.y},
               vector<double>{0.1, 0.1}, 0.1);
  const double bias_temperature = 3000.0;
  const vector<vector<double>>& point_table = m_bias_mtd.getTable();
  // fixed dimensionality
  const size_t N = 2;
  vector<double> grid_pos(N);
  vector<double> hill_gradient(N);
  // project hills
  vector<double>& sum_hills_array = m_mtd_sum_hills.getRawData();
  vector<double>& mtd_bias_force_array = m_bias_mtd.getRawData();
  for (size_t i = 0; i < point_table[0].size(); ++i) {
    for (size_t j = 0; j < point_table.size(); ++j) {
      grid_pos[j] = point_table[j][i];
    }
    size_t addr = 0;
    m_mtd_sum_hills.address(grid_pos, addr);
    double hill_energy = h.hillEnergy(grid_pos, m_bias_mtd.getAxes());
    hill_gradient = h.hillGradients(grid_pos, m_bias_mtd.getAxes());
    // well-tempered
    const double previous_Vbias = sum_hills_array[addr];
    const double well_tempered_factor = std::exp(-1.0 * previous_Vbias / (boltzmann_constant * bias_temperature));
    hill_energy *= well_tempered_factor;
    for (size_t i = 0; i < N; ++i) {
      mtd_bias_force_array[addr + i] += -1.0 * hill_gradient[i] * well_tempered_factor;
    }
    // Work-in-progress: store sum V and its derivatives
    sum_hills_array[addr] += hill_energy;
  }
  // collect info for CZAR
  const vector<double> tmp_real_pos{m_real_positions.x, m_real_positions.y};
  const vector<double> tmp_real_f{-force.x, -force.y};
  m_zcount.add(tmp_real_pos, 1.0);
  m_zgrad.add(tmp_real_pos, tmp_real_f);
  return force;
}

double3 BiasWTMeABF2D::biasForce(const double3& position) {
  double3 force{0, 0, 0};
  const vector<double> tmp_pos{position.x, position.y};
  vector<double> abf_bias_force(2, 0);
  vector<double> mtd_bias_force(2, 0);
  double count = 0;
  m_bias_abf.get(tmp_pos, abf_bias_force);
  m_bias_mtd.get(tmp_pos, mtd_bias_force);
  m_count.get(tmp_pos, count);
  force.x = -1.0 * (abf_bias_force[0] / count + mtd_bias_force[0]);
  force.y = -1.0 * (abf_bias_force[1] / count + mtd_bias_force[1]);
  return force;
}

void BiasWTMeABF2D::writeOutput(const string& filename) const {
  m_count.writeToFile(filename + ".count");
  m_zcount.writeToFile(filename + ".zcount");
  m_mtd_sum_hills.writeToFile(filename + ".mtd");
  m_bias_abf.writeToFile(filename + ".abf.grad");
  m_zgrad.writeToFile(filename + ".zgrad");
  // write CZAR gradients
  const string czar_grad_filename = filename + ".czar.grad";
  std::ofstream ofs_czar_grad(czar_grad_filename.c_str());
  const vector<vector<double>>& point_table = m_zgrad.getTable();
  const vector<double>& zgrad_data = m_zgrad.getRawData();
  const size_t N = 2;
  vector<double> grid_pos(N);
  for (size_t i = 0; i < point_table[0].size(); ++i) {
    for (size_t j = 0; j < point_table.size(); ++j) {
      grid_pos[j] = point_table[j][i];
      ofs_czar_grad << fmt::format(" {:15.10f}", grid_pos[j]);
    }
    size_t addr = 0;
    m_zgrad.address(grid_pos, addr);
    const vector<double> log_deriv = m_zcount.getLogDerivative(grid_pos);
    // merge gradients
    for (size_t j = 0; j < grid_pos.size(); ++j) {
      const double grad_value = -1.0 / beta() * log_deriv[j] + zgrad_data[addr + j];
      ofs_czar_grad << fmt::format(" {:15.10f}", grad_value);
    }
    ofs_czar_grad << '\n';
  }
}
