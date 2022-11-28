#include "Bias.h"
#include <fmt/format.h>
#include <numbers>

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
      if (mAxes[i].isPeriodic()) {
        const double count_next = mValue[addr + mAccu[i]];
        const double count_prev = mValue[addr_last];
        if (count_next > 0 && count_prev > 0) {
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2.0 * bin_width);
        }
      } else {
        const double count_next = mValue[addr + mAccu[i]];
        const double count_next2 = mValue[addr + mAccu[i] * 2];
        if (count_next > 0 && count_next2 > 0 && count_this > 0) {
          result[i] =
            (std::log(count_next2) * (-1.0) + std::log(count_next) * 4.0 -
             std::log(count_this) * 3.0) /
            (2.0 * bin_width);
        }
      }
    } else if (addr == addr_last) {
      if (mAxes[i].isPeriodic()) {
        const double count_prev = mValue[addr - mAccu[i]];
        const double count_next = mValue[addr_first];
        if (count_next > 0 && count_prev > 0) {
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2.0 * bin_width);
        }
      } else {
        const double count_prev = mValue[addr - mAccu[i]];
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
  const vector<Axis>& ax, const vector<Axis>& mtd_ax,
  double tau,
  double kappa, double temperatue,
  double friction, double timestep,
  const std::string& hill_traj_filename):
  m_first_time(true), m_step(0), m_mass{kappa * tau * tau / (4.0 * std::numbers::pi * std::numbers::pi)},
  m_forces{0, 0, 0}, m_velocities{0, 0, 0}, m_positions{0, 0, 0},
  m_kappa{kappa}, m_temperatue{temperatue}, m_friction{friction},
  m_timestep{timestep}, m_real_positions{0, 0, 0},
  m_random_generator(m_random_device()),
  m_bias_abf(ax, ax.size()), m_bias_mtd(mtd_ax, mtd_ax.size()),
  m_mtd_sum_hills(mtd_ax), m_count(ax), m_bias_force{0, 0, 0},
  m_tmp_current_hill(2), m_hill_freq(1000),
  m_hill_initial_height(0.1), m_hill_sigma(2), m_tmp_grid_pos(2),
  m_tmp_hill_gradient(2), m_tmp_pos(2), m_tmp_system_f(2),
  m_zcount(ax), m_zgrad(ax, ax.size()) {
  m_factor1 = std::exp(-1.0 * m_friction * m_timestep);
  m_factor2 = std::sqrt(1.0 / (beta() * m_mass)) *
              std::sqrt(1.0 - std::exp(-2.0 * m_friction * m_timestep));
  m_hill_sigma[0] = 8.0 * 0.05;
  m_hill_sigma[1] = 8.0 * 0.05;
  m_hill_traj.open(hill_traj_filename);
  m_hill_traj << "# step x y sigma_x sigma_y height\n";
}

BiasWTMeABF2D::~BiasWTMeABF2D() {
  m_hill_traj.close();
}

double BiasWTMeABF2D::randGaussian() {
  return m_normal_distribution(m_random_generator);
}

double BiasWTMeABF2D::beta() const {
  return 1.0 / (m_temperatue * boltzmann_constant);
}

void BiasWTMeABF2D::applyBiasForce(double3& force) {
  // apply the biasing force to real position
  force.x += m_kappa * (m_positions.x - m_real_positions.x);
  force.y += m_kappa * (m_positions.y - m_real_positions.y);
  force.z += m_kappa * (m_positions.z - m_real_positions.z);
  m_forces = updateForce(m_positions);
}

void BiasWTMeABF2D::positionCallback(const double3& position) {
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
    m_bias_force = biasForce(m_positions);
    return;
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
  // m_forces = updateForce(m_positions);
  m_bias_force = biasForce(m_positions);
  const auto restraint_force = restraintForce(m_positions);
  // apply bias force
  m_forces.x += m_bias_force.x;
  m_forces.y += m_bias_force.y;
  m_forces.z += m_bias_force.z;
  m_forces.x += restraint_force.x;
  m_forces.y += restraint_force.y;
  m_forces.z += restraint_force.z;
  // update v_{i+1}
  m_velocities.x += 0.5 * m_timestep * m_forces.x / m_mass;
  m_velocities.y += 0.5 * m_timestep * m_forces.y / m_mass;
  m_velocities.z += 0.5 * m_timestep * m_forces.z / m_mass;
}

void BiasWTMeABF2D::positionCallback2(const double3& position) {
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
    return;
  }
  const double c1 = std::sqrt(m_factor1);
  const double c2 = std::sqrt((1.0-c1*c1)/(beta() * m_mass));
  // update f_{i+1}
  // m_forces = updateForce(m_positions);
  m_bias_force = biasForce(m_positions);
  const auto restraint_force = restraintForce(m_positions);
  // apply bias force
  m_forces.x += m_bias_force.x;
  m_forces.y += m_bias_force.y;
  m_forces.z += m_bias_force.z;
  m_forces.x += restraint_force.x;
  m_forces.y += restraint_force.y;
  m_forces.z += restraint_force.z;
  // update v_{i+1/2}
  m_velocities.x += 0.5 * m_timestep * m_forces.x / m_mass;
  m_velocities.y += 0.5 * m_timestep * m_forces.y / m_mass;
  m_velocities.z += 0.5 * m_timestep * m_forces.z / m_mass;
  // Langevin thermostat, half step
  m_velocities.x = c1 * m_velocities.x + c2 * randGaussian();
  m_velocities.y = c1 * m_velocities.y + c2 * randGaussian();
  m_velocities.z = c1 * m_velocities.z + c2 * randGaussian();
  // Langevin thermostat, half step
  m_velocities.x = c1 * m_velocities.x + c2 * randGaussian();
  m_velocities.y = c1 * m_velocities.y + c2 * randGaussian();
  m_velocities.z = c1 * m_velocities.z + c2 * randGaussian();
  // update v_{i+1}
  m_velocities.x += 0.5 * m_timestep * m_forces.x / m_mass;
  m_velocities.y += 0.5 * m_timestep * m_forces.y / m_mass;
  m_velocities.z += 0.5 * m_timestep * m_forces.z / m_mass;
  // update x_{i+1}
  m_positions.x += m_velocities.x * m_timestep;
  m_positions.y += m_velocities.y * m_timestep;
  m_positions.z += m_velocities.z * m_timestep;
}

double BiasWTMeABF2D::sumHistoryHillsAtPosition(const std::vector<double>& pos) const {
  double potential = 0;
  const auto& ax = m_bias_mtd.getAxes();
  for (size_t i = 0; i < m_history_hills.size(); ++i) {
    potential += m_history_hills[i].hillEnergy(pos, ax);
  }
  return potential;
}

// TODO: tune MTD parameters
double3 BiasWTMeABF2D::updateForce(const double3& position) {
  // std::cerr << "Update force at step " << m_step << std::endl;
  // ABF
  // system force acting on position
  const double3 system_force{
    -1.0 * m_kappa * (position.x - m_real_positions.x),
    -1.0 * m_kappa * (position.y - m_real_positions.y),
    -1.0 * m_kappa * (position.z - m_real_positions.z)
  };
  m_tmp_pos[0] = position.x;
  m_tmp_pos[1] = position.y;
  // m_tmp_system_f[0] = system_force.x;
  // m_tmp_system_f[1] = system_force.y;
  // compute the negative average force and store it
  double current_count = 0;
  if (m_count.get(m_tmp_pos, current_count)) {
    std::vector<double> previous_system_force{0, 0};
    m_bias_abf.get(m_tmp_pos, previous_system_force);
    m_tmp_system_f[0] = (-system_force.x / (current_count + 1.0)) + previous_system_force[0] * (current_count / (current_count + 1.0));
    m_tmp_system_f[1] = (-system_force.y / (current_count + 1.0)) + previous_system_force[1] * (current_count / (current_count + 1.0));
    m_bias_abf.set(m_tmp_pos, m_tmp_system_f);
    m_count.add(m_tmp_pos, 1.0);
    // std::cout << fmt::format("system force: {:15.10f} {:15.10f}, previous avg: {:15.10f} {:15.10f}, count: {:15.10f}\n",
    //                          system_force.x, system_force.y, previous_system_force[0], previous_system_force[1], current_count);
  }
  // collect info for CZAR
  const vector<double> tmp_real_pos{m_real_positions.x, m_real_positions.y};
  double current_zcount = 0.0;
  if (m_zcount.get(tmp_real_pos, current_zcount)) {
    std::vector<double> previous_zgrad{0, 0};
    std::vector<double> tmp_real_f{0, 0};
    m_zgrad.get(tmp_real_pos, previous_zgrad);
    tmp_real_f[0] = (system_force.x / (current_zcount + 1.0)) + previous_zgrad[0] * (current_zcount / (current_zcount + 1.0));
    tmp_real_f[1] = (system_force.y / (current_zcount + 1.0)) + previous_zgrad[1] * (current_zcount / (current_zcount + 1.0));
    m_zgrad.set(tmp_real_pos, tmp_real_f);
    m_zcount.add(tmp_real_pos, 1.0);
  }
  // MTD
  // if (m_step % m_hill_freq == 0) {
  //   // setup hill
  //   m_tmp_current_hill.mCenters[0] = position.x;
  //   m_tmp_current_hill.mCenters[1] = position.y;
  //   m_tmp_current_hill.mSigmas[0] = m_hill_sigma[0];
  //   m_tmp_current_hill.mSigmas[1] = m_hill_sigma[1];
  //   // well-tempered MTD requires the previous biasing potential
  //   double previous_bias_V = 0;
  //   vector<double>& sum_hills_array = m_mtd_sum_hills.getRawData();
  //   if (!m_mtd_sum_hills.isInGrid(m_tmp_current_hill.mCenters)) {
  //     // if the position is outside the boundary, then sum it from the history
  //     std::cerr << fmt::format("Warning: at step {:d}, a hill is deposited out-of-bound at {:.10f} {:.10f}\n",
  //                              m_step, m_tmp_current_hill.mCenters[0], m_tmp_current_hill.mCenters[1]);
  //     previous_bias_V = sumHistoryHillsAtPosition(m_tmp_current_hill.mCenters);
  //   } else {
  //     // fast path: if the hill is inside the boundary, then get it from the grid
  //     size_t addr = 0;
  //     m_mtd_sum_hills.address(m_tmp_current_hill.mCenters, addr);
  //     previous_bias_V = sum_hills_array[addr];
  //   }
  //   const double bias_temperature = 3000.0;
  //   const double well_tempered_factor = std::exp(-1.0 * previous_bias_V / (boltzmann_constant * bias_temperature));
  //   m_tmp_current_hill.mHeight = well_tempered_factor * m_hill_initial_height;
  //   // save the current hill
  //   m_history_hills.push_back(m_tmp_current_hill);
  //   // print hill trajectory
  //   m_hill_traj << fmt::format(" {:>15d} {:15.10f} {:15.10f} {:15.10f} {:15.10f} {:15.10f}\n",
  //                              m_step, m_tmp_current_hill.mCenters[0], m_tmp_current_hill.mCenters[1],
  //                              m_tmp_current_hill.mSigmas[0], m_tmp_current_hill.mSigmas[1],
  //                              m_tmp_current_hill.mHeight);
  //   // project hill
  //   const vector<vector<double>>& point_table = m_bias_mtd.getTable();
  //   vector<double>& mtd_bias_force_array = m_bias_mtd.getRawData();
  //   // TODO: should be parallelized for performance
  //   // std::cout << "dimension: " << point_table.size() << std::endl;
  //   for (size_t i = 0; i < point_table[0].size(); ++i) {
  //     for (size_t j = 0; j < point_table.size(); ++j) {
  //       m_tmp_grid_pos[j] = point_table[j][i];
  //     }
  //     // m_tmp_grid_pos: the grid point in the histogram
  //     size_t addr = 0;
  //     double hill_energy = 0;
  //     m_mtd_sum_hills.address(m_tmp_grid_pos, addr);
  //     // compute hill energy and gradients
  //     m_tmp_current_hill.hillEnergyGradients(
  //       m_tmp_grid_pos, m_bias_mtd.getAxes(),
  //       m_tmp_hill_gradient, hill_energy);
  //     for (size_t i = 0; i < point_table.size(); ++i) {
  //       mtd_bias_force_array[addr + i] += -1.0 * m_tmp_hill_gradient[i];
  //     }
  //     sum_hills_array[addr] += hill_energy;
  //   }
  // }
  // std::cerr << "Update force done" << std::endl;
  return system_force;
}

double3 BiasWTMeABF2D::biasForce(const double3& position) {
  double3 force{0, 0, 0};
  const vector<double> tmp_pos{position.x, position.y};
  vector<double> abf_bias_force(2, 0);
  vector<double> mtd_bias_force(2, 0);
  double count = 0;
  m_bias_abf.get(tmp_pos, abf_bias_force);
  // m_bias_mtd.get(tmp_pos, mtd_bias_force);
  m_count.get(tmp_pos, count);
  // abf_bias_force is actually the sum of instantaneous collective force
  const double fullsample = 200.0;
  double abf_force_factor = 0;
  if (count < fullsample / 2) {
    abf_force_factor = 0;
  } else if (count > fullsample) {
    abf_force_factor = 1;
  } else {
    abf_force_factor = count / fullsample;
  }
  force.x = abf_force_factor * abf_bias_force[0] + mtd_bias_force[0];
  force.y = abf_force_factor * abf_bias_force[1] + mtd_bias_force[1];
  if (m_step % 1000 == 0) {
    std::cout << fmt::format("step {:10d}, abf force: {:15.10f} {:15.10f}, "
                             "factor: {:12.7f} "
                             "mtd force: {:15.10f} {:15.10f}\n",
                             m_step, abf_force_factor * abf_bias_force[0], abf_force_factor * abf_bias_force[1],
                             abf_force_factor,
                             mtd_bias_force[0], mtd_bias_force[1]);
  }
  return force;
}

double3 BiasWTMeABF2D::restraintForce(const double3& position) {
  // wall boundaries at -6 and 6
  double3 force{0, 0, 0};
  const double force_constant = 10000.0;
  const double x_lower = -4.0;
  const double x_upper = 4.0;
  const double y_lower = -4.0;
  const double y_upper = 4.0;
  if (position.x < x_lower) {
    force.x = -1.0 * force_constant * (position.x - x_lower);
  } else if (position.x > x_upper) {
    force.x = -1.0 * force_constant * (position.x - x_upper);
  }
  if (position.y < y_lower) {
    force.y = -1.0 * force_constant * (position.y - y_lower);
  } else if (position.y > y_upper) {
    force.y = -1.0 * force_constant * (position.y - y_upper);
  }
  return force;
}

void BiasWTMeABF2D::writeOutput(const string& filename) const {
  m_count.writeToFile(filename + ".abf.count");
  m_zcount.writeToFile(filename + ".zcount");
  m_mtd_sum_hills.writeToFile(filename + ".mtd");
  m_bias_abf.writeToFile(filename + ".abf.grad");
  m_zgrad.writeToFile(filename + ".zgrad");
  // write CZAR gradients
  const string czar_grad_filename = filename + ".czar.grad";
  std::ofstream ofs_czar_grad(czar_grad_filename.c_str());
  ofs_czar_grad << "# " << m_zgrad.getDimension() << '\n';
  for (size_t j = 0; j < m_zgrad.getDimension(); ++j) {
    ofs_czar_grad << m_zgrad.getAxes()[j].infoHeader() << '\n';
  }
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
      const double grad_value = -1.0 / beta() * log_deriv[j] + zgrad_data[addr * grid_pos.size() + j];
      ofs_czar_grad << fmt::format(" {:15.10f}", grad_value);
    }
    ofs_czar_grad << '\n';
  }
}

void BiasWTMeABF2D::writeTrajectory(std::ostream& os) const {
  static bool header = true;
  if (header) {
    os << "# step x r_x y r_y fb_rx fb_ry\n";
    header = false;
  }
  os << fmt::format(" {:>15d} {:15.10f} {:15.10f} {:15.10f} {:15.10f}"
                    " {:15.10f} {:15.10f}\n", m_step,
                    m_real_positions.x, m_positions.x, 
                    m_real_positions.y, m_positions.y,
                    m_bias_force.x, m_bias_force.y);
}
