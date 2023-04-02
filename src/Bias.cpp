#include "Bias.h"
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <numbers>

vector<double> CZARCount::getLogDerivative(const vector<double> &pos) const {
  vector<double> result(dimension(), 0.0);
  bool inGrid = true;
  size_t addr = address(pos, &inGrid);
  if (!inGrid) return result;
  const double count_this = mData[addr];
  for (size_t i = 0; i < dimension(); ++i) {
    const double bin_width = mAxes[i].width();
    const size_t addr_first = addr - mAccu[i] * mAxes[i].index(pos[i]) + 0;
    const size_t addr_last = addr_first + mAccu[i] * (mAxes[i].bin() - 1);
    if (addr == addr_first) {
      if (mAxes[i].periodic()) {
        const double count_next = mData[addr + mAccu[i]];
        const double count_prev = mData[addr_last];
        if (count_next > 0 && count_prev > 0) {
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2.0 * bin_width);
        }
      } else {
        const double count_next = mData[addr + mAccu[i]];
        const double count_next2 = mData[addr + mAccu[i] * 2];
        if (count_next > 0 && count_next2 > 0 && count_this > 0) {
          result[i] =
            (std::log(count_next2) * (-1.0) + std::log(count_next) * 4.0 -
             std::log(count_this) * 3.0) /
            (2.0 * bin_width);
        }
      }
    } else if (addr == addr_last) {
      if (mAxes[i].periodic()) {
        const double count_prev = mData[addr - mAccu[i]];
        const double count_next = mData[addr_first];
        if (count_next > 0 && count_prev > 0) {
          result[i] =
            (std::log(count_next) - std::log(count_prev)) / (2.0 * bin_width);
        }
      } else {
        const double count_prev = mData[addr - mAccu[i]];
        const double count_prev2 = mData[addr - mAccu[i] * 2];
        if (count_prev > 0 && count_this > 0 && count_prev2 > 0) {
          result[i] = (std::log(count_this) * 3.0 - std::log(count_prev) * 4.0 +
                      std::log(count_prev2)) /
                      (2.0 * bin_width);
        }
      }
    } else {
      const double count_prev = mData[addr - mAccu[i]];
      const double count_next = mData[addr + mAccu[i]];
      if (count_next > 0 && count_prev > 0)
        result[i] =
          (std::log(count_next) - std::log(count_prev)) / (2 * bin_width);
    }
  }
  return result;
}

BiasExtendedLagrangianBase::BiasExtendedLagrangianBase(
  size_t num_extended_cvs,
  const std::vector<double>& tau,
  const std::vector<double>& kappa,
  const std::vector<double>& temperature,
  const std::vector<double>& friction,
  double timestep
): m_first_time(true), m_step(0), m_dof(num_extended_cvs),
   m_mass(m_dof), m_forces(m_dof, 0), m_bias_force(m_dof, 0),
   m_applied_forces(m_dof, 0),
   m_velocities(m_dof, 0), m_positions(m_dof),
   m_kappa(kappa), m_temperature(temperature), m_friction(friction),
   m_timestep{timestep}, m_factor1(m_dof),
   m_factor2(m_dof), m_real_positions(m_dof),
   m_random_generator(m_random_device()),
   m_half_timestep(m_dof), m_half_timestep_divide_mass(m_dof) {
  std::cout << "Initialize extended Lagrangian:\n";
  for (size_t i = 0; i < m_dof; ++i) {
    m_mass[i] = kappa[i] * tau[i] * tau[i]  / (4.0 * std::numbers::pi * std::numbers::pi);
    m_factor1[i] = std::exp(-1.0 * m_friction[i] * m_timestep);
    m_factor2[i] = std::sqrt(1.0 / (beta(m_temperature[i]) * m_mass[i])) *
                   std::sqrt(1.0 - std::exp(-2.0 * m_friction[i] * m_timestep));
    std::cout << fmt::format("  variable {}: kappa = {:10.5f}, mass = {:10.5f}, temperature = {:10.5f}, friction = {:10.5f}\n",
                             i, m_kappa[i], m_mass[i], m_temperature[i], m_friction[i]);
  }
  for (size_t i = 0; i < m_dof; ++i) {
    m_half_timestep[i] = 0.5 * m_timestep;
    m_half_timestep_divide_mass[i] = 0.5 * m_timestep / m_mass[i];
  }
}

double BiasExtendedLagrangianBase::randGaussian() {
  return m_normal_distribution(m_random_generator);
}

void BiasExtendedLagrangianBase::updateExtendedLagrangian() {
  // BAOAB
  for (size_t i = 0; i < m_dof; ++i) {
    // update v_{i+1/2}
    m_velocities[i] += m_half_timestep_divide_mass[i] * m_forces[i];
    // update x_{i+1/2}
    m_positions[i] += m_half_timestep[i] * m_velocities[i];
    // Langevin thermostat, full step
    m_velocities[i] = m_factor1[i] * m_velocities[i] + m_factor2[i] * randGaussian();
    // update x_{i+1}
    m_positions[i] += m_half_timestep[i] * m_velocities[i];
  }
  // update f_{i+1}
  updateForce();
  // report potential energy here
  m_bias_force = biasForce(m_positions);
  for (size_t i = 0; i < m_dof; ++i) {
    // update f_{i+1}
    m_forces[i] += m_bias_force[i];
    // update v_{i+1}
    m_velocities[i] += m_half_timestep_divide_mass[i] * m_forces[i];
  }
  // report kinetic energy here
}

void BiasExtendedLagrangianBase::applyBiasForce(std::vector<double>& force) {
  // apply the biasing force to real position
  for (size_t i = 0; i < m_dof; ++i) {
    force[i] += m_applied_forces[i];
  }
}

void BiasExtendedLagrangianBase::updateCV(const std::vector<double>& position) {
  m_real_positions = position;
  if (m_first_time) {
    m_positions = m_real_positions;
    for (size_t i = 0; i < m_dof; ++i) {
      m_velocities[i] = randGaussian() * std::sqrt(1 / (beta(m_temperature[i]) * m_mass[i]));
    }
    m_first_time = false;
  }
  // updateForce();
}

void BiasExtendedLagrangianBase::updateForce() {
  for (size_t i = 0; i < m_dof; ++i) {
    // system force acting on position
    m_forces[i] = -1.0 * m_kappa[i] * (m_positions[i] - m_real_positions[i]);
    // force applied back to real positions
    m_applied_forces[i] = -m_forces[i];
  }
}


BiasWTMeABF::BiasWTMeABF(
  const vector<Axis>& ax, const vector<Axis>& mtd_ax,
  double tau,
  double kappa, double temperatue,
  double friction, double timestep,
  double hill_initial_height, int64_t new_hill_frequency,
  double bias_temperature, double sigma_width_ratio,
  double fullsample):
  BiasExtendedLagrangianBase(
    ax.size(), std::vector(ax.size(), tau), std::vector(ax.size(), kappa),
    std::vector(ax.size(), temperatue), std::vector(ax.size(), friction), timestep),
  m_bias_abf(ax, ax.size()), m_bias_mtd(mtd_ax, mtd_ax.size()),
  m_mtd_sum_hills(mtd_ax), m_count(ax),
  m_tmp_current_hill(ax.size()), m_hill_freq(new_hill_frequency),
  m_hill_initial_height(hill_initial_height), m_bias_temperature(bias_temperature),
  m_hill_sigma(ax.size()),
  m_tmp_hill_gradient(ax.size()), m_abf_bias_force(ax.size(), 0),
  m_mtd_bias_force(ax.size(), 0), m_abf_force_factor(0),
  m_fullsample(200.0), m_zcount(ax), m_zgrad(ax, ax.size()) {
  std::cout << "Initialize WTM-eABF calculation:\n";
  std::cout << "  ABF fullsamples: " << m_fullsample << "\n";
  std::cout << "  ABF grid settings on extended Lagrangian:\n";
  for (size_t i = 0; i < m_dof; ++i) {
    std::cout << fmt::format("    Axis {}: lower = {:10.5f} ; upper = {:10.5f} ; width = {:10.5f} ; bins = {}\n",
                             i, m_bias_abf.axes()[i].lowerBound(), m_bias_abf.axes()[i].upperBound(),
                             m_bias_abf.axes()[i].width(), m_bias_abf.axes()[i].bin());
  }
  std::cout << "  Well-tempered metadynamics settings on extended Lagrangian:\n";
  std::cout << "    Hill initial height: " << m_hill_initial_height << "\n";
  std::cout << "    New hill frequency: " << m_hill_freq << "\n";
  std::cout << "    Bias temperature: " << m_bias_temperature << "\n";
  for (size_t i = 0; i < m_dof; ++i) {
    m_hill_sigma[i] = sigma_width_ratio * ax[i].width();
    std::cout << fmt::format("    Hill sigma at variable {}: {}\n", i, m_hill_sigma[i]);
  }
  std::cout << "  Well-tempered grid:\n";
  for (size_t i = 0; i < m_dof; ++i) {
    std::cout << fmt::format("    Axis {}: lower = {:10.5f} ; upper = {:10.5f} ; width = {:10.5f} ; bins = {}\n",
                             i, m_bias_mtd.axes()[i].lowerBound(), m_bias_mtd.axes()[i].upperBound(),
                             m_bias_mtd.axes()[i].width(), m_bias_mtd.axes()[i].bin());
  }
}

BiasWTMeABF::~BiasWTMeABF() = default;

double BiasWTMeABF::sumHistoryHillsAtPosition(const std::vector<double>& pos) const {
  double potential = 0;
  const auto& ax = m_bias_mtd.axes();
  for (size_t i = 0; i < m_history_hills.size(); ++i) {
    potential += m_history_hills[i].hillEnergy(pos, ax);
  }
  return potential;
}

void BiasWTMeABF::updateForce() {
  // update the system force for extended CVs
  BiasExtendedLagrangianBase::updateForce();
  // ABF on extended Lagrangian
  // compute the negative average force and store it
  if (m_count.isInGrid(m_positions)) {
    const size_t N = m_bias_abf.multiplicity();
    const size_t addr = m_count.address(m_positions);
    const size_t current_count = m_count[addr];
    std::vector<double> previous_system_force = m_bias_abf(m_positions);
    for (size_t i = 0; i < N; ++i) {
      m_bias_abf[addr * N + i] = (-m_forces[i] / (current_count + 1.0)) + previous_system_force[i] * (current_count / (current_count + 1.0));
    }
    m_count[addr] += 1;
  }
  // collect info for CZAR
  if (m_zcount.isInGrid(m_real_positions)) {
    const size_t N = m_zgrad.multiplicity();
    const size_t addr = m_zcount.address(m_real_positions);
    const size_t current_zcount = m_zcount[addr];
    std::vector<double> previous_zgrad = m_zgrad(m_real_positions);
    // std::vector<double> tmp_real_f{0, 0};
    for (size_t i = 0; i < N; ++i) {
      m_zgrad[addr * N + i] = (-m_forces[i] / (current_zcount + 1.0)) + previous_zgrad[i] * (current_zcount / (current_zcount + 1.0));
    }
    m_zcount[addr] += 1;
  }
}

std::vector<double> BiasWTMeABF::biasForce(const std::vector<double>& position) {
  // MTD
  if (m_step % m_hill_freq == 0 && m_step > 0) {
    // setup hill
    m_tmp_current_hill.mCenters = position;
    m_tmp_current_hill.mSigmas = m_hill_sigma;
    // well-tempered MTD requires the previous biasing potential
    double previous_bias_V = 0;
    // vector<double>& sum_hills_array = m_mtd_sum_hills.getRawData();
    if (!m_mtd_sum_hills.isInGrid(m_tmp_current_hill.mCenters)) {
      // if the position is outside the boundary, then sum it from the history
      std::cerr << fmt::format("Warning: at step {:d}, a hill is deposited out-of-bound at {::.10f}\n",
                               m_step, m_tmp_current_hill.mCenters);
      previous_bias_V = sumHistoryHillsAtPosition(m_tmp_current_hill.mCenters);
    } else {
      // fast path: if the hill is inside the boundary, then get it from the grid
      const size_t addr = m_mtd_sum_hills.address(m_tmp_current_hill.mCenters);
      previous_bias_V = -m_mtd_sum_hills[addr];
    }
    const double well_tempered_factor = std::exp(-1.0 * previous_bias_V / (boltzmann_constant * m_bias_temperature));
    m_tmp_current_hill.mHeight = well_tempered_factor * m_hill_initial_height;
    // save the current hill
    m_history_hills.push_back(m_tmp_current_hill);
    // project hill
    const auto& point_table = m_bias_mtd.pointTable();
    const auto& point_table_addr = m_bias_mtd.pointTableAddr();
    // TODO: should be parallelized for performance
    for (size_t i = 0; i < m_bias_mtd.histogramSize(); ++i) {
      const size_t& addr = point_table_addr[i];
      double hill_energy = 0;
      // compute hill energy and gradients
      m_tmp_current_hill.hillEnergyGradients(
          point_table[i], m_bias_mtd.axes(),
          m_tmp_hill_gradient, hill_energy);
      for (size_t j = 0; j < m_bias_mtd.dimension(); ++j) {
        m_bias_mtd[addr * m_bias_mtd.dimension() + j] += -m_tmp_hill_gradient[j];
      }
      m_mtd_sum_hills[addr] += -hill_energy;
    }
  }
  std::vector<double> force(m_dof, 0);
  double count = 0;
  if (m_bias_mtd.isInGrid(position)) m_mtd_bias_force = m_bias_mtd(position);
  if (m_count.isInGrid(position)) count = m_count(position);
  // abf_bias_force is actually the sum of instantaneous collective force
//  double abf_force_factor = 0;
  if (count < m_fullsample / 2) {
    m_abf_force_factor = 0;
  } else if (count > m_fullsample) {
    m_abf_force_factor = 1;
  } else {
    m_abf_force_factor = count / m_fullsample;
  }
  if (m_bias_abf.isInGrid(position)) {
    m_abf_bias_force = m_bias_abf(position);
    for (auto& f : m_abf_bias_force) {
      f *= m_abf_force_factor;
    }
  }
  for (size_t i = 0; i < m_dof; ++i) {
    force[i] = m_abf_bias_force[i] + m_mtd_bias_force[i];
  }
  // report the current step
  if (m_step % 1000 == 0) {
    std::cout << fmt::format("step {:10d}: abf force: {:15.10f}, "
                             "abf scaling factor: {:10.5f} ; "
                             "mtd force: {:15.10f}\n",
                             m_step, fmt::join(m_abf_bias_force, " "),
                             m_abf_force_factor, fmt::join(m_mtd_bias_force, " "));
  }
  return force;
}

void BiasWTMeABF::writeHills(std::ostream& os) const {
  if (m_step > 0 && m_step % m_hill_freq == 0) {
    const auto& last_hill = m_tmp_current_hill;
    os << fmt::format(" {:>15d} {:15.10f} {:15.10f} {:15.10f}\n",
                      m_step, fmt::join(last_hill.mCenters, " "),
                      fmt::join(last_hill.mSigmas, " "),
                      last_hill.mHeight);
  }
}

void BiasWTMeABF::writeOutput(string filename, size_t freq) const {
  if (m_step > 0 && m_step % freq == 0) {
    filename += "_" + std::to_string(m_step);
    m_count.writeToFile(filename + ".abf.count");
    m_zcount.writeToFile(filename + ".zcount");
    m_mtd_sum_hills.writeToFile(filename + ".mtd");
    m_bias_mtd.writeToFile(filename + ".mtd.grad");
    m_bias_abf.writeToFile(filename + ".abf.grad");
    m_zgrad.writeToFile(filename + ".zgrad");
    // write CZAR gradients
    const string czar_grad_filename = filename + ".czar.grad";
    std::ofstream ofs_czar_grad(czar_grad_filename.c_str());
    ofs_czar_grad << "# " << m_zgrad.dimension() << '\n';
    for (size_t j = 0; j < m_zgrad.dimension(); ++j) {
      ofs_czar_grad << m_zgrad.axes()[j].infoHeader() << '\n';
    }
    const auto& point_table = m_zgrad.pointTable();
    const auto& point_table_addr = m_zgrad.pointTableAddr();
    for (size_t i = 0; i < m_zcount.histogramSize(); ++i) {
      ofs_czar_grad << fmt::format(" {:15.10f}", fmt::join(point_table[i], " "));
      const size_t& addr = point_table_addr[i];
      const vector<double> log_deriv = m_zcount.getLogDerivative(point_table[i]);
      // merge gradients
      for (size_t j = 0; j < m_zcount.dimension(); ++j) {
        const double grad_value = -1.0 / beta(m_temperature[j]) * log_deriv[j] + m_zgrad[addr * m_zcount.dimension() + j];
        ofs_czar_grad << fmt::format(" {:15.10f}", grad_value);
      }
      ofs_czar_grad << '\n';
    }
  }
}

void BiasWTMeABF::writeTrajectory(std::ostream& os, size_t freq) const {
  if (m_step % freq == 0) {
    os << fmt::format(" {:>15d} {:15.10f} {:15.10f} {:15.10f}\n",
                      m_step, fmt::join(m_real_positions, " "),
                      fmt::join(m_positions, " "),
                      fmt::join(m_bias_force, " "));
  }
}

void BiasWTMeABF::recordStep(const int64_t &step) {
  // then update
  BiasExtendedLagrangianBase::recordStep(step);
}

HarmonicWalls::HarmonicWalls(
    std::vector<double> lower, std::vector<double> upper, double force_constant):
    m_lower(std::move(lower)), m_upper(std::move(upper)), m_constants(m_lower.size(), force_constant) {
  if (m_lower.size() != m_upper.size() || m_lower.size() != m_constants.size()) {
    throw std::runtime_error("Invalid config of harmonic walls.");
  }
  // maybe need more checks
  m_position.assign(m_lower.size(), 0);
  m_force.assign(m_lower.size(), 0);
  m_energy = 0;
}

void HarmonicWalls::update_value(const std::vector<double> &position) {
  m_position = position;
  m_energy = 0;
  for (size_t i = 0; i < m_lower.size(); ++i) {
    if (m_position[i] < m_lower[i]) {
      m_force[i] = -1.0 * m_constants[i] * (m_position[i] - m_lower[i]);
      m_energy += 0.5 * m_constants[i] * (m_position[i] - m_lower[i]) * (m_position[i] - m_lower[i]);
    } else if (m_position[i] > m_upper[i]) {
      m_force[i] = -1.0 * m_constants[i] * (m_position[i] - m_upper[i]);
      m_energy += 0.5 * m_constants[i] * (m_position[i] - m_upper[i]) * (m_position[i] - m_upper[i]);
    } else {
      m_force[i] = 0;
      m_energy += 0;
    }
  }
}

HarmonicWalls::HarmonicWalls(
    std::vector<double> lower, std::vector<double> upper, std::vector<double> force_constants):
    m_lower(std::move(lower)), m_upper(std::move(upper)), m_constants(std::move(force_constants))  {
  if (m_lower.size() != m_upper.size() || m_lower.size() != m_constants.size()) {
    throw std::runtime_error("Invalid config of harmonic walls.");
  }
  // maybe need more checks
  m_position.assign(m_lower.size(), 0);
  m_force.assign(m_lower.size(), 0);
  m_energy = 0;
}
