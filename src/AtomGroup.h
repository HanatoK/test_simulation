#ifndef ATOMGROUP_H
#define ATOMGROUP_H

#include <stdexcept>
#include <vector>

struct AtomGroup {
  size_t m_num_atoms;
  std::vector<double> m_pos_x;
  std::vector<double> m_pos_y;
  std::vector<double> m_pos_z;
  std::vector<double> m_vel_x;
  std::vector<double> m_vel_y;
  std::vector<double> m_vel_z;
  std::vector<double> m_mass;
  AtomGroup(size_t num_atoms, std::vector<double> mass):
    m_num_atoms(num_atoms),
    m_pos_x(num_atoms), m_pos_y(num_atoms), m_pos_z(num_atoms),
    m_vel_x(num_atoms), m_vel_y(num_atoms), m_vel_z(num_atoms),
    m_mass(std::move(mass)) {
      if (m_mass.size() != num_atoms) {
        throw std::runtime_error("The mass vector has a size of " +
        std::to_string(m_mass.size()) + ", which is different from the "
        "number of atoms (" + std::to_string(m_num_atoms) + ").\n");
      }
    }
};

#endif // ATOMGROUP_H
