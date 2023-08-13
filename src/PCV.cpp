#include "PCV.h"
#include "Common.h"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace ArithmeticPathCV;

CVBasedPath::CVBasedPath(const std::string& path_filename) {
  // read path definition
  std::ifstream ifs_path(path_filename.c_str());
  if (!ifs_path.is_open()) {
    throw std::runtime_error("Failed to read " + path_filename);
  }
  std::string line;
  std::vector<std::string> fields;
  while (std::getline(ifs_path, line)) {
    fields.clear();
    splitString(line, " ", fields);
    if (fields.empty()) continue;
    std::vector<double> points(m_num_cvs, 0);
    for (size_t i = 0; i < m_num_cvs; ++i) {
      points[i] = std::stod(fields[i]);
    }
    m_path_centers.push_back(std::move(points));
  }
}

void CVBasedPath::updateDistanceToReferenceFrames(const std::vector<double>& point, std::vector<std::vector<double>>& distances) {
  size_t num_elements = point.size();
  for (size_t i_frame = 0; i_frame < distances.size(); ++i_frame) {
    for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
      distances[i_frame][j_elem] = point[j_elem] - m_path_centers[i_frame][j_elem];
    }
  }
}

std::vector<std::vector<double>> CVBasedPath::get_centers() const {
  return m_path_centers;
}

std::vector<double> computeDistances(const std::vector<std::vector<double>>& points) {
  if (points.size() < 2) {
    throw std::runtime_error("Too less points for calculating the distances.\n");
  }
  std::vector<double> distances;
  for (size_t i = 1; i < points.size(); ++i) {
    double sum = 0.0;
    for (size_t j = 0; j < points[i].size(); ++j) {
      const double diff = points[i][j] - points[i-1][j];
      sum += diff * diff;
    }
    // fmt::print("points: {} {}\n", points[i][0], points[i][1]);
    distances.push_back(std::sqrt(sum));
  }
  return distances;
}

PathCV::PathCV(const std::string& path_filename): ArithmeticPathBase(), CVBasedPath(path_filename), m_point(m_num_cvs, double()) {
  // call base class for initialization
  ArithmeticPathBase::initialize(m_num_cvs, m_path_centers.size(), -1.0, std::vector<double>(m_num_cvs, 1.0));
  ArithmeticPathBase::reComputeLambda(computeDistances(m_path_centers));
  m_frame_element_distances.resize(m_path_centers.size(), std::vector<double>(m_num_cvs, 0));
  m_s = 0;
  m_z = 0;
  m_grad_s.resize(m_num_cvs, 0);
  m_grad_z.resize(m_num_cvs, 0);
  m_dsdx.resize(m_path_centers.size(), std::vector<double>(m_num_cvs, 0));
  m_dzdx.resize(m_path_centers.size(), std::vector<double>(m_num_cvs, 0));
}

void PathCV::update_value(double x, double y) {
  m_point[0] = x;
  m_point[1] = y;
  updateDistanceToReferenceFrames(m_point, m_frame_element_distances);
  computeValue(m_frame_element_distances, &m_s, &m_z);
  computeDerivatives(m_frame_element_distances, &m_dsdx, &m_dzdx);
  m_grad_s.assign(m_num_cvs, 0);
  m_grad_z.assign(m_num_cvs, 0);
  for (size_t i = 0; i < m_path_centers.size(); ++i) {
    for (size_t j = 0; j < m_num_cvs; ++j) {
      m_grad_s[j] += m_dsdx[i][j];
      m_grad_z[j] += m_dzdx[i][j];
    }
  }
}
