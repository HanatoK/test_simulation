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
  std::string line;
  std::vector<std::string> fields;
  while (std::getline(ifs_path, line)) {
    fields.clear();
    splitString(line, " ", fields);
    if (fields.size() == 0) continue;
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
  ArithmeticPathBase::initialize(m_num_cvs, m_path_centers.size(), -1.0, std::vector<double>(m_num_cvs, 0), std::vector<double>(m_num_cvs, 1.0));
  ArithmeticPathBase::reComputeLambda(computeDistances(m_path_centers));
}

void PathCV::updateDistanceToReferenceFrames() {
  CVBasedPath::updateDistanceToReferenceFrames(m_point, frame_element_distances);
}

void PathCV::update_value(double x, double y) {
  m_point[0] = x;
  m_point[1] = y;
  compute();
}
