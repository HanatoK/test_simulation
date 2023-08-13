#ifndef PCV_H
#define PCV_H

#include <vector>
#include <string>
#include "colvar_arithmeticpath.h"


std::vector<double> computeDistances(const std::vector<std::vector<double>>& points);


class CVBasedPath {
public:
  CVBasedPath(const std::string& path_filename);
  std::vector<std::vector<double>> get_centers() const;
protected:
  const size_t m_num_cvs = 2;
  std::vector<std::vector<double>> m_path_centers;
  void updateDistanceToReferenceFrames(const std::vector<double>& point, std::vector<std::vector<double>>& distances);
};


class PathCV: public ArithmeticPathCV::ArithmeticPathBase<double>, public CVBasedPath {
public:
  PathCV(const std::string& path_filename);
  void update_value(double x, double y);
  double get_s() const {return m_s;};
  double get_z() const {return m_z;};
  const auto& get_dsdx() const {return m_grad_s;};
  const auto& get_dzdx() const {return m_grad_z;};
  const auto& get_point() const {return m_point;}
private:
  std::vector<double> m_point;
  std::vector<std::vector<double>> m_frame_element_distances;
  std::vector<double> m_grad_s;
  std::vector<double> m_grad_z;
  double m_s;
  double m_z;
  std::vector<std::vector<double>> m_dsdx;
  std::vector<std::vector<double>> m_dzdx;
};

#endif // PCV_H
