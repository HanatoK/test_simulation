#ifndef BOX_H
#define BOX_H

#include <array>
#include <cmath>

class Simulation;

class Box {
public:
  using BoxVectorT = std::array<double, 3>;
  friend class Simulation;
  Box(BoxVectorT vec_x, BoxVectorT vec_y, BoxVectorT vec_z,
      BoxVectorT origin): m_vec_x(vec_x), m_vec_y(vec_y),
      m_vec_z(vec_z), m_origin(origin) {}
  double volume() const {
    const BoxVectorT x_cross_y{
      m_vec_x[1] * m_vec_y[2] - m_vec_x[2] * m_vec_y[1],
      m_vec_x[2] * m_vec_y[0] - m_vec_x[0] * m_vec_y[2],
      m_vec_x[0] * m_vec_y[1] - m_vec_x[1] * m_vec_y[0]
    };
    const double vol = x_cross_y[0] * m_vec_z[0] +
                       x_cross_y[1] * m_vec_z[1] +
                       x_cross_y[2] * m_vec_z[2];
    return std::abs(vol);
  }
  void wrap(double x, double y, double z,
            double& w_x, double& w_y, double& w_z) {
    // TODO
  }
private:
  BoxVectorT m_vec_x;
  BoxVectorT m_vec_y;
  BoxVectorT m_vec_z;
  BoxVectorT m_origin;
};

#endif // BOX_H
