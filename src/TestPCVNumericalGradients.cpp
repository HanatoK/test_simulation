#include "PCV.h"
#include <iostream>

void testPCVNumericalGradients() {
  PathCV pcv("../data/path_new3.txt");
  const double x = -5.0;
  const double y = 3.0;
  pcv.update_value(x, y);
  std::vector<double> ds = pcv.get_dsdx();
  std::vector<double> dz = pcv.get_dzdx();
  // check gradients
  std::cout << "Check gradient of PCV s with respect to x\n";
  double epsilon = 0.1;
  for (int i = 0; i < 3; ++i) {
    pcv.update_value(x + epsilon, y);
    const double ds_next = pcv.get_s();
    pcv.update_value(x - epsilon, y);
    const double ds_prev = pcv.get_s();
    const double num_grad = (ds_next - ds_prev) / (2.0 * epsilon);
    const double diff = num_grad - ds[0];
    std::cout << "  Error = " << diff * diff / (epsilon * epsilon) << std::endl;
    epsilon /= 10.0;
  }
  std::cout << "Check gradient of PCV s with respect to y\n";
  epsilon = 0.1;
  for (int i = 0; i < 3; ++i) {
    pcv.update_value(x, y + epsilon);
    const double ds_next = pcv.get_s();
    pcv.update_value(x, y - epsilon);
    const double ds_prev = pcv.get_s();
    const double num_grad = (ds_next - ds_prev) / (2.0 * epsilon);
    const double diff = num_grad - ds[1];
    std::cout << "  Error = " << diff * diff / (epsilon * epsilon) << std::endl;
    epsilon /= 10.0;
  }
  std::cout << "Check gradient of PCV z with respect to x\n";
  epsilon = 0.1;
  for (int i = 0; i < 3; ++i) {
    pcv.update_value(x + epsilon, y);
    const double ds_next = pcv.get_z();
    pcv.update_value(x - epsilon, y);
    const double ds_prev = pcv.get_z();
    const double num_grad = (ds_next - ds_prev) / (2.0 * epsilon);
    const double diff = num_grad - dz[0];
    std::cout << "  Error = " << diff * diff / (epsilon * epsilon) << std::endl;
    epsilon /= 10.0;
  }
  std::cout << "Check gradient of PCV z with respect to y\n";
  epsilon = 0.1;
  for (int i = 0; i < 3; ++i) {
    pcv.update_value(x, y + epsilon);
    const double ds_next = pcv.get_z();
    pcv.update_value(x, y - epsilon);
    const double ds_prev = pcv.get_z();
    const double num_grad = (ds_next - ds_prev) / (2.0 * epsilon);
    const double diff = num_grad - dz[1];
    std::cout << "  Error = " << diff * diff / (epsilon * epsilon) << std::endl;
    epsilon /= 10.0;
  }
}

int main() {
  testPCVNumericalGradients();
}