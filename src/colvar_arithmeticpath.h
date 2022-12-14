#ifndef ARITHMETICPATHCV_H
#define ARITHMETICPATHCV_H

#include <type_traits>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <fmt/printf.h>

namespace ArithmeticPathCV {

using std::vector;

template <typename element_type, typename scalar_type>
class ArithmeticPathBase {
public:
    ArithmeticPathBase() {}
    virtual ~ArithmeticPathBase() {}
    virtual void initialize(size_t p_num_elements, size_t p_total_frames, double p_lambda, const vector<element_type>& p_element, const vector<double>& p_weights);
    virtual void updateDistanceToReferenceFrames() = 0;
    virtual void computeValue();
    virtual void computeDerivatives();
    virtual void compute();
    virtual void reComputeLambda(const vector<scalar_type>& rmsd_between_refs);
protected:
    scalar_type lambda;
    vector<scalar_type> weights;
    size_t num_elements;
    size_t total_frames;
    vector< vector<element_type> > frame_element_distances;
    scalar_type s;
    scalar_type z;
    vector<element_type> dsdx;
    vector<element_type> dzdx;
private:
    // intermediate variables
    vector<scalar_type> s_numerator_frame;
    vector<scalar_type> s_denominator_frame;
    scalar_type numerator_s;
    scalar_type denominator_s;
    scalar_type normalization_factor;
};

template <typename element_type, typename scalar_type>
void ArithmeticPathBase<element_type, scalar_type>::initialize(size_t p_num_elements, size_t p_total_frames, double p_lambda, const vector<element_type>& p_element, const vector<double>& p_weights) {
    lambda = p_lambda;
    weights = p_weights;
    num_elements = p_num_elements;
    total_frames = p_total_frames;
    frame_element_distances.resize(total_frames, p_element);
    for (size_t i_frame = 0; i_frame < frame_element_distances.size(); ++i_frame) {
        for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
            frame_element_distances[i_frame][j_elem] = element_type();
        }
    }
    s = scalar_type(0);
    z = scalar_type(0);
    dsdx = p_element;
    dzdx = p_element;
    s_numerator_frame.resize(total_frames, scalar_type(0));
    s_denominator_frame.resize(total_frames, scalar_type(0));
    numerator_s = scalar_type(0);
    denominator_s = scalar_type(0);
    normalization_factor = 1.0 / static_cast<scalar_type>(total_frames - 1);
}

template <typename element_type, typename scalar_type>
void ArithmeticPathBase<element_type, scalar_type>::computeValue() {
    updateDistanceToReferenceFrames();
    numerator_s = scalar_type(0);
    denominator_s = scalar_type(0);
    for (size_t i_frame = 0; i_frame < frame_element_distances.size(); ++i_frame) {
        scalar_type exponent_tmp = scalar_type(0);
        for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
            exponent_tmp += weights[j_elem] * frame_element_distances[i_frame][j_elem] * weights[j_elem] * frame_element_distances[i_frame][j_elem];
        }
        exponent_tmp = std::exp(-1.0 * lambda * exponent_tmp);
        numerator_s += static_cast<scalar_type>(i_frame) * exponent_tmp;
        denominator_s += exponent_tmp;
        s_numerator_frame[i_frame] = static_cast<scalar_type>(i_frame) * exponent_tmp;
        s_denominator_frame[i_frame] = exponent_tmp;
    }
    s = numerator_s / denominator_s * normalization_factor;
    z = -1.0 / lambda * std::log(denominator_s);
}

template <typename element_type, typename scalar_type>
void ArithmeticPathBase<element_type, scalar_type>::compute() {
    computeValue();
    computeDerivatives();
}

template <typename element_type, typename scalar_type>
void ArithmeticPathBase<element_type, scalar_type>::computeDerivatives() {
    for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
        element_type dsdxj_numerator_part1 = std::decay_t<decltype(dsdx[j_elem])>(0.0);
        element_type dsdxj_numerator_part2 = std::decay_t<decltype(dsdx[j_elem])>(0.0);
        element_type dzdxj_numerator = std::decay_t<decltype(dzdx[j_elem])>(0.0);
        for (size_t i_frame = 0; i_frame < frame_element_distances.size(); ++i_frame) {
            const element_type derivative_tmp = -2.0 * lambda * weights[j_elem] * weights[j_elem] * frame_element_distances[i_frame][j_elem];
            dsdxj_numerator_part1 += s_numerator_frame[i_frame] * derivative_tmp;
            dsdxj_numerator_part2 += s_denominator_frame[i_frame] * derivative_tmp;
            dzdxj_numerator += s_denominator_frame[i_frame] * derivative_tmp;
        }
        dsdxj_numerator_part1 *= denominator_s;
        dsdxj_numerator_part2 *= numerator_s;
        dsdx[j_elem] = (dsdxj_numerator_part1 - dsdxj_numerator_part2) / (denominator_s * denominator_s) * normalization_factor;
        dzdx[j_elem] = -1.0 / lambda * dzdxj_numerator / denominator_s;
    }
}

template <typename element_type, typename scalar_type>
void ArithmeticPathBase<element_type, scalar_type>::reComputeLambda(const vector<scalar_type>& rmsd_between_refs) {
    scalar_type mean_square_displacements = 0.0;
    for (size_t i_frame = 1; i_frame < total_frames; ++i_frame) {
        fmt::print("Distance between frame {:3d} and {:3d} is {:12.7f}\n", i_frame, i_frame + 1, rmsd_between_refs[i_frame - 1]);
        mean_square_displacements += rmsd_between_refs[i_frame - 1] * rmsd_between_refs[i_frame - 1];
    }
    mean_square_displacements /= scalar_type(total_frames - 1);
    lambda = 1.0 / mean_square_displacements;
    fmt::print("lambda is updated to {}\n", lambda);
}
}

#endif // ARITHMETICPATHCV_H
