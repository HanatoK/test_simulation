#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "Axis.h"
#include "Grid.h"
#include "MetaDynamics.h"

#include <vector>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <cstdio>

using namespace std;

class potential {
private:
    // GaMD parameters
    double Vmin;
    double Vmax;
    double Vavg;
    double sigmaV;
    double sigma0;
    double k0p;
    double k0;
    double k;
    double E;
    ofstream ofs_traj;
public:
    MetaDynamics meta;
public:
    potential() {
        Axis phi_ax(-10, 10, 200);
        Axis psi_ax(-10, 10, 200);
        vector<Axis> ax{phi_ax, psi_ax};
        meta = MetaDynamics(ax);
        vector<Hill> hills;
        hills.push_back(Hill({6, 6}, {2, 2}, -3.5));
        hills.push_back(Hill({-4, -6}, {2, 2}, -5));
        hills.push_back(Hill({-9, 9}, {2.5, 2.5}, 12));
        hills.push_back(Hill({-5, 5}, {4, 4}, 4));
        hills.push_back(Hill({-5, 3}, {4, 6}, 4));
        hills.push_back(Hill({9.5, -9.5}, {3, 3}, 10));
        hills.push_back(Hill({4, -4}, {5, 5}, 7));
        hills.push_back(Hill({2, -2}, {4, 4}, 4));
        hills.push_back(Hill({0, 0}, {3, 3}, 3));
        for (const auto& h : hills) {
            meta.addHill(h);
        }
        meta.compute();
        sigma0 = 2.0;
        generateGaMDParameters();
        ofs_traj.open("GaMD.traj");
        ofs_traj << "# x y grad_x grad_y V biasV" << endl;
    }
    potential(const potential& p2) {
        Vmin = p2.Vmin;
        Vmax = p2.Vmax;
        Vavg = p2.Vavg;
        sigmaV = p2.sigmaV;
        sigma0 = p2.sigma0;
        k0p = p2.k0p;
        k0 = p2.k0;
        k = p2.k;
        E = p2.E;
        meta = p2.meta;
        ofs_traj.open("GaMD.traj");
    }
    potential& operator=(const potential& p2) {
        Vmin = p2.Vmin;
        Vmax = p2.Vmax;
        Vavg = p2.Vavg;
        sigmaV = p2.sigmaV;
        sigma0 = p2.sigma0;
        k0p = p2.k0p;
        k0 = p2.k0;
        k = p2.k;
        E = p2.E;
        meta = p2.meta;
        ofs_traj.open("GaMD.traj");
        ofs_traj << "# x y grad_x grad_y V biasV" << endl;
        return *this;
    }
    ~potential() {
        ofs_traj.close();
    }
    void generateGaMDParameters() {
        const vector<double>& raw_data = meta.getRawData();
        Vmin = *std::min_element(raw_data.begin(), raw_data.end());
        Vmax = *std::max_element(raw_data.begin(), raw_data.end());
        Vavg = 0;
        double V2avg = 0;
        for (const auto& i : raw_data) {
            Vavg += i;
            V2avg += i * i;
        }
        Vavg /= raw_data.size();
        V2avg /= raw_data.size();
        sigmaV = std::sqrt(V2avg - Vavg * Vavg);
        // k0'
        k0p = sigma0 / sigmaV * (Vmax - Vmin) / (Vmax - Vavg);
        // k0 = min(1.0, k0')
        k0 = k0p > 1.0 ? 1.0 : k0p;
        // k
        k = k0 * 1.0 / (Vmax - Vmin);
        // E
        E = Vmax;
        printf("GaMD: E = %lf ; sigma0 = %lf ; sigmaV = %lf ; k0' = %lf ; k = %lf\n", E, sigma0, sigmaV, k0p, k);
    }
    tuple<double, double> getGaMDBiasPotentialAndConstant(double potential) {
        // bias potential
        const double bias_potential = 0.5 * k * (E - potential) * (E - potential);
        return make_tuple(bias_potential, k);
    }
    vector<double> getGaMDForces(const double x, const double y, const double z) {
        double potential, bias_potential, k;
        vector<double> grad(2);
        meta.computeEnergyGradients(vector<double>{x, y}, potential, grad);
        tie(bias_potential, k) = getGaMDBiasPotentialAndConstant(potential);
//         const double E = Vmax;
        if (potential < E) {
            for (auto& i : grad) {
                // multiply -1.0 to force
                i *= -1.0 * (1.0 - k * (E - potential));
            }
        } else {
            bias_potential = 0;
        }
        grad.push_back(0.0);
        ofs_traj.setf(ios::fixed);
        const int precision = 12;
        const int width = 18;
        ofs_traj << setprecision(precision);
        ofs_traj << setw(width) << x << setw(width) << y;
        ofs_traj << setw(width) << grad[0] << setw(width) << grad[1];
        ofs_traj << setw(width) << potential << setw(width) << bias_potential << '\n';
        return grad;
    }
};

extern potential p;

double getPotential(const double x, const double y, const double z = 0);
vector<double> getGradients(const double x, const double y, const double z = 0);
vector<double> getGaMDForces(const double x, const double y, const double z = 0);

#endif // POTENTIAL_H
