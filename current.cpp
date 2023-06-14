#include "parameters.h"
#include "leads_self_energy.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
#include "current.h"
using namespace std;

double get_current_mw(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_int_r, std::vector<Eigen::MatrixXcd> &gf_int_l,
    std::vector<Eigen::MatrixXcd> &embedding_self_energy, int left_right, int voltage_step) {

    std::vector<Eigen::MatrixXcd> gamma(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

    get_gamma(parameters, embedding_self_energy, gamma);
    dcomp integrand = 0.0;
    if (left_right == 0) {
        for (int r = 0; r < parameters.steps; r++) {
            Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
            Eigen::MatrixXcd spectral_function = parameters.j1 * (gf_int_r.at(r) - gf_int_a);
            Eigen::MatrixXcd trace = fermi_function(parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) *
                    gamma.at(r) * spectral_function + parameters.j1 * gamma.at(r) * gf_int_l.at(r);
            integrand += trace.trace();

			//cout << fermi_function(parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) * spectral_function << "\n \n" << 
            //    parameters.j1 * gf_int_l.at(r) <<  "\n  \n \n";
        }
    } else {
        for (int r = 0; r < parameters.steps; r++) {
            Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
            Eigen::MatrixXcd spectral_function = parameters.j1 * (gf_int_r.at(r) - gf_int_a);
            Eigen::MatrixXcd trace = fermi_function(parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters) *
                    gamma.at(r) * spectral_function + parameters.j1 * gamma.at(r) * gf_int_l.at(r);

            integrand += trace.trace();
        }
    }
    cout << integrand << endl;
    return integrand.real() * parameters.delta_energy * 0.5;
}

double get_current_transmission(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_retarded, 
    std::vector<Eigen::MatrixXcd> &self_energy_left, std::vector<Eigen::MatrixXcd> &self_energy_right, int voltage_step) {
    double integrand = 0;
    std::vector<Eigen::MatrixXcd> gamma_l(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
        gamma_r(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

    get_gamma(parameters, self_energy_left, gamma_l);
    get_gamma(parameters, self_energy_right, gamma_r);

    std::ofstream my_file;
	my_file.open("transmission.dat");
    for (int r =0; r < parameters.steps; r++) {
        Eigen::MatrixXcd gf_int_a = gf_retarded.at(r).adjoint();
        dcomp transmission = (gamma_l.at(r) * gf_retarded.at(r) * gamma_r.at(r) * gf_int_a).trace();
        my_file << parameters.energy.at(r) << "  " << transmission.real() << "  " << transmission.imag() << "\n";
        integrand += transmission.real() * (fermi_function(parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) - 
            fermi_function(parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters));
    }
    my_file.close();

    return (integrand) * parameters.delta_energy * 0.5;
    }