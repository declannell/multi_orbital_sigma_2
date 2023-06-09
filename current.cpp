#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
#include "current.h"
using namespace std;

void get_current(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_int_r, std::vector<Eigen::MatrixXcd> &gf_int_l,
    std::vector<Eigen::MatrixXcd> &gamma, double &current, int left_right) {
        dcomp integrand = 0.0;
        if (left_right == 0) {
            for (int r = 0; r < parameters.steps; r++) {
                Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
                Eigen::MatrixXcd spectral_function = parameters.j1 * (gf_int_r.at(r) - gf_int_a);
                Eigen::MatrixXcd trace = fermi_function(parameters.energy.at(r) - parameters.voltage, parameters) *
                        gamma.at(r) * spectral_function + parameters.j1 * gamma.at(r) * gf_int_l.at(r);
                //dcomp test = 0;
  
                //cout << trace.trace() << " " << test << endl;
                //for (int i = 0; i < parameters.num_orbitals; i++) {
                //    integrand += trace(i, i);
                //}
                integrand += trace.trace();
                //cout << integrand << endl;
            }
        } else {
            for (int r = 0; r < parameters.steps; r++) {
                Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
                Eigen::MatrixXcd spectral_function = parameters.j1 * (gf_int_r.at(r) - gf_int_a);
                Eigen::MatrixXcd trace = fermi_function(parameters.energy.at(r) + parameters.voltage, parameters) *
                        gamma.at(r) * spectral_function + parameters.j1 * gamma.at(r) * gf_int_l.at(r);
                integrand += trace.trace();
                //cout << integrand << endl;
            }
        }
        //cout << integrand << endl;
        current = (integrand.real()) * (parameters.energy.at(1) - parameters.energy.at(0)) * 0.000038740155;
    }