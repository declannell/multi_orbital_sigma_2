#include "parameters.h"
#include <iostream>
#include "leads_self_energy.h"
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>


void get_self_energies_wba(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &self_energy) {
    for (int r = 0; r < parameters.steps; r++) {
        for (int i = 0; i < parameters.num_orbitals; i++) {
            self_energy.at(r)(i, i) = - parameters.j1 * parameters.gamma;
        }
    }
}

void get_embedding_lesser(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &self_energy,
    std::vector<Eigen::MatrixXcd> &self_energy_lesser, int left_right, int voltage_step) {
    if (left_right == 0) {//corresponds to the left self energy
        for (int r = 0; r < parameters.steps; r++) {
            self_energy_lesser.at(r) = - fermi_function(parameters.energy.at(r) - parameters.voltage_l[voltage_step], parameters) *
                (self_energy.at(r) - self_energy.at(r).adjoint());
        }
    } else {
        for (int r = 0; r < parameters.steps; r++) {
            self_energy_lesser.at(r) = - fermi_function(parameters.energy.at(r) - parameters.voltage_r[voltage_step], parameters) *
                (self_energy.at(r) - self_energy.at(r).adjoint());
        }
    }
}  

void get_gamma(const Parameters &parameters, const  std::vector<Eigen::MatrixXcd> &self_energy, std::vector<Eigen::MatrixXcd> &gamma) {
    for (int r = 0; r < parameters.steps; r++) {
        gamma.at(r) = parameters.j1 * (self_energy.at(r) - self_energy.at(r).adjoint());
    }
}