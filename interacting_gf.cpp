#include "parameters.h"
#include "leads_self_energy.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>
#include "interacting_gf.h"
#include <limits>
using namespace std;


void get_gf_retarded(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &self_energy_left,
     const std::vector<Eigen::MatrixXcd> &self_energy_right, const std::vector<Eigen::MatrixXcd> &self_energy_mb_r, const Eigen::MatrixXcd hamiltonian) {

	Eigen::MatrixXcd inverse_gf(parameters.num_orbitals, parameters.num_orbitals),
        energy(parameters.num_orbitals, parameters.num_orbitals);

    for (int r = 0; r < parameters.steps; r++) {
        for (int i = 0; i < parameters.num_orbitals; i++) {
            energy(i, i) = parameters.energy.at(r) + parameters.delta_gf * parameters.j1;
        }
        inverse_gf = energy - hamiltonian - self_energy_left.at(r) - self_energy_right.at(r) - self_energy_mb_r.at(r);
        gf_retarded.at(r) = inverse_gf.inverse();
    }
}

void get_gf_lesser(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, std::vector<Eigen::MatrixXcd> &gf_lesser, 
    const std::vector<Eigen::MatrixXcd> &self_energy_left_lesser, const std::vector<Eigen::MatrixXcd> &self_energy_right_lesser, 
    const std::vector<Eigen::MatrixXcd> &self_energy_mb_l) {

    for (int r = 0; r < parameters.steps; r++) {
        Eigen::MatrixXcd gf_advanced = gf_retarded.at(r).adjoint();
        gf_lesser.at(r) = gf_retarded.at(r) * (self_energy_left_lesser.at(r) + self_energy_right_lesser.at(r) + self_energy_mb_l.at(r)) * 
            gf_advanced;
    }
}


void get_hamiltonian(Parameters const &parameters, Eigen::MatrixXcd &hamiltonian){
        
    for (int i = 0; i < parameters.num_orbitals; i++) {
        for (int j = 0; j < parameters.num_orbitals; j++) {
            hamiltonian(i, j) = parameters.hopping;
        }
        hamiltonian(i, i) = parameters.onsite;
    }
    cout << "The Hamiltonain is \n" << hamiltonian << endl;
}

         

