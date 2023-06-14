#pragma once
#include "parameters.h"
#include "leads_self_energy.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>

void get_gf_retarded(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &self_energy_left,
     const std::vector<Eigen::MatrixXcd> &self_energy_right, const std::vector<Eigen::MatrixXcd> &self_energy_mb_r, const Eigen::MatrixXcd hamiltonian);

void get_gf_lesser(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, std::vector<Eigen::MatrixXcd> &gf_lesser, 
    const std::vector<Eigen::MatrixXcd> &self_energy_left_lesser, const std::vector<Eigen::MatrixXcd> &self_energy_right_lesser, 
    const std::vector<Eigen::MatrixXcd> &self_energy_mb_l);

void get_hamiltonian(Parameters const &parameters, Eigen::MatrixXcd &hamiltonian);
