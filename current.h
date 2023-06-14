#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
 
double get_current_mw(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_int_r, std::vector<Eigen::MatrixXcd> &gf_int_l,
    std::vector<Eigen::MatrixXcd> &embedding_self_energy, int left_right, int voltage_step);

double get_current_transmission(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_retarded, 
    std::vector<Eigen::MatrixXcd> &self_energy_left, std::vector<Eigen::MatrixXcd> &self_energy_right, int voltage_step);