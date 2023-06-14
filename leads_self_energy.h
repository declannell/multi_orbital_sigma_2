#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>


void get_self_energies_wba(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &self_energy);

void get_embedding_lesser(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &self_energy,
    std::vector<Eigen::MatrixXcd> &self_energy_lesser, int left_right, int voltage_step);

void get_gamma(const Parameters &parameters, const  std::vector<Eigen::MatrixXcd> &self_energy, std::vector<Eigen::MatrixXcd> &gamma);
