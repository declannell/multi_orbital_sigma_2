#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>



void read_non_interacting_gf(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_non_int_r, int spin);

void get_interacting_gf(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_non_int_r, std::vector<std::vector<dcomp>> &sigma_mb_r, 
    std::vector<Eigen::MatrixXcd> &gf_int_r);

void read_gamma(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gamma_left, std::vector<Eigen::MatrixXcd> &gamma_right, int spin);

void get_lesser_se(Parameters &parameters, std::vector<Eigen::MatrixXcd> &se_lesser, std::vector<Eigen::MatrixXcd> &sgamma, int left_right);

void get_lesser_gf(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_int_r, std::vector<Eigen::MatrixXcd> &se_left,
     std::vector<Eigen::MatrixXcd> &se_right, std::vector<std::vector<dcomp>> &sigma_mb_l, std::vector<Eigen::MatrixXcd> &gf_int_l);
