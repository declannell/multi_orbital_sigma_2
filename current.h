#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
 

void get_current(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_int_r, std::vector<Eigen::MatrixXcd> &gf_int_l,
    std::vector<Eigen::MatrixXcd> &gamma, double &current, int left_right);