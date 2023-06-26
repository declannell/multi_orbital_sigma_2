#pragma once
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


void get_se_mb(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
    std::vector<Eigen::MatrixXcd> &self_energy_mb_r, std::vector<Eigen::MatrixXcd> &self_energy_mb_l);

void get_self_energy(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_retarded, std::vector<Eigen::MatrixXcd> &gf_lesser,
    const std::vector<Eigen::MatrixXcd> &self_energy_left, const std::vector<Eigen::MatrixXcd> &self_energy_right, std::vector<Eigen::MatrixXcd> &self_energy_mb_r,
    std::vector<Eigen::MatrixXcd> &self_energy_mb_l, int voltage_step);

void get_difference_self_energy(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &self_energy_mb_up,
	 std::vector<Eigen::MatrixXcd> &old_self_energy_mb_up, double &difference);

double absolute_value(double num1);

void restructure_gf(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf, std::vector<std::vector<std::vector<dcomp>>> &gf_restruct);

void get_se_mb_diagonal(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
    std::vector<Eigen::MatrixXcd> &self_energy_mb_r, std::vector<Eigen::MatrixXcd> &self_energy_mb_l, int count);

int get_closest_index(const Parameters &parameters, double E_t);

double get_delta(const Parameters &parameters, int i, int j);


//void get_f_function(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
//    std::vector<Eigen::MatrixXd> &f_function);

//void get_f_function(const Parameters &parameters, const std::vector<std::vector<std::vector<dcomp>>>  &gf_retarded,
//    const std::vector<std::vector<std::vector<dcomp>>>  &gf_lesser, 
//    std::vector<std::vector<std::vector<dcomp>>>  &f_function);

//void get_se_mb_kk(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
//    std::vector<Eigen::MatrixXcd> &self_energy_mb_r, std::vector<Eigen::MatrixXcd> &self_energy_mb_l);

//double kramer_kronig_relation(const Parameters& parameters, std::vector<double>& impurity_self_energy_imag, int r);
