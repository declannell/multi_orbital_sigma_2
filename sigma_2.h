
#include </usr/include/eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstdlib>
#include "interacting_gf.h"
#include "leads_self_energy.h"
#include "parameters.h"
#include "mb_self_energy.h"



#include </usr/include/eigen3/Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstdlib>
#include "interacting_gf.h"
#include "leads_self_energy.h"
#include "parameters.h"
#include "mb_self_energy.h"



dcomp integrate(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, const int r);

void self_energy_2nd_order(const Parameters& parameters, std::vector<Eigen::MatrixXcd> &gf_lesser, std::vector<Eigen::MatrixXcd> &gf_retarded, 
	std::vector<Eigen::MatrixXcd> &se_retarded, std::vector<Eigen::MatrixXcd> &se_lesser, int voltage_step, int count);


double get_prefactor(const int i, const int j, const int r, const int voltage_step, const Parameters &parameters,
	const std::vector<dcomp> &fermi_up_1, const std::vector<dcomp> &fermi_up_2, const std::vector<dcomp> &fermi_up_3);

double integrate_equilibrium(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, 
	const int r, const std::vector<dcomp> &fermi_up_1, const std::vector<dcomp> &fermi_up_2, const std::vector<dcomp> &fermi_up_3, int voltage_step);

void self_energy_2nd_order_kramers_kronig(const Parameters &parameters, const std::vector<std::vector<std::vector<dcomp>>>  &gf_lesser, 
	const std::vector<std::vector<std::vector<dcomp>>>  &gf_retarded, std::vector<Eigen::MatrixXcd> &se_retarded,
	 std::vector<Eigen::MatrixXcd> &se_lesser, int voltage_step, int count);

void impurity_solver_sigma_2(const Parameters &parameters, const int voltage_step, 
    std::vector<Eigen::MatrixXcd> &gf_lesser, std::vector<Eigen::MatrixXcd> &gf_retarded, std::vector<Eigen::MatrixXcd> &se_lesser,
	std::vector<Eigen::MatrixXcd> &se_retarded);

void self_energy_2nd_order_kk_diag_approx(const Parameters &parameters, const std::vector<std::vector<std::vector<dcomp>>>  &gf_lesser, 
	const std::vector<std::vector<std::vector<dcomp>>>  &gf_retarded, std::vector<Eigen::MatrixXcd> &se_retarded,
	std::vector<Eigen::MatrixXcd> &se_lesser, int voltage_step, int count);


