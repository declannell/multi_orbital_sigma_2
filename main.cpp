#include "parameters.h"
#include "gf.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip>  
#include "read_sigma.h"
#include </usr/include/eigen3/Eigen/Dense>
#include "current.h"
using namespace std;

int main(int argc, char **argv)
{
    //std::cout << argv[1] << std::endl;
	Parameters parameters = Parameters::from_file();
    std::cout << "created parameters \n";
	print_parameters(parameters);
	//std::vector<double> kx(parameters.num_kx_points, 0);
	//std::vector<double> ky(parameters.num_ky_points, 0);
    //get_momentum_vectors(kx, ky, parameters);
	parameters.voltage = 0.5 * 1.375 / 13.6057039763;
	parameters.num_orbitals = 30;
    std::cout << parameters.num_orbitals  << "  " << parameters.voltage << " " << parameters.chemical_potential << std::endl;

	std::vector<std::vector<dcomp>> sigma_mb_r_up, sigma_mb_r_down, sigma_mb_l_up, sigma_mb_l_down;

	get_sigma(parameters, sigma_mb_r_up, sigma_mb_r_down, sigma_mb_l_up, sigma_mb_l_down);	


    std::vector<Eigen::MatrixXcd> gf_non_int_r_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		gf_non_int_r_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

	read_non_interacting_gf(parameters, gf_non_int_r_up, 1);
	read_non_interacting_gf(parameters, gf_non_int_r_down, 2);

    std::vector<Eigen::MatrixXcd> gf_int_r_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		gf_int_r_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

	get_interacting_gf(parameters, gf_non_int_r_up, sigma_mb_r_up, gf_int_r_up);
	get_interacting_gf(parameters, gf_non_int_r_down, sigma_mb_r_down, gf_int_r_down);	

	cout << "Got the interacting GF \n";
	
    std::vector<Eigen::MatrixXcd> gamma_left_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		gamma_right_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));
    std::vector<Eigen::MatrixXcd> gamma_left_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		gamma_right_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

	read_gamma(parameters, gamma_left_up, gamma_right_up, 1);
	read_gamma(parameters, gamma_left_down, gamma_right_down, 2);

	cout << "Got the gamma matrices \n";

    std::vector<Eigen::MatrixXcd> se_left_lesser_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		se_right_lesser_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));
    std::vector<Eigen::MatrixXcd> se_left_lesser_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		se_right_lesser_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

	get_lesser_se(parameters, se_left_lesser_up, gamma_left_up, 0);
	get_lesser_se(parameters, se_left_lesser_down, gamma_left_down, 0);
	get_lesser_se(parameters, se_right_lesser_up, gamma_right_up, 1);
	get_lesser_se(parameters, se_right_lesser_down, gamma_right_down, 1);

	cout << "Got the lesser self energy \n";

    std::vector<Eigen::MatrixXcd> gf_int_l_up(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
		gf_int_l_down(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

	get_lesser_gf(parameters, gf_int_r_up, se_left_lesser_up, se_right_lesser_up, sigma_mb_l_up, gf_int_l_up);
	get_lesser_gf(parameters, gf_int_r_down, se_left_lesser_down, se_right_lesser_down, sigma_mb_l_down, gf_int_l_down);


	double current_up_l, current_up_r, current_down_l, current_down_r;

	get_current(parameters, gf_int_r_up, gf_int_l_up, gamma_left_up, current_up_l, 0);
	get_current(parameters, gf_int_r_up, gf_int_l_up, gamma_right_up, current_up_r, 1);
	get_current(parameters, gf_int_r_down, gf_int_l_down, gamma_left_down, current_down_l, 0);
	get_current(parameters, gf_int_r_down, gf_int_l_down, gamma_right_down, current_down_r, 1);

	cout << "The spin up left current is " << current_up_l << endl;
	cout << "The spin up right current is " << current_up_r << endl;
	cout << "The spin down left current is " << current_down_l << endl;
	cout << "The spin down right current is " << current_down_r << endl;

	//for (int i = 0; i < parameters.num_orbitals; i++) {
	//	std::ofstream my_file;
	//	std::ostringstream oss;
    //    oss << "se_left_up" << i + 1 <<  "_" << i + 1 << ".dat";
    //    std::string var = oss.str();
	//	my_file.open(var);
	//	for(int r = 0; r < parameters.steps; r++) {
	//		my_file << parameters.energy.at(r) << "  " << se_left_up.at(r)(i, i).real() << " " << se_left_up.at(r)(i, i).imag() << "\n";
	//	}
    //	my_file.close();
	//}
//
	//for (int i = 0; i < parameters.num_orbitals; i++) {
	//	std::ofstream my_file;
	//	std::ostringstream oss;
    //    oss << "se_left_lesser_up" << i + 1 <<  "_" << i + 1 << ".dat";
    //    std::string var = oss.str();
	//	my_file.open(var);
	//	for(int r = 0; r < parameters.steps; r++) {
	//		my_file << parameters.energy.at(r) << "  " << se_left_lesser_up.at(r)(i, i).real() << " " << se_left_lesser_up.at(r)(i, i).imag() << "\n";
	//	}
    //	my_file.close();
	//}
//
	for (int i = 0; i < parameters.num_orbitals; i++) {
		std::ofstream my_file;
		std::ostringstream oss;
        oss << "gf_int_l_up" << i + 1 <<  "_" << i + 1 << ".dat";
        std::string var = oss.str();
		my_file.open(var);
		for(int r = 0; r < parameters.steps; r++) {
			my_file << parameters.energy.at(r) << "  " << gf_int_l_up.at(r)(i, i).real() << " " << gf_int_l_up.at(r)(i, i).imag() << "\n";
		}
    	my_file.close();
	}
	
	for (int i = 0; i < parameters.num_orbitals; i++) {
		std::ofstream my_file;
		std::ostringstream oss;
        oss << "gf_int_l_down" << i + 1 <<  "_" << i + 1 << ".dat";
        std::string var = oss.str();
		my_file.open(var);
		for(int r = 0; r < parameters.steps; r++) {
			my_file << parameters.energy.at(r) << "  " << gf_int_l_down.at(r)(i, i).real() << " " << gf_int_l_down.at(r)(i, i).imag() << "\n";
		}
    	my_file.close();
	}
    ////std::ofstream my_file;
	////my_file.open("gf_int.dat");
	////for(int r = 0; r < 1; r++) {
	////	cout <<  gf_non_int_r_up.at(r) << "\n";
	////}
    ////my_file.close();
	//std::cout << "wrote file \n";

}