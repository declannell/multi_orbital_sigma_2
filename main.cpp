#include "parameters.h"
#include "interacting_gf.h"
#include "current.h"
#include "mb_self_energy.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip>  
#include </usr/include/eigen3/Eigen/Dense>
using namespace std;

void print_to_file(Parameters &parameters, string filename, std::vector<Eigen::MatrixXcd> &quantity, int orbital_1, int orbital_2, int voltage_step) {
	std::ostringstream ossgf;
	ossgf  << voltage_step << "_" << filename << "_" << orbital_1 << "_" << orbital_2 << ".dat";
	std::string var = ossgf.str();
	std::ofstream file;
	std::cout << "Printing file " << var <<  "\n";
	file.open(var);
	for (int r = 0; r < parameters.steps; r++) {
		file << parameters.energy.at(r) << "  " << quantity.at(r)(orbital_1, orbital_1).real() << "   " << quantity.at(r)(orbital_1, orbital_2).imag() << "\n";
	}
	file.close();	
}


int main(int argc, char **argv)
{
    //std::cout << argv[1] << std::endl;
	Parameters parameters = Parameters::from_file();
    if (parameters.kk_se != 0) {
        parameters.self_consistent_steps = 1;
        std::cout << "The impurity solver can't be done self consistently with the Kramer-kronig yet \n";
    }
    std::cout << "created parameters \n";
	print_parameters(parameters);
	//std::vector<double> kx(parameters.num_kx_points, 0);
	//std::vector<double> ky(parameters.num_ky_points, 0);
    //get_momentum_vectors(kx, ky, parameters);
	std::vector<double> current_coherent(parameters.niv_points, 0), current_left(parameters.niv_points, 0), current_right(parameters.niv_points, 0);

	//for (int r = 0; r < parameters.steps; r++) {
	//	std::cout << parameters.energy.at(r) << "  " << r << std::endl;
	//}

	for (int m = parameters.niv_start; m < parameters.niv_points; m++) {
		std::vector<Eigen::MatrixXcd> self_energy_mb_r(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
			self_energy_mb_l(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
			self_energy_left(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
			self_energy_right(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
			gf_lesser(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
			gf_retarded(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));		

		get_self_energies_wba(parameters, self_energy_left);
		get_self_energies_wba(parameters, self_energy_right);
			
		get_self_energy(parameters, gf_retarded, gf_lesser, self_energy_left, self_energy_right, self_energy_mb_r, self_energy_mb_l, m);

		
		current_coherent.at(m) = get_current_transmission(parameters, gf_retarded, self_energy_left, self_energy_right, m);
		if (parameters.multiple_grids == 1) {
			current_left.at(m) = get_current_mw(parameters, gf_retarded, gf_lesser, self_energy_left, 0, m);
			current_right.at(m) = get_current_mw(parameters, gf_retarded, gf_lesser, self_energy_right, 1, m);
			std::cout << current_left.at(m) << "  " << current_right.at(m) << "\n";
			current_left.at(m) = get_current_mw_multiple_grids(parameters, gf_retarded, gf_lesser, self_energy_left, 0, m);
			current_right.at(m) = get_current_mw_multiple_grids(parameters, gf_retarded, gf_lesser, self_energy_right, 1, m);
		} else {
			current_left.at(m) = get_current_mw(parameters, gf_retarded, gf_lesser, self_energy_left, 0, m);
			current_right.at(m) = get_current_mw(parameters, gf_retarded, gf_lesser, self_energy_right, 1, m);
		}


		std::cout << "The voltage is " << parameters.voltage_l[m] - parameters.voltage_r[m] << ". The coherent current is " << current_coherent.at(m) 
			<< ". The left current is " << current_left.at(m) << ". The right current is " << current_right.at(m) << endl;

		for (int i = 0; i < parameters.num_orbitals; i++) {
			for (int j = 0; j < parameters.num_orbitals; j++) {
				print_to_file(parameters, "gf_retarded", gf_retarded, i, j , m);
				print_to_file(parameters, "gf_lesser", gf_lesser, i, j , m);
				print_to_file(parameters, "se_retarded", self_energy_mb_r, i, j , m);
				print_to_file(parameters, "se_lesser", self_energy_mb_l, i, j , m);
			}
		}


	}
}