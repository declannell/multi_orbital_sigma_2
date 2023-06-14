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

int main(int argc, char **argv)
{
    //std::cout << argv[1] << std::endl;
	Parameters parameters = Parameters::from_file();
    std::cout << "created parameters \n";
	print_parameters(parameters);
	//std::vector<double> kx(parameters.num_kx_points, 0);
	//std::vector<double> ky(parameters.num_ky_points, 0);
    //get_momentum_vectors(kx, ky, parameters);
	std::vector<double> current_coherent(parameters.niv_points, 0), current_left(parameters.niv_points, 0), current_right(parameters.niv_points, 0);

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
		current_left.at(m) = get_current_mw(parameters, gf_retarded, gf_lesser, self_energy_left, 0, m);
		current_right.at(m) = get_current_mw(parameters, gf_retarded, gf_lesser, self_energy_right, 1, m);

		std::cout << "The voltage is " << parameters.voltage_l[m] - parameters.voltage_r[m] << ". The coherent current is " << current_coherent.at(m) 
			<< ". The left current is " << current_left.at(m) << ". The right current is " << current_right.at(m) << endl;

		std::ostringstream ossgf_lesser;
		ossgf_lesser << m << ".gf_l.dat";
		std::string var_gf_lesser = ossgf_lesser.str();
		std::ofstream gf_l_file;
		gf_l_file.open(var_gf_lesser);
		for (int r = 0; r < parameters.steps; r++) {
			gf_l_file << parameters.energy.at(r) << "  " << gf_lesser.at(r)(1, 0).real() << "   " << gf_lesser.at(r)(1, 0).imag() << "\n";
		}
		gf_l_file.close();	

		std::ostringstream ossgf_retarded;
		ossgf_retarded << m << ".gf_r.dat";
		std::string var_gf_retarded = ossgf_retarded.str();
		std::ofstream gf_retarded_file;
		gf_retarded_file.open(var_gf_retarded);
		for (int r = 0; r < parameters.steps; r++) {
			gf_retarded_file << parameters.energy.at(r) << "  " << gf_retarded.at(r)(1, 0).real() << "   " << gf_retarded.at(r)(1, 0).imag() << "\n";
		}
		gf_retarded_file.close();	

		std::ostringstream ossgf;
		ossgf << m << ".mb_se_r.dat";
		std::string var = ossgf.str();
		std::ofstream se_mb_r_file;
		se_mb_r_file.open(var);
		for (int r = 0; r < parameters.steps; r++) {
			se_mb_r_file << parameters.energy.at(r) << "  " << self_energy_mb_r.at(r)(0, 0).real() << "   " << self_energy_mb_r.at(r)(0, 0).imag() << "\n";
		}
		se_mb_r_file.close();		

		std::ostringstream ossgf_l;
		ossgf_l << m << ".mb_se_l.dat";
		std::string var_l = ossgf_l.str();
		std::ofstream se_mb_l_file;
		se_mb_l_file.open(var_l);
		for (int r = 0; r < parameters.steps; r++) {
			se_mb_l_file << parameters.energy.at(r) << "  " << self_energy_mb_l.at(r)(0, 0).real() << "   " << self_energy_mb_l.at(r)(0, 0).imag() << "\n";
		}
		se_mb_l_file.close();	
	}




    //std::vector<Eigen::Matrix
    ////std::ofstream my_file;
	////my_file.open("gf_int.dat");
	////for(int r = 0; r < 1; r++) {
	////	cout <<  gf_non_int_r_up.at(r) << "\n";
	////}
    ////my_file.close();
	//std::cout << "wrote file \n";

}