#include "parameters.h"
#include "gf.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include </usr/include/eigen3/Eigen/Dense>
using namespace std;


void read_non_interacting_gf(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_non_int_r, int spin) {

    for (int i = 0; i < parameters.num_orbitals; i++) {
        for (int j = 0; j < parameters.num_orbitals; j++) {
            fstream my_file;
            std::ostringstream oss;
            oss << "textfiles/Av-k_ReImGF_ReEne_K_" << i + 1 <<  "_" << j + 1 << "_" << spin << ".dat";
            std::string var = oss.str();
            //std::cout << var << std::endl;
	        my_file.open(var, ios::in);
            int count = 0, energy_index = 0;
		    std::string line;
            //cout << gf_non_int_r.at(i).size() << endl;
		    while (!my_file.eof()) {
		    	my_file >> line;
                energy_index = count / 3;
                //std::cout << energy_index << endl;
                //std::cout << line << " " << count << endl;
                if (count % 3 == 1) {
                    gf_non_int_r.at(energy_index)(i, j) = stod(line);
                } else if (count % 3 == 2) {
                    gf_non_int_r.at(energy_index)(i, j) += parameters.j1 * stod(line);
                }
                count++;
		    }
            my_file.close();            
        }
    }
}


void get_interacting_gf(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_non_int_r, std::vector<std::vector<dcomp>> &sigma_mb_r, 
    std::vector<Eigen::MatrixXcd> &gf_int_r) {

    Eigen::MatrixXcd gf_non_int_r_inverse(parameters.num_orbitals, parameters.num_orbitals);
    for(int r = 0; r < parameters.steps; r++) {
        gf_non_int_r_inverse = gf_non_int_r.at(r).inverse();
        for (int i = 0; i < parameters.num_orbitals; i++) {
            gf_non_int_r_inverse(i, i) -= sigma_mb_r.at(i).at(r);
        }
        gf_int_r.at(r) = gf_non_int_r_inverse.inverse();
    }

}


void read_gamma(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gamma_left, std::vector<Eigen::MatrixXcd> &gamma_right, int spin) {
    for (int i = 0; i < parameters.num_orbitals; i++) {
        for (int j = 0; j < parameters.num_orbitals; j++) {
            fstream my_file_left;
            std::ostringstream oss_left;
            oss_left << "textfiles/Av-k_ReImGammaL_ReEne_K_" << i + 1 <<  "_" << j + 1 << "_" << spin << ".dat";
            std::string var_left = oss_left.str();
            //std::cout << var_left << std::endl;
	        my_file_left.open(var_left, ios::in);
            int count = 0, energy_index = 0;
		    std::string line;
            //cout << gf_non_int_r.at(i).size() << endl;
		    while (!my_file_left.eof()) {
		    	my_file_left >> line;
                //cout << se_left.at(0).size() << "  " << energy_index <<  "  " << se_left.at(energy_index)(i, j) << endl;
                //std::cout << energy_index << endl;
                //std::cout << line << " " << count << endl;
                if (count % 3 == 1) {
                    energy_index = count / 3;
                    gamma_left.at(energy_index)(i, j) = stod(line);
                } 
                count++;
		    }
            my_file_left.close();   

            fstream my_file_right;
            std::ostringstream oss_right;
            oss_right << "textfiles/Av-k_ReImGammaR_ReEne_K_" << i + 1 <<  "_" << j + 1 << "_" << spin << ".dat";
            std::string var_right = oss_right.str();
            //std::cout << var_right << std::endl;
	        my_file_right.open(var_right, ios::in);
            count = 0, energy_index = 0;
            //cout << gf_non_int_r.at(i).size() << endl;
		    while (!my_file_right.eof()) {
		    	my_file_right >> line;
                //std::cout << energy_index << endl;
                //std::cout << line << " " << count << endl;
                if (count % 3 == 1) {
                    energy_index = count / 3;
                    gamma_right.at(energy_index)(i, j) = stod(line);
                } 
                count++;
		    }
            my_file_right.close();            
        }
    }   
}


void get_lesser_se(Parameters &parameters, std::vector<Eigen::MatrixXcd> &se_lesser, std::vector<Eigen::MatrixXcd> &gamma, int left_right) {
    if (left_right == 0) {
        cout << "getting a left self energy \n";
        for (int r = 0; r < parameters.steps; r++) {
            se_lesser.at(r) = parameters.j1 * fermi_function(parameters.energy.at(r) - parameters.voltage, parameters) * gamma.at(r);
        }
    } else {
        cout << "getting a right self energy \n";
        for (int r = 0; r < parameters.steps; r++) {
            se_lesser.at(r) = parameters.j1 * fermi_function(parameters.energy.at(r) + parameters.voltage, parameters) * gamma.at(r);
        }
    }
}

void get_lesser_gf(Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_int_r, std::vector<Eigen::MatrixXcd> &se_left,
     std::vector<Eigen::MatrixXcd> &se_right, std::vector<std::vector<dcomp>> &sigma_mb_l, std::vector<Eigen::MatrixXcd> &gf_int_l) {

    for(int r = 0; r < parameters.steps; r++) {
        Eigen::MatrixXcd sigma_lesser_mat(parameters.num_orbitals, parameters.num_orbitals);
        Eigen::MatrixXcd gf_int_a = gf_int_r.at(r).adjoint(); 
        for (int i = 0; i < parameters.num_orbitals; i++) {
            sigma_lesser_mat(i, i) = sigma_mb_l.at(i).at(r);
            //for (int j = 0; j < parameters.num_orbitals; j++) {
            //    gf_int_a(i, j) = std::conj(gf_int_r.at(r)(j, i));
            //}
        }
        gf_int_l.at(r) = gf_int_r.at(r) * (se_left.at(r) + se_right.at(r) + sigma_lesser_mat) * gf_int_a;
    } 
}