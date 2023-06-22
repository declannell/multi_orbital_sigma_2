#include "parameters.h"
#include "leads_self_energy.h"
#include "mb_self_energy.h"
#include "sigma_2.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include </usr/include/eigen3/Eigen/Dense>
#include "interacting_gf.h"
#include <limits>
#include <iomanip> 
#include <sstream>

double absolute_value(double num1) {
	return std::sqrt((num1 ) * (num1));
}

dcomp integrate(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, const int r)
{
	double delta_energy = (parameters.e_upper_bound - parameters.e_lower_bound) / (double)parameters.steps;
	dcomp result = 0;
	for (int i = 0; i < parameters.steps; i++) {
		for (int j = 0; j < parameters.steps; j++) {
			if (((i + j - r) > 0) && ((i + j - r) < parameters.steps)) {
				//this integrates the equation in PHYSICAL REVIEW B 74, 155125 2006
				//I say the green function is zero outside e_lower_bound and e_upper_bound. This means I need the final green function in the integral to be within an energy of e_lower_bound
				//and e_upper_bound. The index of 0 corresponds to e_lower_bound. Hence we need i+J-r>0 but in order to be less an energy of e_upper_bound we need i+j-r<steps.
				//These conditions ensure the enrgy of the gf3 greens function to be within (e_upper_bound, e_lower_bound)
				result += (delta_energy / (2.0 * M_PI)) * (delta_energy / (2.0 * M_PI)) * gf_1.at(i) * gf_2.at(j) * gf_3.at(i + j - r);
			}
		}
	}
	return result;
}

void get_difference_self_energy(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &self_energy_mb_up,
	 std::vector<Eigen::MatrixXcd> &old_self_energy_mb_up, double &difference){
	difference = - std::numeric_limits<double>::infinity();
	double old_difference = 0;
	double real_difference = 0, imag_difference = 0;
	for (int r = 0; r < parameters.steps; r++) {
		for (int i = 0; i < parameters.num_orbitals; i++) {
			real_difference = absolute_value(self_energy_mb_up.at(r)(i, i).real() - old_self_energy_mb_up.at(r)(i, i).real());
			imag_difference = absolute_value(self_energy_mb_up.at(r)(i, i).imag() - old_self_energy_mb_up.at(r)(i, i).imag());
			//std::cout << gf_local_up.at(r)(i, j).real() << " " << old_green_function.at(r)(i, j).real() << std::endl;
			//std::cout << real_difference << "  " << imag_difference << "  "  << difference << "\n";
			difference = std::max(difference, std::max(real_difference, imag_difference));
			old_self_energy_mb_up.at(r)(i, i) = self_energy_mb_up.at(r)(i, i);
			old_difference = difference;
		}
	}
}

void get_se_mb(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
    std::vector<Eigen::MatrixXcd> &self_energy_mb_r, std::vector<Eigen::MatrixXcd> &self_energy_mb_l, int count) {
    
    std::vector<Eigen::MatrixXcd> gf_advanced(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
        gf_greater(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));
    for (int r = 0; r < parameters.steps; r++) {
        gf_advanced.at(r) = gf_retarded.at(r).adjoint();
        gf_greater.at(r) = gf_retarded.at(r) - gf_advanced.at(r) + gf_lesser.at(r);
    }

   
    for (int r = 0; r < parameters.steps; r++) {
        dcomp result_retarded = 0, result_lesser;
        //std::cout << "The value of r is " << r << std::endl;
        for (int t = 0; t < parameters.num_orbitals; t++) {
            for (int p = 0; p < parameters.num_orbitals; p++) {
                for (int n = 0; n < parameters.num_orbitals; n++) {
                    for (int m = 0; m < parameters.num_orbitals; m++) {
                        for (int q = 0; q < parameters.num_orbitals; q++) {
                            for (int s = 0; s < parameters.num_orbitals; s++) {//spin degenerate so i multiple the first term by 2

	                            for (int i = 0; i < parameters.steps; i++) {
		                            for (int j = 0; j < parameters.steps; j++) {
			                            if (((i + j - r) > 0) && ((i + j - r) < parameters.steps)) {
				                            result_retarded += (2.0 * (gf_lesser.at(i)(n, m) * gf_lesser.at(j)(p, q) * gf_advanced.at(i + j - r)(s, t)
                                                 + gf_lesser.at(i)(n, m) * gf_retarded.at(j)(p, q) * gf_lesser.at(i + j - r)(s, t)
                                                 + gf_retarded.at(i)(n, m) * gf_lesser.at(j)(p, q) * gf_lesser.at(i + j - r)(s, t)
                                                 + gf_retarded.at(i)(n, m) * gf_retarded.at(j)(p, q) * gf_lesser.at(i + j - r)(s, t)) 
                                            - gf_lesser.at(i)(n, q) * gf_lesser.at(j)(p, m) * gf_advanced.at(i + j - r)(s, t)
                                            - gf_lesser.at(i)(n, q) * gf_retarded.at(j)(p, m) * gf_lesser.at(i + j - r)(s, t)
                                            - gf_retarded.at(i)(n, q) * gf_lesser.at(j)(p, m) * gf_lesser.at(i + j - r)(s, t)
                                            - gf_retarded.at(i)(n, q) * gf_retarded.at(j)(p, m) * gf_lesser.at(i + j - r)(s, t));

                                            result_lesser += (2.0 * (gf_lesser.at(i)(n, m) * gf_lesser.at(j)(p, q) * gf_greater.at(i + j - r)(s, t)) 
                                            - gf_lesser.at(i)(n, q) * gf_lesser.at(j)(p, m) * gf_greater.at(i + j - r)(s, t));
			                            }
		                            }
	                            }
                            }
                        }
                    }
                }
            }
        }

        if (count == 0) {
            for (int i = 0; i < parameters.num_orbitals; i++) {
                for (int j = 0; j < parameters.num_orbitals; j++) {
                    self_energy_mb_r.at(r)(i, j) = result_retarded * (parameters.delta_energy / (2.0 * M_PI)) * (parameters.delta_energy / (2.0 * M_PI)) * 
                        parameters.hubbard_interaction * parameters.hubbard_interaction;
                    self_energy_mb_l.at(r)(i, j) = result_lesser * (parameters.delta_energy / (2.0 * M_PI)) * (parameters.delta_energy / (2.0 * M_PI)) * 
                        parameters.hubbard_interaction * parameters.hubbard_interaction;
                }
            }
        } else if (count < 5) {
            //std::cout << "0.5 convergence method \n";
            for (int i = 0; i < parameters.num_orbitals; i++) {
                for (int j = 0; j < parameters.num_orbitals; j++) {
                    self_energy_mb_r.at(r)(i, j) = 0.125 * result_retarded * (parameters.delta_energy / M_PI) * (parameters.delta_energy / (M_PI)) * 
                        parameters.hubbard_interaction * parameters.hubbard_interaction + 0.5 * self_energy_mb_r.at(r)(i, j);
                    self_energy_mb_l.at(r)(i, j) = 0.125 * result_lesser * (parameters.delta_energy / (M_PI)) * (parameters.delta_energy / (M_PI)) * 
                        parameters.hubbard_interaction * parameters.hubbard_interaction + 0.5 * self_energy_mb_l.at(r)(i, j);
                }
            }
        } else {
            //std::cout << "0.9 convergence method \n";
            for (int i = 0; i < parameters.num_orbitals; i++) {
                for (int j = 0; j < parameters.num_orbitals; j++) {
                    self_energy_mb_r.at(r)(i, j) = 0.025 * result_retarded * (parameters.delta_energy / M_PI) * (parameters.delta_energy / (M_PI)) * 
                        parameters.hubbard_interaction * parameters.hubbard_interaction + 0.9 * self_energy_mb_r.at(r)(i, j);
                    self_energy_mb_l.at(r)(i, j) = 0.025 * result_lesser * (parameters.delta_energy / (M_PI)) * (parameters.delta_energy / (M_PI)) * 
                        parameters.hubbard_interaction * parameters.hubbard_interaction + 0.9 * self_energy_mb_l.at(r)(i, j);
                }
            }
        }
    }
}

void get_f_function(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
    std::vector<Eigen::MatrixXd> &f_function){

        for (int r = 0; r < parameters.steps; r++) {
            for (int i = 0; i < parameters.num_orbitals; i++) {
                for (int j = 0; j < parameters.num_orbitals; j++) {
                    f_function.at(r)(i, j) = - gf_lesser.at(r)(i, j).imag() * 0.5 / gf_retarded.at(r)(i, j).imag();
                    if (f_function.at(r)(i, j) != f_function.at(r)(i, j)) {
                        f_function.at(r)(i, j) = 1.0;
                    }
                    std::cout << f_function.at(r)(i, j) << " " << i << " " << j << std::endl;
                }
            }
        }
}

void get_f_function(const Parameters &parameters, const std::vector<std::vector<std::vector<dcomp>>>  &gf_retarded,
    const std::vector<std::vector<std::vector<dcomp>>>  &gf_lesser, 
    std::vector<std::vector<std::vector<dcomp>>>  &f_function){

        for (int r = 0; r < parameters.steps; r++) {
            for (int i = 0; i < parameters.num_orbitals; i++) {
                for (int j = 0; j < parameters.num_orbitals; j++) {
                    f_function.at(i).at(j).at(r) = - gf_lesser.at(i).at(j).at(r).imag() * 0.5 / gf_retarded.at(i).at(j).at(r).imag();
                    if (f_function.at(i).at(j).at(r) != f_function.at(i).at(j).at(r)) {
                        f_function.at(i).at(j).at(r) = 0.5;
                    }
                    //std::cout << f_function.at(i).at(j).at(r) << " " << i << " " << j << std::endl;
                }
            }
        }
}



void get_se_mb_kk(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf_retarded, const std::vector<Eigen::MatrixXcd> &gf_lesser, 
    std::vector<Eigen::MatrixXcd> &self_energy_mb_r, std::vector<Eigen::MatrixXcd> &self_energy_mb_l) {
    
    std::vector<Eigen::MatrixXd> f_function(parameters.steps, Eigen::MatrixXd::Zero(parameters.num_orbitals, parameters.num_orbitals));
    std::vector<Eigen::MatrixXcd> gf_advanced(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
        gf_greater(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

    for (int r = 0; r < parameters.steps; r++) {
        gf_advanced.at(r) = gf_retarded.at(r).adjoint();
        gf_greater.at(r) = gf_retarded.at(r) - gf_advanced.at(r) + gf_lesser.at(r);
    }

    std::vector<double> result_retarded(parameters.steps);
    std::vector<dcomp> result_lesser(parameters.steps);
    get_f_function(parameters, gf_retarded, gf_lesser, f_function);

    if (parameters.num_orbitals == 0) {
        //for (int r = 0; r < parameters.steps; r++) {
        //    for (int i = 0; i < parameters.steps; i++) {
        //        for (int j = 0; j < parameters.steps; j++) {
        //            if (((i + j - r) > 0) && ((i + j  - r) < parameters.steps)) {
        //                result_retarded.at(r) += gf_retarded.at(i)(0, 0).imag() * gf_retarded.at(j)(0, 0).imag() * gf_retarded.at(i + j - r)(0, 0).imag()
        //                    * (f_function.at(i)(0, 0) * f_function.at(j)(0, 0) + f_function.at(i + j - r)(0, 0)) * (1.0 - f_function.at(i)(0, 0) 
        //                    - f_function.at(j)(0, 0));
        //            }
        //        }
        //    }
        //}

        for (int r = 0; r < parameters.steps; r++) {
            for (int t = 0; t < parameters.num_orbitals; t++){
                for (int p = 0; p < parameters.num_orbitals; p++) {
                    for (int n = 0; n < parameters.num_orbitals; n++) {
                        for (int m = 0; m < parameters.num_orbitals; m++) {
                            for (int q = 0; q < parameters.num_orbitals; q++) {
                                for (int s = 0; s < parameters.num_orbitals; s++) {

                                    for (int i = 0; i < parameters.steps; i++) {
                                        for (int j = 0; j < parameters.steps; j++) {
                                            if (((i + j - r) > 0) && ((i + j  - r) < parameters.steps)) {
                                                result_retarded.at(r) += 2.0 * gf_retarded.at(i)(n, m).imag() * gf_retarded.at(j)(p, q).imag() * gf_retarded.at(i + j - r)(s, t).imag()
                                                    * (f_function.at(i)(n, m) * f_function.at(j)(p, q) + f_function.at(i + j - r)(s, t)) * (1.0 - f_function.at(i)(n, m) 
                                                    - f_function.at(j)(p, q)) - gf_retarded.at(i)(n, q).imag() * gf_retarded.at(j)(p, m).imag() * gf_retarded.at(i + j - r)(s, t).imag()
                                                    * (f_function.at(i)(n, q) * f_function.at(j)(p, m) + f_function.at(i + j - r)(s, t)) * (1.0 - f_function.at(i)(n, q) 
                                                    - f_function.at(j)(p, m));
                                            }
                                        }
                                    }

                                }
                            }
                        }
                    }
                }

            }
        }

        for (int r = 0; r < parameters.steps; r++) {
            double real = kramer_kronig_relation(parameters, result_retarded, r);
            self_energy_mb_r.at(r)(0, 0) = (real + parameters.j1 * result_retarded.at(r)) * parameters.delta_energy * parameters.delta_energy;
        }

    } else {
        for (int r = 0; r < parameters.steps; r++) {
            for (int t = 0; t < parameters.num_orbitals; t++){
                for (int p = 0; p < parameters.num_orbitals; p++) {
                    for (int n = 0; n < parameters.num_orbitals; n++) {
                        for (int m = 0; m < parameters.num_orbitals; m++) {
                            for (int q = 0; q < parameters.num_orbitals; q++) {
                                for (int s = 0; s < parameters.num_orbitals; s++) {

                                    for (int i = 0; i < parameters.steps; i++) {
                                        for (int j = 0; j < parameters.steps; j++) {
                                            if (((i + j - r) > 0) && ((i + j  - r) < parameters.steps)) {
                                                result_retarded.at(r) -= 2.0 * gf_retarded.at(i)(n, m).imag() * gf_retarded.at(j)(p, q).imag() * gf_retarded.at(i + j - r)(s, t).imag()
                                                    * (f_function.at(i)(n, m) * f_function.at(j)(p, q) + f_function.at(i + j - r)(s, t)) * (1.0 - f_function.at(i)(n, m) 
                                                    - f_function.at(j)(p, q)) - gf_retarded.at(i)(n, q).imag() * gf_retarded.at(j)(p, m).imag() * gf_retarded.at(i + j - r)(s, t).imag()
                                                    * (f_function.at(i)(n, q) * f_function.at(j)(p, m) + f_function.at(i + j - r)(s, t)) * (1.0 - f_function.at(i)(n, q) 
                                                    - f_function.at(j)(p, m));



                                                //2.0 * (gf_retarded.at(i)(n, m).imag() * gf_retarded.at(j)(p, q).imag() 
                                                //    * gf_retarded.at(i + j - r)(s, t).imag()) * (f_function.at(i + j -r)(s, t) * (f_function.at(i)(n, m) 
                                                //    + f_function.at(j)(p, q) - 1.0) - f_function.at(i)(n, m) * f_function.at(j)(p, q))
                                                //    - gf_retarded.at(i)(n, q).imag() * gf_retarded.at(j)(p, m).imag() 
                                                //    * gf_retarded.at(i + j - r)(s, t).imag() * 
                                                //    (f_function.at(i + j -r)(s, t) * (f_function.at(i)(n, q) + f_function.at(j)(p, m) - 1.0) 
                                                //    - f_function.at(i)(n, q) * f_function.at(j)(p, m));

                                                result_lesser.at(r) += (2.0 * (gf_lesser.at(i)(n, m) * gf_lesser.at(j)(p, q) * gf_greater.at(i + j - r)(s, t)) 
                                                - gf_lesser.at(i)(n, q) * gf_lesser.at(j)(p, m) * gf_greater.at(i + j - r)(s, t));
	    		                            }
	    	                            }
	                                }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int r = 0; r < parameters.steps; r++) {
        double impurity_self_energy_real = kramer_kronig_relation(parameters, result_retarded, r);
        for (int i = 0; i < parameters.num_orbitals; i++) {
            for (int j = 0; j < parameters.num_orbitals; j++) {
                self_energy_mb_r.at(r)(i, j) = (impurity_self_energy_real + parameters.j1 * result_retarded.at(r)) * parameters.delta_energy * 
                    parameters.delta_energy * parameters.hubbard_interaction * parameters.hubbard_interaction / (4.0 * M_PI);
                self_energy_mb_l.at(r)(i, j) = result_lesser.at(r) * (parameters.delta_energy / (2.0 * M_PI)) * (parameters.delta_energy / (2.0 * M_PI)) * 
                    parameters.hubbard_interaction * parameters.hubbard_interaction;
            }
        }
    }

    }


double kramer_kronig_relation(const Parameters& parameters, std::vector<double>& impurity_self_energy_imag, int r)
{
	double real_self_energy = 0;
	double delta_energy = (parameters.e_upper_bound - parameters.e_lower_bound) / (double)parameters.steps;
	for (int q = 0; q < parameters.steps; q++) {
		if (q != r) {
			real_self_energy += delta_energy * impurity_self_energy_imag.at(q) / (parameters.energy.at(q) - parameters.energy.at(r));
        }
    }
	return real_self_energy / M_PI;
}

        
void restructure_gf(const Parameters &parameters, const std::vector<Eigen::MatrixXcd> &gf, std::vector<std::vector<std::vector<dcomp>>> &gf_restruct) {

    for (int r = 0; r < parameters.steps; r++) {
        for(int i = 0; i < parameters.num_orbitals; i++) {
            for (int j = 0; j < parameters.num_orbitals; j++) {
                gf_restruct.at(i).at(j).at(r) = gf.at(r)(i, j);
            }
        }
    }
}

void get_self_energy(const Parameters &parameters, std::vector<Eigen::MatrixXcd> &gf_retarded, std::vector<Eigen::MatrixXcd> &gf_lesser,
    const std::vector<Eigen::MatrixXcd> &self_energy_left, const std::vector<Eigen::MatrixXcd> &self_energy_right, std::vector<Eigen::MatrixXcd> &self_energy_mb_r,
    std::vector<Eigen::MatrixXcd> &self_energy_mb_l, int voltage_step) {
        	
    std::vector<Eigen::MatrixXcd> old_self_energy_mb(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
        self_energy_left_lesser(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)), 
        self_energy_right_lesser(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));

    get_embedding_lesser(parameters, self_energy_left, self_energy_left_lesser, 0, voltage_step);           
    get_embedding_lesser(parameters, self_energy_right, self_energy_right_lesser, 1, voltage_step);               

	Eigen::MatrixXcd hamiltonian(parameters.num_orbitals, parameters.num_orbitals);
    get_hamiltonian(parameters, hamiltonian);

    old_self_energy_mb.at(0)(0, 0) = 100;// this is just so the difference is not zero on the first step
    double difference = std::numeric_limits<double>::infinity();
    int count = 0;

    get_gf_retarded(parameters, gf_retarded, self_energy_left, self_energy_right, self_energy_mb_r, hamiltonian);
    get_gf_lesser(parameters, gf_retarded, gf_lesser, self_energy_left_lesser, self_energy_right_lesser, self_energy_mb_l);
    
	while (difference > parameters.convergence && count < parameters.self_consistent_steps) {        
		
        //if (parameters.kk_se == true) {
        //    get_se_mb_kk(parameters, gf_retarded, gf_lesser, self_energy_mb_r, self_energy_mb_l);
        //} else {
        //    get_se_mb(parameters, gf_retarded, gf_lesser, self_energy_mb_r, self_energy_mb_l);
        //}
        std::vector<std::vector<std::vector<dcomp>>> gf_lesser_restruct(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
            std::vector<dcomp>(parameters.steps))),
         gf_retarded_restruct(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
            std::vector<dcomp>(parameters.steps)));

        restructure_gf(parameters, gf_lesser, gf_lesser_restruct);
        restructure_gf(parameters, gf_retarded, gf_retarded_restruct);


	    if (parameters.kk_se == 1) {// kramer kronig relation.
	    	std::cout << "using the kramer-kronig relation for second order perturbation theory\n";
	    	self_energy_2nd_order_kramers_kronig(parameters, gf_lesser_restruct, gf_retarded_restruct, self_energy_mb_r, self_energy_mb_l, voltage_step, count);
        } else if (parameters.kk_se == 2) {
	    	std::cout << "using the diagonal kramer-kronig approximation for second order perturbation theory\n";
	    	self_energy_2nd_order_kk_diag_approx(parameters, gf_lesser_restruct, gf_retarded_restruct, self_energy_mb_r, self_energy_mb_l, voltage_step, count);            
        } 
        else {
            std::cout << "Using the brute force relations \n";
	    	get_se_mb(parameters, gf_retarded, gf_lesser, self_energy_mb_r, self_energy_mb_l, count);
	    }

        get_gf_retarded(parameters, gf_retarded, self_energy_left, self_energy_right, self_energy_mb_r, hamiltonian);
        get_gf_lesser(parameters, gf_retarded, gf_lesser, self_energy_left_lesser, self_energy_right_lesser, self_energy_mb_l);

        get_difference_self_energy(parameters, self_energy_mb_r, old_self_energy_mb, difference);
        std::cout << std::setprecision(15) << "The difference is " << difference << ". The count is " << count << std::endl;
		if (difference < parameters.convergence) {
			break;
		}
        count++;
    }
}