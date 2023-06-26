//
//#include </usr/include/eigen3/Eigen/Dense>
//#include <iomanip>
//#include <iostream>
//#include <vector>
//#include <cstdlib>
//#include "interacting_gf.h"
//#include "leads_self_energy.h"
//#include "parameters.h"
//#include "mb_self_energy.h"
//
//
//
//dcomp integrate_brute_force(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, const int r)
//{
//
//	dcomp result = 0;
//	for (int i = 0; i < parameters.steps; i++) {
//		for (int j = 0; j < parameters.steps; j++) {
//			if (((i + j - r) > 0) && ((i + j - r) < parameters.steps)) {
//				//this integrates the equation in PHYSICAL REVIEW B 74, 155125 2006
//				//I say the green function is zero outside e_lower_bound and e_upper_bound. This means I need the final green function in the integral to be within an energy of e_lower_bound
//				//and e_upper_bound. The index of 0 corresponds to e_lower_bound. Hence we need i+J-r>0 but in order to be less an energy of e_upper_bound we need i+j-r<steps.
//				//These conditions ensure the enrgy of the gf3 greens function to be within (e_upper_bound, e_lower_bound)
//				result += (parameters.delta_energy / (2.0 * M_PI)) * (parameters.delta_energy / (2.0 * M_PI)) * gf_1.at(i) * gf_2.at(j) * gf_3.at(i + j - r);
//			}
//		}
//	}
//	return result;
//}
//
//void self_energy_2nd_order(const Parameters& parameters, std::vector<Eigen::MatrixXcd> &gf_lesser, std::vector<Eigen::MatrixXcd> &gf_retarded, 
//	std::vector<Eigen::MatrixXcd> &se_retarded, std::vector<Eigen::MatrixXcd> &se_lesser, int voltage_step, int count)
//{
//
//	std::vector<dcomp> gf_retarded_up(parameters.steps), advanced_down(parameters.steps), 
//		gf_lesser_up(parameters.steps), gf_greater_down(parameters.steps);
//
//    std::vector<Eigen::MatrixXcd> gf_advanced(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals)),
//        gf_greater(parameters.steps, Eigen::MatrixXcd::Zero(parameters.num_orbitals, parameters.num_orbitals));
//
//    for (int r = 0; r < parameters.steps; r++) {
//        gf_advanced.at(r) = gf_retarded.at(r).adjoint();
//        gf_greater.at(r) = gf_retarded.at(r) - gf_advanced.at(r) + gf_lesser.at(r);
//		gf_retarded_up.at(r) = (gf_retarded.at(r)(0, 0));
//		advanced_down.at(r) = (gf_advanced.at(r)(0, 0));	
//		gf_lesser_up.at(r) = gf_lesser.at(r)(0, 0);
//		gf_greater_down.at(r) = gf_greater.at(r)(0, 0);			
//    }
//
//
//
//	for (int r = 0; r < parameters.steps; r++){
//
//		se_retarded.at(r)(0, 0) = parameters.hubbard_interaction * parameters.hubbard_interaction
//		    * integrate_brute_force(parameters, gf_retarded_up, gf_retarded_up, gf_lesser_up, r); //this resets the self energy	
//		se_retarded.at(r)(0, 0) += parameters.hubbard_interaction * parameters.hubbard_interaction
//		    * integrate_brute_force(parameters, gf_retarded_up, gf_lesser_up, gf_lesser_up, r); 
//		se_retarded.at(r)(0, 0) += parameters.hubbard_interaction * parameters.hubbard_interaction
//		    * integrate_brute_force(parameters, gf_lesser_up, gf_retarded_up, gf_lesser_up, r); 
//		se_retarded.at(r)(0, 0) += parameters.hubbard_interaction * parameters.hubbard_interaction
//		    * integrate_brute_force(parameters, gf_lesser_up, gf_lesser_up, advanced_down, r); 
//
//		se_lesser.at(r)(0, 0) = parameters.j1 * (parameters.hubbard_interaction * parameters.hubbard_interaction
//		    * integrate_brute_force(parameters, gf_lesser_up, gf_lesser_up, gf_greater_down, r)).imag(); 	
//	}
//}
//
//
//
//double get_prefactor(const int i, const int j, const int r, const int voltage_step, const Parameters &parameters,
//	const std::vector<dcomp> &fermi_up_1, const std::vector<dcomp> &fermi_up_2, const std::vector<dcomp> &fermi_up_3)
//{
//	double prefactor;
//	if (voltage_step == 10){
//		prefactor = fermi_function(parameters.energy.at(i), parameters) * fermi_function(parameters.energy.at(j), parameters) 
//			+ (1 - fermi_function(parameters.energy.at(i), parameters) - fermi_function(parameters.energy.at(j), parameters)) *
//			 fermi_function(parameters.energy.at(j) + parameters.energy.at(i) - parameters.energy.at(r), parameters); 
//	} else {
//		int a = i + j - r;
//		dcomp prefactor_complex = fermi_up_1.at(i) * fermi_up_2.at(j) + 
//			fermi_up_3.at(a) * (1.0 - fermi_up_1.at(i) - fermi_up_2.at(j));
//		prefactor = prefactor_complex.real();
//	}		
//	return prefactor;
//}
//
//double integrate_equilibrium(const Parameters& parameters, const std::vector<dcomp>& gf_1, const std::vector<dcomp>& gf_2, const std::vector<dcomp>& gf_3, 
//	const int r, const std::vector<dcomp> &fermi_up_1, const std::vector<dcomp> &fermi_up_2, const std::vector<dcomp> &fermi_up_3, int voltage_step)
//{
//	double result = 0;
//	for (int i = 0; i < parameters.steps; i++) {
//		for (int j = 0; j < parameters.steps; j++) {
//			if (((i + j - r) > 0) && ((i + j - r) < parameters.steps)) {
//				double prefactor = get_prefactor(i, j, r, voltage_step, parameters, fermi_up_1, fermi_up_2, fermi_up_3);
//				//this integrates the equation in PHYSICAL REVIEW B 74, 155125 2006
//				//I say the green function is zero outside e_lower_bound and e_upper_bound. This means I need the final green function in the integral to be within an energy of e_lower_bound
//				//and e_upper_bound. The index of 0 corresponds to e_lower_bound. Hence we need i+J-r>0 but in order to be less an energy of e_upper_bound we need i+j-r<steps.
//				//These conditions ensure the enrgy of the gf3 greens function to be within (e_upper_bound, e_lower_bound)
//				result += prefactor * (parameters.delta_energy / (M_PI)) * (parameters.delta_energy / (M_PI))
//					 * gf_1.at(i).imag() * gf_2.at(j).imag() * gf_3.at(i + j - r).imag();
//			}
//		}
//	}
//	return result;
//}
//
//void self_energy_2nd_order_kramers_kronig(const Parameters &parameters, const std::vector<std::vector<std::vector<dcomp>>>  &gf_lesser, 
//	const std::vector<std::vector<std::vector<dcomp>>>  &gf_retarded, std::vector<Eigen::MatrixXcd> &se_retarded,
//	std::vector<Eigen::MatrixXcd> &se_lesser, int voltage_step, int count)
//{
//	std::vector<double> impurity_self_energy_imag(parameters.steps); //this is for the kramer-kronig relation. 
//	std::vector<double> impurity_self_energy_real(parameters.steps);//, impurity_self_energy_imag_myid(parameters.steps_myid);
//
//    std::vector<std::vector<std::vector<dcomp>>> fermi_eff(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
//            std::vector<dcomp>(parameters.steps))), gf_advanced(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
//            std::vector<dcomp>(parameters.steps))),  gf_greater(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
//            std::vector<dcomp>(parameters.steps)));
//
//	get_f_function(parameters, gf_retarded, gf_lesser, fermi_eff);
//
//    for (int i = 0; i < parameters.num_orbitals; i++) {
//    	for (int j = 0; j < parameters.num_orbitals; j++) {
//	    	for (int r = 0; r < parameters.steps; r++) {
//        		gf_advanced.at(i).at(j).at(r) = std::conj(gf_retarded.at(j).at(i).at(r));
//        		gf_greater.at(i).at(j).at(r) = gf_retarded.at(i).at(j).at(r) - gf_advanced.at(i).at(j).at(r) + gf_lesser.at(i).at(j).at(r);
//    		}
//		}
//	}
//
//
//	for (int r = 0; r < parameters.steps; r++) {
//		std::cout << "The value of r is " << r << std::endl;
//		for (int t = 0; t < parameters.num_orbitals; t++) {
//			for (int p = 0; p < parameters.num_orbitals; p++){
//				for (int n = 0; n < parameters.num_orbitals; n++) {
//					for (int m = 0; m < parameters.num_orbitals; m++) {
//						for (int q = 0; q < parameters.num_orbitals; q++) {
//							for (int s = 0; s < parameters.num_orbitals; s++) {
//								
//								impurity_self_energy_imag.at(r) += 2.0 * parameters.hubbard_interaction * parameters.hubbard_interaction
//								    * integrate_equilibrium(parameters, gf_retarded.at(n).at(m), gf_retarded.at(p).at(q), gf_retarded.at(s).at(t), r,
//									 fermi_eff.at(n).at(m), fermi_eff.at(p).at(q), fermi_eff.at(s).at(t), voltage_step)
//									- parameters.hubbard_interaction * parameters.hubbard_interaction
//								    * integrate_equilibrium(parameters, gf_retarded.at(n).at(q), gf_retarded.at(p).at(m), gf_retarded.at(s).at(t), r,
//									 fermi_eff.at(n).at(q), fermi_eff.at(p).at(m), fermi_eff.at(s).at(t), voltage_step); 	
//
//								se_lesser.at(r)(0, 0) += 2.0 * parameters.j1 * (parameters.hubbard_interaction * parameters.hubbard_interaction * (integrate_brute_force(parameters,
//								gf_lesser.at(n).at(m), gf_lesser.at(p).at(q), gf_greater.at(s).at(t) , r))).imag()
//								 - parameters.j1 * (parameters.hubbard_interaction * parameters.hubbard_interaction * (integrate_brute_force(parameters,
//								gf_lesser.at(n).at(q), gf_lesser.at(p).at(m), gf_greater.at(s).at(t) , r))).imag(); 
//							}
//						}
//					}
//				}
//			}
//		}
//
//	}
//
//	if (count == 0) {
//		for (int r = 0; r < parameters.steps; r++) {
//			int i = 0, j = 0;
//			impurity_self_energy_real.at(r) = kramer_kronig_relation(parameters, impurity_self_energy_imag, r);
//			for (int i = 0; i <parameters.num_orbitals; i++) {
//				for (int j = 0; j < parameters.num_orbitals; j++) {
//					se_retarded.at(r)(i, j) = impurity_self_energy_real.at(r) + parameters.j1 * impurity_self_energy_imag.at(r);
//					se_lesser.at(r)(i, j) = se_lesser.at(r)(0, 0);
//				}
//			}	
//		}
//	} else {
//		for (int r = 0; r < parameters.steps; r++) {
//			int i = 0, j = 0;
//			impurity_self_energy_real.at(r) = kramer_kronig_relation(parameters, impurity_self_energy_imag, r);
//			for (int i = 0; i <parameters.num_orbitals; i++) {
//				for (int j = 0; j < parameters.num_orbitals; j++) {
//					se_retarded.at(r)(i, j) = 0.5 * (impurity_self_energy_real.at(r) + parameters.j1 * impurity_self_energy_imag.at(r)) + 0.5 * se_retarded.at(r)(i, j);
//					se_lesser.at(r)(i, j) = se_lesser.at(r)(0, 0);
//				}
//			}	
//		}		
//	}
//
//
//	std::cout << "calculated using restructed code \n";
//
//
//}
//
//void self_energy_2nd_order_kk_diag_approx(const Parameters &parameters, const std::vector<std::vector<std::vector<dcomp>>>  &gf_lesser, 
//	const std::vector<std::vector<std::vector<dcomp>>>  &gf_retarded, std::vector<Eigen::MatrixXcd> &se_retarded,
//	std::vector<Eigen::MatrixXcd> &se_lesser, int voltage_step, int count)
//{
//	std::vector<double> impurity_self_energy_imag(parameters.steps); //this is for the kramer-kronig relation. 
//	std::vector<double> impurity_self_energy_real(parameters.steps);//, impurity_self_energy_imag_myid(parameters.steps_myid);
//
//    std::vector<std::vector<std::vector<dcomp>>> fermi_eff(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
//            std::vector<dcomp>(parameters.steps))), gf_advanced(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
//            std::vector<dcomp>(parameters.steps))),  gf_greater(parameters.num_orbitals, std::vector<std::vector<dcomp>>(parameters.num_orbitals, 
//            std::vector<dcomp>(parameters.steps)));
//
//	get_f_function(parameters, gf_retarded, gf_lesser, fermi_eff);
//
//    for (int i = 0; i < parameters.num_orbitals; i++) {
//    	for (int j = 0; j < parameters.num_orbitals; j++) {
//	    	for (int r = 0; r < parameters.steps; r++) {
//        		gf_advanced.at(i).at(j).at(r) = std::conj(gf_retarded.at(j).at(i).at(r));
//        		gf_greater.at(i).at(j).at(r) = gf_retarded.at(i).at(j).at(r) - gf_advanced.at(i).at(j).at(r) + gf_lesser.at(i).at(j).at(r);
//    		}
//		}
//	}
//
//
//	for (int r = 0; r < parameters.steps; r++) {
//		std::cout << "The value of r is " << r << std::endl;
//		for (int p = 0; p < parameters.num_orbitals; p++){
//			for (int n = 0; n < parameters.num_orbitals; n++) {
//				for (int s = 0; s < parameters.num_orbitals; s++) {
//					impurity_self_energy_imag.at(r) += 2.0 * parameters.hubbard_interaction * parameters.hubbard_interaction
//					    * integrate_equilibrium(parameters, gf_retarded.at(n).at(n), gf_retarded.at(p).at(p), gf_retarded.at(s).at(s), r,
//						 fermi_eff.at(n).at(n), fermi_eff.at(p).at(p), fermi_eff.at(s).at(s), voltage_step)
//						- parameters.hubbard_interaction * parameters.hubbard_interaction
//					    * integrate_equilibrium(parameters, gf_retarded.at(n).at(n), gf_retarded.at(p).at(p), gf_retarded.at(s).at(s), r,
//						 fermi_eff.at(n).at(n), fermi_eff.at(p).at(p), fermi_eff.at(s).at(s), voltage_step); 
//					se_lesser.at(r)(0, 0) += 2.0 * parameters.j1 * (parameters.hubbard_interaction * parameters.hubbard_interaction * (integrate_brute_force(parameters,
//					gf_lesser.at(n).at(n), gf_lesser.at(p).at(p), gf_greater.at(s).at(s) , r))).imag()
//					 - parameters.j1 * (parameters.hubbard_interaction * parameters.hubbard_interaction * (integrate_brute_force(parameters,
//					gf_lesser.at(n).at(n), gf_lesser.at(p).at(p), gf_greater.at(s).at(s) , r))).imag(); 
//				}
//			}
//		}
//	}
//
//	for (int r = 0; r < parameters.steps; r++) {
//		int i = 0, j = 0;
//		impurity_self_energy_real.at(r) = kramer_kronig_relation(parameters, impurity_self_energy_imag, r);
//		for (int i = 0; i <parameters.num_orbitals; i++) {
//			se_retarded.at(r)(i, i) = impurity_self_energy_real.at(r) + parameters.j1 * impurity_self_energy_imag.at(r);
//			se_lesser.at(r)(i, i) = se_lesser.at(r)(0, 0);
//		}
//		
//	}
//
//	std::cout << "calculated using approximate kramer-kronig code \n";
//}



//void impurity_solver_sigma_2(const Parameters &parameters, const int voltage_step, 
//    std::vector<Eigen::MatrixXcd> &gf_lesser, std::vector<Eigen::MatrixXcd> &gf_retarded, std::vector<Eigen::MatrixXcd> &se_lesser,
//	std::vector<Eigen::MatrixXcd> &se_retarded)
//{
//	if (parameters.kk_se == true) {// kramer kronig relation.
//		std::cout << "using the kramer-kronig relation for second order perturbation theory\n";
//		self_energy_2nd_order_kramers_kronig(parameters, gf_lesser, gf_retarded, se_retarded, se_lesser, voltage_step);
//	} else {
//		self_energy_2nd_order(parameters, gf_lesser, gf_retarded, se_retarded, se_lesser, voltage_step);
//	}
//}
//