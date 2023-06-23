#include "parameters.h"
#include <complex>  //this contains complex numbers and trig functions
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void simple_tokenizer(std::string s, std::string &variable, std::string &value)
{
    std::stringstream ss(s);
    std::string word;
	int count = 0;
    while (ss >> word) {
        //std::cout << word << std::endl;
		if (count == 0) {
			variable = word;
		} else if (count == 1) {
			value = word;
		}
		count++;
    }
}


Parameters Parameters::from_file()
{
	Parameters parameters;
	std::string line, variable, value;
	std::ifstream input_file;

	input_file.open("input_file");

	if(!input_file.is_open())
	{
		std::cout << "The input_file doesn't exist \n" << std::endl;
		std::exit(1);
	} 
	
	while (getline(input_file, line)) {
			//std::cout << line << '\n';
			simple_tokenizer(line, variable, value);
			//std::cout << "The variable name is " << variable << " The value is  " << value << std::endl;
			if (variable == "onsite") {
				parameters.onsite = std::stod(value);
			} else if (variable == "hopping") {
				parameters.hopping = std::stod(value);
			} else if (variable == "num_orbitals") {
				parameters.num_orbitals = std::stoi(value);
			} else if (variable == "chemical_potential") {
				parameters.chemical_potential = std::stod(value);
			} else if (variable == "temperature") {
				parameters.temperature = std::stod(value);
			} else if (variable == "e_upper_bound") {
				parameters.e_upper_bound = std::stod(value);
			} else if (variable == "e_lower_bound") {
				parameters.e_lower_bound = std::stod(value);
			} else if (variable == "hubbard_interaction") {
				parameters.hubbard_interaction = std::stod(value);
			} else if (variable == "voltage_step") {
				parameters.voltage_step = std::stoi(value);
			} else if (variable == "self_consistent_steps") {
				parameters.self_consistent_steps = std::stod(value);
			} else if (variable == "niv_points") {
				parameters.niv_points = std::stoi(value);
			} else if (variable == "niv_start") {
                parameters.niv_start = std::stoi(value);
            } else if (variable == "delta_v") {
				parameters.delta_v = std::stod(value);
			} else if (variable == "delta_gf") {
				parameters.delta_gf = std::stod(value);
			} else if (variable == "spin_up_occup") {
				parameters.convergence = std::stod(value);
			} else if (variable == "gamma") {
				parameters.gamma = std::stod(value);
			} else if (variable == "steps") {
				parameters.steps = std::stoi(value);
			} else if (variable == "print_gf") {
				std::istringstream(value) >> parameters.print_gf;
			} else if (variable == "spin_polarised") {
				std::istringstream(value) >> parameters.spin_polarised;
            } else if (variable == "kk_se") {
				parameters.kk_se = std::stoi(value);
			} else if (variable == "convergence") {
				parameters.convergence = std::stod(value);
			} else if (variable == "e_upper_bound_tilde") {
				parameters.e_upper_bound_tilde = std::stod(value);
			} else if (variable == "e_lower_bound_tilde") {
				parameters.e_lower_bound_tilde = std::stod(value);
			} else if (variable == "multiple") {
				parameters.multiple = std::stoi(value);
			}  else if (variable == "multiple_grids") {
				parameters.multiple_grids = std::stoi(value);
			}
	}
	input_file.close();
	
	parameters.voltage_l.resize(parameters.niv_points);
	parameters.voltage_r.resize(parameters.niv_points);
	for (int i = 0; i < parameters.niv_points; i++) {
		parameters.voltage_l.at(i) = parameters.delta_v * (double)(i);
		parameters.voltage_r.at(i) = -parameters.delta_v * (double)(i);
	}


	parameters.energy.resize(parameters.steps);

	parameters.j1 = -1;
	parameters.j1 = sqrt(parameters.j1);

	parameters.delta_energy =
	    (parameters.e_upper_bound - parameters.e_lower_bound) / (double)parameters.steps;

	for (int i = 0; i < parameters.steps; i++) {
		parameters.energy.at(i) = parameters.e_lower_bound + parameters.delta_energy * (double)i;
	}
	return parameters;
}

double fermi_function(double energy, const Parameters& parameters) {
	if (parameters.temperature == 0) {
		if (energy < parameters.chemical_potential) {
			return 1.0;
		} else {
			return 0.0;
		}
	} else {
		return 1.0 / (1.0 + exp((energy - parameters.chemical_potential) / parameters.temperature));
	}
}
//Parameters params = Parameters::from_file();

void print_parameters(Parameters& parameters)
{
	std::cout << "chemical_potential = " << parameters.chemical_potential << std::endl;
	//std::cout << "parameters.steps = " << parameters.steps << std::endl;
	std::cout << "num_orbitals = " << parameters.num_orbitals << std::endl;
	std::cout << "hopping = " << parameters.hopping << std::endl;
	std::cout << "onsite = " << parameters.onsite << std::endl;
	std::cout << "gamma = " << parameters.gamma << std::endl;
	std::cout << "steps = " << parameters.steps << std::endl;
	std::cout << "e_upper_bound = " << parameters.e_upper_bound << std::endl;
	std::cout << "e_lower_bound = " << parameters.e_lower_bound << std::endl;
	std::cout << "num_orbitals = " << parameters.num_orbitals << std::endl;
	std::cout << "temperature = " << parameters.temperature << std::endl;
	std::cout << "delta_gf = " << parameters.delta_gf << std::endl;
	std::cout << "parameters.convergence is " << parameters.convergence << std::endl;
	std::cout << "e_upper_bound_tilde = " << parameters.e_upper_bound_tilde << std::endl;
	std::cout << "e_lower_bound_tilde = " << parameters.e_lower_bound_tilde << std::endl;
	std::cout << "multiple = " << parameters.multiple << std::endl;
	std::cout << "multiple_grids = " << parameters.multiple_grids << std::endl;
}
