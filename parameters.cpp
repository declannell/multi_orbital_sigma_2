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
			if (variable == "chemical_potential") {
				parameters.chemical_potential = std::stod(value);
			} else if (variable == "num_orbitals") {
				parameters.num_orbitals = std::stoi(value);
			} else if (variable == "voltage") {
				parameters.voltage = std::stod(value);
			} 
	}
	input_file.close();
	
	//parameters.energy.resize(parameters.steps);

	parameters.j1 = -1;
	parameters.j1 = sqrt(parameters.j1);

	//double delta_energy =
	//    (parameters.e_upper_bound - parameters.e_lower_bound) / (double)parameters.steps;
//
	//for (int i = 0; i < parameters.steps; i++) {
	//	parameters.energy.at(i) = parameters.e_lower_bound + delta_energy * (double)i;
	//}
	return parameters;
}

double fermi_function(double energy, const Parameters& parameters)
{
	return 1.0 / (1.0 + exp((energy - parameters.chemical_potential) / 0.00190006954));
}
//Parameters params = Parameters::from_file();

void print_parameters(Parameters& parameters)
{
	std::cout << "chemical_potential = " << parameters.chemical_potential << std::endl;
	//std::cout << "parameters.steps = " << parameters.steps << std::endl;
	std::cout << "num_orbitals = " << parameters.num_orbitals << std::endl;
}
