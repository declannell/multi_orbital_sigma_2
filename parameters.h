#pragma once
#include <complex>
#include <cmath>
#include <vector>

typedef std::complex<double> dcomp;


struct Parameters
{
    double chemical_potential;
    int steps; // number of energy points we take
    std::vector<double> energy;
    dcomp j1; // this is a complex number class defined within the complex library
    int num_orbitals;
    double voltage;
    static Parameters from_file();
};




double fermi_function(double energy, const Parameters &parameters);
void print_parameters(Parameters& parameters);
