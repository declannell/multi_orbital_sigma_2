#pragma once
#include <complex>
#include <cmath>
#include <vector>
//#include <mpi.h>

typedef std::complex<double> dcomp;


struct Parameters
{
    double onsite;         // onsite energy in the correlated region
    double hopping;        // the hopping the z direction of the scattering region
    int num_orbitals;
    double chemical_potential;
    double temperature;
    double e_upper_bound;       // this is the max energy value
    double e_lower_bound;       // this is the min energy value
    double hubbard_interaction; // this is the hubbard interaction
    int voltage_step;        // voltage step of zero is equilibrium. This is an integer and higher values correspond to a higher potential difference between the two leads.
    double self_consistent_steps; // this is the number of self consistent steps my code needs
    int niv_points;//number of IV points
	int niv_start; //starting bias for the calculation. 0 for equilibrium
    double delta_v; //the voltage step between IV points
    double delta_gf; //delta in the scattering region
    bool leads_3d;//if true, this will attach 3d leads to a 1d scattering region
    double convergence;
    double gamma; //this is the value of the imag part of the self energy in the WBL.
    std::vector<double> voltage_r;
    std::vector<double> voltage_l;
    int steps; // number of energy points we take. The code will change to this a multiple of 10.
    std::vector<double> energy;
    static Parameters from_file();
    dcomp j1; // this is a complex number class defined within the complex library
    int size; //the size of the communicator
    int myid; //the id of the process
    std::vector<int> start; //the starting index of the energy array for each process
    std::vector<int> end; //the ending index of the energy array for each process
    int steps_myid; //this is the number of steps the process has
    std::vector<int> steps_proc; //this is the number of steps the other processes have
    //MPI_Comm comm;
    bool print_gf;
    bool spin_polarised;
    double delta_energy;
    int kk_se; //if true then the code uses the kramer-kronig relation for the self energy
    double e_upper_bound_tilde;       // this is the max energy value
    double e_lower_bound_tilde;
    int multiple; // this is the multiple of the wider energy range.
    int multiple_grids;    
};



double fermi_function(double energy, const Parameters &parameters);
void print_parameters(Parameters& parameters);