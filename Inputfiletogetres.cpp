
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "utils.h"
#include "CNDO2.h"



using namespace std;


int main() {

    AO INPUT_ao("H2.txt");

    cout << "Overlap Matrix for H2: " << endl;
    vector<BasisFunction> basis_set = INPUT_ao.basis_set;

    arma::mat S = overlap_matrix(basis_set);
    S.print();
    CNDO CNDO_INPUT;
    // compute the gamma matrix
    cout << "Gamma Matrix for H2: " << endl;
    int natoms = INPUT_ao.get_natoms();
    arma::mat gamma = CNDO_INPUT.computeGammaMatrix(natoms, basis_set);
    gamma.print();

    cout << "Compute Core Hamiltonian Matrix for H: " << endl;
    vector<string> atom_types = INPUT_ao.get_atom_types();
    arma::mat Hcore = CNDO_INPUT.computeCoreHamiltonianMatrix(atom_types, basis_set);
    Hcore.print();

    cout << "Compute Fock Matrix for H: " << endl;
    arma::mat density_matrix = arma::zeros(basis_set.size(), basis_set.size());
    arma::vec Ptotal = arma::zeros(natoms);
    arma::mat Fock = CNDO_INPUT.computeFockMatrix(atom_types, basis_set, S, Hcore, density_matrix, Ptotal);
    Fock.print();

    CNDO_INPUT.updateDensityMatrix(INPUT_ao, "zhigongH2.txt");


    return 0;
}