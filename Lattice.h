#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <string>

using namespace Eigen;

typedef std::complex<double> dcomp;
typedef std::tuple<dcomp, VectorXi, int, int> elem;

class Lattice{
    private:
    // Fundamental variables
    std::string name;
    int dim;
    int quasi_dim;
    std::vector<int> n_cells;
    int n_orbs;
    std::vector<VectorXd> basis_vecs;
    std::vector<VectorXd> orb_locs;
    std::vector<elem> ham_elems;

    bool PBC;  // Use periodic boundary condition

    // Derived variables
    std::vector<VectorXd> recip_basis;
    int n_sites;
    std::vector<int> dim_incr; // Increments of each dimension in idx

    std::vector<MatrixXcd> BZ_ham;

    public:
    Lattice(std::string geometry_path, bool PBC, bool precalc_ham);
    void geometry_from_file(std::string geometry_path);
    void tight_binding_from_file(std::string geometry_path);

    std::string get_name(){return name;}
    int get_dim(){return dim;}
    int get_n_cells(int dim){return n_cells[dim];}
    std::vector<int> get_n_cells(){return n_cells;}
    int get_n_orbs(){return n_orbs;}
    int get_n_sites(){return n_sites;}

    int site_to_idx(std::vector<int> site);
    std::vector<int> idx_to_site(int idx);
    int sum_idx(int idx1, int idx2);

    VectorXd get_BZ_vector(std::vector<int> k);
    VectorXd get_BZ_vector(int idx);
    MatrixXcd get_Hamiltonian(VectorXd k);

    void precalculate_BZ_ham();

    MatrixXcd get_BZ_Hamiltonian(int idx){
        if(BZ_ham.size() != 0)
            return BZ_ham[idx];
        else
            return this->get_Hamiltonian(this->get_BZ_vector(idx));
    }
    MatrixXcd get_BZ_Hamiltonian(std::vector<int> site){
        if(BZ_ham.size() != 0){
            int idx = this->site_to_idx(site);
            return BZ_ham[idx];
        }
        else
            return this->get_Hamiltonian(this->get_BZ_vector(site));
    }




};
#endif
