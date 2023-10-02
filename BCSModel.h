#ifndef BCSMODEL_H
#define BCSMODEL_H

#include <vector>
#include <complex>
#include <Eigen/Dense>

#include "Lattice.h"

typedef std::complex<double> dcomp;
typedef std::tuple<dcomp, VectorXi, int, int> elem;

using namespace std;
using namespace Eigen;

class BCSModel{
    private:
    Lattice lattice;
    double U;
    double T;
    double mu;
    VectorXd Hartree;
    VectorXcd Delta;

    public:
    BCSModel(Lattice lattice, double U, double T,double mu);
    MatrixXcd get_HBdG(VectorXd k);
    MatrixXcd get_HBdG(int idx);
    MatrixXcd get_HBdG(int idx, VectorXd new_Hartree, VectorXcd new_Delta);
    VectorXcd solve_Hartree_and_Delta(const VectorXcd& X);
    VectorXcd solve_Delta(const VectorXcd& X);
    //VectorXcd solve_Hartree_and_Delta();
    VectorXd get_Hartree(){return Hartree;}
    VectorXcd get_Delta(){return Delta;}
    void self_consistent_loop(VectorXd Hartree0, VectorXcd Delta0, string scf_cfg_path,string save_path);
    void self_consistent_loop(VectorXcd Delta0, string scf_cfg_path,string save_path);


    MatrixXcd get_Gl();

    void set_U(double U){this->U = U;}
    void set_mu(double mu){this->mu = mu;}
    void set_T(double T){this->T = T;}
    void set_Hartree(VectorXd Hartree){this->Hartree = Hartree;}
    void set_Delta(VectorXcd Delta){this->Delta = Delta;}
};

#endif
