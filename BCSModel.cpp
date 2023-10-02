#include <iostream>
#include <complex>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "BCSModel.h"
#include "Lattice.h"
#include "fd_dist.h"
#include "ScfMethod.h"
#include "ScfSolver.h"

using namespace Eigen;
using namespace std;

BCSModel::BCSModel(Lattice lattice, double U, double T, double mu):lattice{lattice},U{U},T{T},mu{mu}
{
    Hartree = VectorXd::Zero(lattice.get_n_orbs());
    Delta = VectorXcd::Zero(lattice.get_n_orbs());
}

MatrixXcd BCSModel::get_HBdG(VectorXd k)
{
    int n_orbs = lattice.get_n_orbs();
    MatrixXcd HBdG = MatrixXcd::Zero(2*n_orbs,2*n_orbs);
    MatrixXcd H0_up = lattice.get_Hamiltonian(k);
    MatrixXcd H0_down = lattice.get_Hamiltonian(-k);
    for(int i = 0; i < n_orbs; i++){
        for(int j = 0; j < n_orbs; j++){
            if( i == j){
                HBdG(2*i,2*i) += Hartree(i)-mu;
                HBdG(2*i,2*i+1) += Delta(i);
                HBdG(2*i+1,2*i) += conj(Delta(i));
                HBdG(2*i+1,2*i+1) -= Hartree(i)-mu;
            }
            HBdG(2*i,2*j) += H0_up(i,j);
            HBdG(2*i+1,2*j+1) -= H0_down(j,i);
        }
    }
    return HBdG;
}

MatrixXcd BCSModel::get_HBdG(int idx)
{
    int n_orbs = lattice.get_n_orbs();
    MatrixXcd HBdG = MatrixXcd::Zero(2*n_orbs,2*n_orbs);

    vector<int> site_up = lattice.idx_to_site(idx);
    vector<int> site_down;
    for(int i = 0; i < site_up.size();++i)
        site_down.push_back(-site_up[i]);
    int idx_down = lattice.site_to_idx(site_down);

    MatrixXcd H0_up = lattice.get_BZ_Hamiltonian(idx);
    MatrixXcd H0_down = lattice.get_BZ_Hamiltonian(idx_down);

    for(int i = 0; i < n_orbs; i++){
        for(int j = 0; j < n_orbs; j++){
            if( i == j){
                HBdG(2*i,2*i) += Hartree(i)-mu;
                HBdG(2*i,2*i+1) += Delta(i);
                HBdG(2*i+1,2*i) += conj(Delta(i));
                HBdG(2*i+1,2*i+1) -= Hartree(i)-mu;
            }
            HBdG(2*i,2*j) += H0_up(i,j);
            HBdG(2*i+1,2*j+1) -= H0_down(j,i);
        }
    }
    return HBdG;
}

MatrixXcd BCSModel::get_HBdG(int idx, VectorXd new_Hartree,
        VectorXcd new_Delta)
{
    int n_orbs = lattice.get_n_orbs();
    MatrixXcd HBdG = MatrixXcd::Zero(2*n_orbs,2*n_orbs);

    vector<int> site_up = lattice.idx_to_site(idx);
    vector<int> site_down;
    for(int i = 0; i < site_up.size();++i)
        site_down.push_back(-site_up[i]);
    int idx_down = lattice.site_to_idx(site_down);

    MatrixXcd H0_up = lattice.get_BZ_Hamiltonian(idx);
    MatrixXcd H0_down = lattice.get_BZ_Hamiltonian(idx_down);
    for(int i = 0; i < n_orbs; i++){
        for(int j = 0; j < n_orbs; j++){
            if( i == j){
                HBdG(2*i,2*i) += new_Hartree(i)-mu;
                HBdG(2*i,2*i+1) += new_Delta(i);
                HBdG(2*i+1,2*i) += conj(new_Delta(i));
                HBdG(2*i+1,2*i+1) -= new_Hartree(i)-mu;
            }
            HBdG(2*i,2*j) += H0_up(i,j);
            HBdG(2*i+1,2*j+1) -= H0_down(j,i);
        }
    }
    return HBdG;
}

VectorXcd BCSModel::solve_Hartree_and_Delta(const VectorXcd& X)
{
    int n_orbs = lattice.get_n_orbs();
    vector<int> n_cells = lattice.get_n_cells();
    int n_sites = lattice.get_n_sites();

    VectorXd new_Hartree = X.head(n_orbs).real();
    VectorXcd new_Delta = X.tail(n_orbs);

    VectorXcd output = VectorXcd::Zero(2*n_orbs);
    if( abs(U) < 1e-6) // If interaction is too small, return zero
        return output;

    for(int n = 0; n < n_sites ;n++){
        MatrixXcd HBdG = this->get_HBdG(n,new_Hartree,new_Delta);
        SelfAdjointEigenSolver<MatrixXcd> es;
        es.compute(HBdG,ComputeEigenvectors);
        MatrixXcd W = es.eigenvectors();
        VectorXd eigE = es.eigenvalues();
        VectorXd fd_eigE = fd_dist(eigE,0.0,T);

        MatrixXcd Gl = 1.0i*W*fd_eigE.asDiagonal()*W.adjoint();
        for(int m = 0; m < n_orbs;m++){
            output(m) += -1.0i*Gl(2*m,2*m);
            output(m+n_orbs) += -1.0i*Gl(2*m,2*m+1);
        }

        //cout << "n" << "\n";
        //cout << n << "\n";
        //cout << "H0" << "\n";
        //cout << lattice.get_BZ_Hamiltonian(n) << "\n";
        //cout << "HBdG: " << "\n";
        //cout << HBdG << "\n";
        //cout << "W: " << "\n";
        //cout << W << "\n";
        //cout << "eigE: " << "\n";
        //cout << eigE << "\n";
        //cout << "fd_eigE: " << "\n";
        //cout << fd_eigE << "\n";

        //for(int m = 0; m < n_orbs; m++){
        //    VectorXcd fd_eigE_Wadj = fd_eigE.cwiseProduct(W.row(2*m).adjoint());
        //    //cout << "fd_eigE_Wadj" << "\n";
        //    //cout << fd_eigE_Wadj << "\n";
        //    output(m) += real(W.row(2*m).adjoint().dot(fd_eigE_Wadj)); // Hartree from up,up block
        //    output(m+n_orbs) += conj(W.row(2*m+1).adjoint().dot(fd_eigE_Wadj)); // Delta from down,up block
        //}
    }
    output *= U/n_sites;

    return output;
}

VectorXcd BCSModel::solve_Delta(const VectorXcd& X)
{
    int n_orbs = lattice.get_n_orbs();
    vector<int> n_cells = lattice.get_n_cells();
    int n_sites = lattice.get_n_sites();

    VectorXd new_Hartree = VectorXd::Zero(n_orbs);

    VectorXcd output = VectorXcd::Zero(n_orbs);
    if( abs(U) < 1e-6) // If interaction is too small, return zero
        return output;

    for(int n = 0; n < n_sites ;n++){
        MatrixXcd HBdG = this->get_HBdG(n,new_Hartree,X);
        SelfAdjointEigenSolver<MatrixXcd> es;
        es.compute(HBdG,ComputeEigenvectors);
        MatrixXcd W = es.eigenvectors();
        VectorXd eigE = es.eigenvalues();
        VectorXd fd_eigE = fd_dist(eigE,0.0,T);

        MatrixXcd Gl_k = 1.0i*W*fd_eigE.asDiagonal()*W.adjoint();
        for(int m = 0; m < n_orbs;m++)
            output(m) += -1.0i*Gl_k(2*m,2*m+1);

        //for(int m = 0; m < n_orbs; m++){
        //    VectorXcd fd_eigE_Wadj = fd_eigE.cwiseProduct(W.row(2*m).adjoint());
        //    output(m) += conj(W.row(2*m+1).adjoint().dot(fd_eigE_Wadj)); // Delta from down,up block
        //}
    }
    output *= U/n_sites;

    return output;
}


void BCSModel::self_consistent_loop(
        VectorXd Hartree0, VectorXcd Delta0, string scf_cfg_path,string save_path){
    int n_orbs = lattice.get_n_orbs();
    using namespace placeholders;

    VectorXcd X(2*n_orbs);

    X.head(n_orbs) = Hartree0;
    X.tail(n_orbs) = Delta0;
    ScfSolver solver(&X,scf_cfg_path,"");
    solver.iterate(std::bind(
                    static_cast<VectorXcd(BCSModel::*)
                    (const VectorXcd&)>
                    (&BCSModel::solve_Hartree_and_Delta),
                    this,_1));

    int total_iter = solver.get_iterations();
    bool converged = solver.get_converged();

    cout << "Scf loop converged: " << converged << "\n";
    cout << "Total iterations: " << total_iter << "\n";

    this->Hartree = X.head(n_orbs).real();
    this->Delta = X.tail(n_orbs);
}

void BCSModel::self_consistent_loop(
        VectorXcd Delta0, string scf_cfg_path,string save_path){
    int n_orbs = lattice.get_n_orbs();
    using namespace placeholders;

    VectorXcd X(n_orbs);
    X = Delta0;
    ScfSolver solver(&X,scf_cfg_path,"");
    solver.iterate(std::bind(
                    static_cast<VectorXcd(BCSModel::*)
                    (const VectorXcd&)>
                    (&BCSModel::solve_Delta),
                    this,_1));

    int total_iter = solver.get_iterations();
    bool converged = solver.get_converged();

    cout << "Scf loop converged: " << converged << "\n";
    cout << "Total iterations: " << total_iter << "\n";

    this->Delta = X;
}

MatrixXcd BCSModel::get_Gl()
{
    int n_orbs = lattice.get_n_orbs();
    int n_sites = lattice.get_n_sites();
    MatrixXcd Gl = MatrixXcd::Zero(2*n_orbs,2*n_orbs);
    for(int n = 0; n < n_sites ;n++){
        MatrixXcd HBdG = this->get_HBdG(n);
        SelfAdjointEigenSolver<MatrixXcd> es;
        es.compute(HBdG,ComputeEigenvectors);
        MatrixXcd W = es.eigenvectors();
        VectorXd eigE = es.eigenvalues();
        VectorXd fd_eigE = fd_dist(eigE,0.0,T);
        Gl += 1.0i*W*fd_eigE.asDiagonal()*W.adjoint();
    }
    Gl /= n_sites;
    return Gl;
}

