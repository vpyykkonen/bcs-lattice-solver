#include <iostream>
#include <complex>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>




#include "file_io.h"
#include "Lattice.h"
#include "BCSModel.h"
#include "H5Cpp.h"
using namespace std;
using namespace Eigen;


int main(void){
    string scf_path("./scf_params.cfg");
    Lattice lattice("./geometry.cfg",true,true);

    int n_sites = lattice.get_n_sites();
    int dim = lattice.get_dim();
    int n_orbs = lattice.get_n_orbs();

    double U = -1.0;
    double mu = -2.0;
    double T = 0.0;
    int n_mus = 50;
    double mu_min = 1.9;
    double mu_max = 2.1;
    string note = "flat_band_zoom";
    VectorXd mus = VectorXd::LinSpaced(n_mus,mu_min,mu_max);


    BCSModel bcs(lattice,U,T,0.0);
    VectorXd Hartree0 = 0.4*U*VectorXd::Constant(n_orbs,1.0);
    Hartree0(1) = 0.3*U;
    //Hartree0(2) = 0.3*U;
    VectorXcd Delta0 = 0.1*U*VectorXcd::Constant(n_orbs,1.0);

    VectorXd n_particles = VectorXd::Zero(n_mus);
    VectorXd Delta_max = VectorXd::Zero(n_mus);

    string output_name = lattice.get_name()+"_";
    for(int n = 0; n < dim; n++){
        output_name += to_string(lattice.get_n_cells(n));
        if(n != dim-1)
            output_name += string("x");
    }
    output_name += "_"+note+"_U"+to_string(U)+"_T"+to_string(T)+".h5";

    vector<MatrixXcd> Gls;

    for(int n = 0; n < n_mus; n++){
        bcs.set_mu(mus(n));
        bcs.self_consistent_loop(
                Hartree0,Delta0,scf_path,"");
        //bcs.self_consistent_loop(
        //        Delta0,scf_path,"");
        
        Gls.push_back(bcs.get_Gl());
        //for(int m = 0; m < n_orbs; m++)
        //    n_particles(n) += real(-1.0i*Gl(2*m,2*m));
        //Delta_max(n) =bcs.get_Delta().cwiseAbs().maxCoeff();
    }

    H5File* save_file = new H5File("./Test/"+output_name,H5F_ACC_TRUNC);
    save_data_to_h5attr(save_file,"U","double",&U);
    save_data_to_h5attr(save_file,"T","double",&T);
    write_MatrixXd_to_group(save_file,mus,"mus");
    for(int n = 0; n < n_mus;n++)
        write_MatrixXcd_to_group(save_file,Gls[n],string("Gl")+to_string(n+1));
    delete save_file;


    //MatrixXcd Gl = bcs.get_Gl();
    //cout << Gl << endl;
    //VectorXcd X = VectorXcd::Zero(2*n_orbs);
    //X.head(n_orbs)=Hartree0;
    //X.tail(n_orbs)=Delta0;
    //VectorXcd FX = bcs.solve_Hartree_and_Delta(X);
    //cout << FX << "\n";

    //bcs.self_consistent_loop(Hartree0,Delta0,scf_path,"");
    //cout << bcs.get_Hartree() << "\n";
    //cout << bcs.get_Delta() << "\n";





    //for(int n = 0; n < n_sites; n++){
    //    //for(int m = 0; m < n_cells2; m++){
    //    MatrixXcd ham_k = lattice.get_BZ_Hamiltonian(n);
    //    vector<int> site = lattice.idx_to_site(n);
    //    cout << "Idx: " << n << "\n";
    //    cout << "Site: ";
    //    for(int m = 0; m < dim; m++){
    //        cout << site[m] ;
    //        if(m != dim-1)
    //            cout << ", ";
    //    }
    //    cout << "\n";
    //    cout << "ham_k: " << "\n";
    //    cout << ham_k.cwiseAbs() << "\n";

    //    SelfAdjointEigenSolver<MatrixXcd> es;
    //    es.compute(ham_k,ComputeEigenvectors);
    //    MatrixXcd W = es.eigenvectors();
    //    VectorXd eigE = es.eigenvalues();
    //    cout << n << ", \n";
    //    cout << eigE << "\n";
    //}
}
