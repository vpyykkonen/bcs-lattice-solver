#include <vector>
#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <math.h>

#include "config_parser.h"
#include "Lattice.h"

using namespace Eigen;
using namespace std;


typedef std::complex<double> dcomp;
typedef std::tuple<dcomp, VectorXi, int, int> elem;

Lattice::Lattice(std::string geometry_path, bool PBC, bool precalc_ham)
{
    geometry_from_file(geometry_path);
    tight_binding_from_file(geometry_path);
    if(precalc_ham)
        precalculate_BZ_ham();
}

void Lattice::geometry_from_file(std::string geometry_path)
{
    string tmp,tmp1,tmp2;
    int n = 0;
    map<string,string> geometry_config = load_config_file(geometry_path, "=");

    name = geometry_config["lattice_name"];
    dim = stoi(geometry_config["dim"]);
    quasi_dim = stoi(geometry_config["quasi_dim"]);
    n_orbs = stoi(geometry_config["n_orbs"]);
    
    // Basis vectors
    for(int i = 0; i < dim; i++){ 
        tmp = geometry_config["basis"+to_string(i+1)];
        VectorXd basis_vec = VectorXd::Zero(dim);
        n = 0;
        while(tmp.size() != 0){
            size_t sep_pos = tmp.find(",");
            basis_vec(n) = stod(tmp.substr(0,sep_pos));
            if(sep_pos == string::npos){
                tmp.erase();
            }
            else
                tmp = tmp.erase(0,sep_pos+1);
            n++;
        }
        basis_vecs.push_back(basis_vec);
    }
    // Orbital locations within the unit cell
    for(int i = 0; i < n_orbs; i++){ 
        tmp = geometry_config["orb_loc"+to_string(i+1)];
        VectorXd orb_loc_vec = VectorXd::Zero(quasi_dim);
        n = 0;
        while(tmp.size() > 0){
            size_t sep_pos = tmp.find(",");
            orb_loc_vec(n) = stod(tmp.substr(0,sep_pos));
            if(sep_pos == string::npos)
                tmp = tmp.erase();
            else
                tmp = tmp.erase(0,sep_pos+1);
            n++;
        }
        orb_locs.push_back(orb_loc_vec);
    }

    // To do: impelement closed system
    PBC = true;

    // Reciprocal lattice basis vectors
    for(n = 0; n < dim; n++){

        VectorXd b = VectorXd::Zero(dim);
        if( dim == 1){
            b = basis_vecs[n];
        }
        if( dim == 2){
            Vector2d vec = basis_vecs[0];
            if( n == 0)
                vec = basis_vecs[1];
            b(0) = vec(1);
            b(1) = -vec(0);
        }
        if( dim == 3){
            Vector3d vec1 = basis_vecs[1];
            Vector3d vec2 = basis_vecs[2];
            if(n == 0){
                vec1 = basis_vecs[1];
                vec2 = basis_vecs[2];
            }
            if(n == 1){
                vec1 = basis_vecs[2];
                vec2 = basis_vecs[0];
            }
            if(n == 2){
                vec1 = basis_vecs[0];
                vec2 = basis_vecs[1];
            }
            b = vec1.cross(vec2);
        }
        b *= 2*M_PI/basis_vecs[n].norm();
        recip_basis.push_back(b);
    }
    // Fix handedness of basis
    //if(dim == 2){
    //    if(recip_basis[0](0)*recip_basis[1](1)-recip_basis[0](1)*recip_basis[1](0) < 0) 
    //        recip_basis[1] *= -1.0;
    //}
    //if(dim == 3){
    //    if(recip_basis[0].dot(recip_basis[1].cross(recip_basis[2])) < 0)
    //        recip_basis[3] *= -1.0;
    //}

    // Size of the lattice, number of unit cells in each dimension
    tmp = geometry_config["n_cells"];
    while(tmp.length() > 0){
        size_t pos = tmp.find(",");
        n_cells.push_back(stoi(tmp.substr(0,pos)));
        if(pos == string::npos)
            tmp = tmp.erase();
        else
            tmp = tmp.erase(0,pos+1);

    }

    dim_incr.push_back(1);
    for(n = 0; n < dim; n++)
        dim_incr.push_back(n_cells[n]*dim_incr[n]);
    n_sites = dim_incr[dim];
}

void Lattice::tight_binding_from_file(std::string geometry_path)
{
    map<string,string> geometry_config = load_config_file(geometry_path, "=");

    string tmp,tmp1,tmp2;
    int n = 0;

    // Tight-binding parameters
    // on-site
    tmp = geometry_config["on_site"];
    n = 0;
    while(tmp.length() > 0){
        size_t pos = tmp.find(",");
        dcomp on_site;
        istringstream strs(tmp.substr(0,pos));
        strs >> on_site;
        elem on_site_term = std::make_tuple(on_site,VectorXi::Zero(dim),n,n);
        n++;
        ham_elems.push_back(on_site_term);
        if(pos == string::npos)
            tmp = tmp.erase();
        else
            tmp = tmp.erase(0,pos+1);
    }


    // hoppings
    // format: 
    // target_vector;target_orb;source_orb = amplitude
    int n_hoppings = stoi(geometry_config["n_hoppings"]);
    for(int i = 0; i < n_hoppings; i++){
        tmp = geometry_config["hopping"+to_string(i+1)]; // store whole str
        // split source+target from amplitude, separated by "="
        size_t pos_eq = tmp.find("=");
        size_t pos;
        // handle first the source+target vector
        // format: 
        tmp1 = tmp.substr(0,pos_eq); // store target+orbital
        dcomp amp;
        int target_orb;
        int source_orb;
        VectorXi target_cell = VectorXi::Zero(dim);
        // first part the target location vector
        pos = tmp1.find(";");
        tmp2 = tmp1.substr(0,pos); // store position vector
        n = 0;
        while(tmp2.length() > 0){
            size_t pos2 = tmp2.find(",");
            target_cell(n) = stoi(tmp2.substr(0,pos2));
            if(pos2 == string::npos)
                tmp2 = tmp2.erase();
            else
                tmp2 = tmp2.erase(0,pos2+1);
            n++;
        }
        tmp1.erase(tmp1.begin(),tmp1.begin()+pos+1);
        // Target orbital index
        pos = tmp1.find(";");
        target_orb = stoi(tmp1.substr(0,pos));
        tmp1.erase(tmp1.begin(),tmp1.begin()+pos+1); 
        // Source orbital index
        source_orb = stoi(tmp1);
        
        // Amplitude
        tmp2 = tmp.substr(pos_eq+1,tmp.length());
        istringstream strs(tmp2);
        strs >> amp;

        cout << amp << " " << target_cell << " " << target_orb << " " << source_orb << "\n";

        elem hopping = std::make_tuple(amp,target_cell,target_orb,source_orb);
        ham_elems.push_back(hopping);
    }
}

// site to 1BZ index
int Lattice::site_to_idx(vector<int> site){
    int output = 0;
    for(int n = 0; n < dim; n++){
        int BZ_index = site[n]%n_cells[n];
        if(BZ_index < 0)
            BZ_index += n_cells[n];
        output += dim_incr[n]*BZ_index;
    }
    return output;
}

// 1BZ index to site
vector<int> Lattice::idx_to_site(int idx){
    vector<int> output;
    if(idx >= n_sites)
        idx %= n_sites;
    for(int n = 0; n < dim;n++)
        output.push_back((idx%dim_incr[n+1])/dim_incr[n]);
    return output;
}

int Lattice::sum_idx(int idx1, int idx2)
{
    std::vector<int> site1 = this->idx_to_site(idx1);
    std::vector<int> site2 = this->idx_to_site(idx2);
    std::vector<int> sum;
    for(int n = 0; n< dim; n++)
        sum.push_back(site1[n]+site2[n]);
    return site_to_idx(sum);
}

VectorXd Lattice::get_BZ_vector(std::vector<int> k_int)
{
    VectorXd k = VectorXd::Zero(dim);
    for(int n = 0; n < dim; n++){
        int n_idx = k_int[n]%n_cells[n]; // Project to 1 BZ
        k += double(n_idx)*recip_basis[n]/double(n_cells[n]);
    }
    return k;
}

VectorXd Lattice::get_BZ_vector(int idx){
    VectorXd k = VectorXd::Zero(dim);
    std::vector<int> site = this->idx_to_site(idx);
    return get_BZ_vector(site);
}


MatrixXcd Lattice::get_Hamiltonian(VectorXd k)
{
    MatrixXcd output = MatrixXcd::Zero(n_orbs,n_orbs);
    int n_elems = ham_elems.size();
    for(int n = 0; n < n_elems; n++){
        VectorXi target_pos = std::get<1>(ham_elems[n]);
        int target_orb = std::get<2>(ham_elems[n]);
        int source_orb = std::get<3>(ham_elems[n]);
        dcomp amp = std::get<0>(ham_elems[n]);
        // k-space phase term
        for(int m =  0; m < dim; m++)
            amp *= std::exp(-1.0i*k.dot(target_pos(m)*basis_vecs[m]));
        // Orbital position effect on the amplitude phase
        amp *= std::exp(-1.0i*k.dot(orb_locs[target_orb]-orb_locs[source_orb]));
        output(target_orb,source_orb) += amp;
        // Add complex conjugate as well
        output(source_orb,target_orb) += conj(amp);
    }
    return output;
}



void Lattice::precalculate_BZ_ham(){
    for(int n = 0; n < n_sites; n++){
        VectorXd k = this->get_BZ_vector(n);
        BZ_ham.push_back(this->get_Hamiltonian(k));
    }
}




