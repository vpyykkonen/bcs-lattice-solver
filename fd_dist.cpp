#include "fd_dist.h"

#include <Eigen/Dense>
#include <cmath>


double fd_dist(double E, double mu,  double T)
{
    double dist;
    double exp_term;

    if (std::abs(T)<1.0e-13){
        if ( E - mu <= 0.0)
            dist = 1.0;
        else
            dist = 0.0;
        return dist;
    }
    exp_term = std::exp((E-mu)/T);

    if (std::isnan(exp_term)){
        dist = 0.0;
        return dist;
    }
    dist = 1.0/(1.0 + exp_term);
    return dist;
}

Eigen::VectorXd fd_dist(Eigen::VectorXd Es, double mu, double  T)
{
    Eigen::VectorXd fd_Es = Eigen::VectorXd::Zero(Es.size());
    for(int n = 0; n < Es.size(); n++)
        fd_Es(n) = fd_dist(Es(n),mu,T);
    return fd_Es;

}

double fd_deriv(double E, double mu,  double T)
{
    double dist_deriv;

    if ( T == 0 ){
        if(std::abs(E) < 1e-10 || E == 0)
            dist_deriv = 1.0;
        else
            dist_deriv = 0.0;
        return dist_deriv;
    }

    dist_deriv = - (1.0/T)/(2.0*std::pow(cosh((E-mu)/(2*T)),2));
    return dist_deriv;
}
