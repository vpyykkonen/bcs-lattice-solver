#ifndef FD_DIST_H
#define FD_DIST_H

#include <Eigen/Dense>

double fd_dist(double E, double mu, double T);
Eigen::VectorXd fd_dist(Eigen::VectorXd Es, double mu, double T);

double fd_deriv(double E, double mu, double T);

#endif
