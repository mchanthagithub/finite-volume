//
// Created by maytee on 4/30/18.
//

#ifndef PROJECT_CALCULUSFUNCTIONS_H
#define PROJECT_CALCULUSFUNCTIONS_H

#include "Eigen/Core"

Eigen::VectorXd divergenceVector(Eigen::MatrixXd inputVect, int Nx, int Ny, int delX, int delY);

// Gradient of a vector of scalar values (like pressure)
Eigen::MatrixXd gradientScalar(Eigen::VectorXd inputVect, int Nx, int Ny, int delX, int delY);

#endif //PROJECT_CALCULUSFUNCTIONS_H
