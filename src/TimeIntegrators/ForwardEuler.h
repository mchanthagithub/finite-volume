//
// Created by maytee on 4/28/18.
//

#ifndef PROJECT_FORWARDEULER_H
#define PROJECT_FORWARDEULER_H


#include "TimeIntegrator.h"

class ForwardEuler : public TimeIntegrator
{
public:
    void integrate(double delT, Grid& grid, Eigen::VectorXd& old_velocity, Eigen::VectorXd& newVelocity,
              Eigen::VectorXd& pressureGradient, Eigen::VectorXd& fluxTerms, double rho);

};


#endif //PROJECT_FORWARDEULER_H
