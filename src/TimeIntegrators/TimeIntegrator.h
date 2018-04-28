//
// Created by maytee on 4/28/18.
//

#ifndef PROJECT_TIMEINTEGRATOR_H
#define PROJECT_TIMEINTEGRATOR_H

#include "../Grids/Grid.h"
class TimeIntegrator {
public:
    virtual void integrate(double delT, Grid& grid, Eigen::VectorXd& old_velocity, Eigen::VectorXd& newVelocity,
              Eigen::VectorXd& oldPressure, Eigen::VectorXd& newPressure, Eigen::VectorXd& fluxTerms) = 0;
};


#endif //PROJECT_TIMEINTEGRATOR_H
