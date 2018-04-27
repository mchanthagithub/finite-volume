//
// Created by maytee on 4/27/18.
//

#ifndef PROJECT_ADVECTIONOPERATORS_H
#define PROJECT_ADVECTIONOPERATORS_H

#include "Eigen/Core"
#include "../Grids/Grid.h"
#include "../Grids/CartesianGrid.h"

class AdvectionOperator {
public:
    // Calculate the advection portion of the flux term
    virtual void calculateAdvectionFluxesCartesian(CartesianGrid &grid) = 0;

    // Vectors below are as long as the number of faces
    Eigen::VectorXd advectionFluxValues;
    Eigen::VectorXd uVelInterp;
    Eigen::VectorXd vVelInterp;
};


#endif //PROJECT_ADVECTIONOPERATORS_H
