//
// Created by maytee on 4/27/18.
//

#ifndef PROJECT_ADVECTIONOPERATORS_H
#define PROJECT_ADVECTIONOPERATORS_H

#include "Eigen/Core"
#include "../Grids/Grid.h"
#include "../Grids/CartesianGrid.h"
#include "../InterpolationOperators/InterpolateUpwind.h"
class AdvectionOperator {
public:
    // Calculate the advection portion of the flux term
    virtual void calculateAdvectionFluxesCartesian(CartesianGrid &grid, InterpolateUpwind& interp) = 0;
    virtual void clearData() = 0;
    // (numFacesPerCell x nDim) in size to account for x and y components
    // Note that flux is usually a scalar, but because momentum is a vector, there is a vector flux
    std::vector<Eigen::MatrixXd> advectionFluxValues;

    // Total cell flux (numCells x nDim) in size
    Eigen::MatrixXd totalCellAdvectionFlux;

    // Matrices are (numCells x numFacesPerCell) in size
    Eigen::MatrixXd uVelInterp;
    Eigen::MatrixXd vVelInterp;
};


#endif //PROJECT_ADVECTIONOPERATORS_H
