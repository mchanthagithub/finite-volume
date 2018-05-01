//
// Created by maytee on 4/30/18.
//

#ifndef PROJECT_INTERPOLATIONOPERATOR_H
#define PROJECT_INTERPOLATIONOPERATOR_H

#include "Eigen/Core"
#include "../Grids/CartesianGrid.h"
class InterpolationOperator {
public:
    virtual void interpolateVelocities(CartesianGrid& grid) = 0;
    virtual void clearData() = 0;
    Eigen::MatrixXd uVelInterp;
    Eigen::MatrixXd vVelInterp;
};


#endif //PROJECT_INTERPOLATIONOPERATOR_H
