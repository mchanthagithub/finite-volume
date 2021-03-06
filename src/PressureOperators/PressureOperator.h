//
// Created by maytee on 4/30/18.
//

#ifndef PROJECT_PRESSUREOPERATOR_H
#define PROJECT_PRESSUREOPERATOR_H

#include "Eigen/Core"
#include "../Grids/CartesianGrid.h"

class PressureOperator {
public:
    // Note that H is the sum of the advective and diffusive flux
    virtual void calculatePressure(CartesianGrid& grid, Eigen::MatrixXd H) = 0;
    virtual void calculatePressureGradient(CartesianGrid& grid, Eigen::MatrixXd H) = 0;
    virtual void clearData(CartesianGrid& grid) = 0;


    Eigen::VectorXd pressure;
    Eigen::MatrixXd pressureGradient;
};


#endif //PROJECT_PRESSUREOPERATOR_H
