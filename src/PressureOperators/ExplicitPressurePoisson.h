//
// Created by maytee on 5/2/18.
//

#ifndef PROJECT_EXPLICTPRESSUREPOISSON_H
#define PROJECT_EXPLICTPRESSUREPOISSON_H

#include "PressureOperator.h"
class ExplicitPressurePoisson : public PressureOperator {
public:
    // Note that H is the sum of the advective and diffusive flux
    void calculatePressure(CartesianGrid& grid, Eigen::MatrixXd H);
    void calculatePressureGradient(CartesianGrid& grid, Eigen::MatrixXd H);
    void clearData();
    void createMappings(CartesianGrid& grid);

    Eigen::VectorXd pressure;
    Eigen::MatrixXd pressureGradient;

    Eigen::VectorXi mappingGlobalToActive;
    Eigen::VectorXi mappingActiveToGlobal;
    int totalDOF;
};


#endif //PROJECT_EXPLICTPRESSUREPOISSON_H
