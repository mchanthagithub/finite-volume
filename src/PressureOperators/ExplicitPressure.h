//
// Created by maytee on 4/30/18.
//

#ifndef PROJECT_EXPLICITPRESSURE_H
#define PROJECT_EXPLICITPRESSURE_H


#include "PressureOperator.h"

class ExplicitPressure : public PressureOperator {
public:
    void calculatePressure(CartesianGrid& grid, Eigen::MatrixXd H);
    void calculatePressureGradient(CartesianGrid& grid, Eigen::MatrixXd H);
    void clearData();
};


#endif //PROJECT_EXPLICITPRESSURE_H
