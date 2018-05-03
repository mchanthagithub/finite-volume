//
// Created by maytee on 5/1/18.
//

#ifndef PROJECT_INTERPOLATEQUICK_H
#define PROJECT_INTERPOLATEQUICK_H

#include "InterpolationOperator.h"
class InterpolateQUICK : public InterpolationOperator {
public:
    void interpolateVelocities(CartesianGrid& grid);
    void clearData();
    void applyDirichletBCs(CartesianGrid& grid);
};


#endif //PROJECT_INTERPOLATEQUICK_H
