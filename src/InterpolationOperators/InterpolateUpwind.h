//
// Created by maytee on 5/1/18.
//

#ifndef PROJECT_INTERPOLATEUPWIND_H
#define PROJECT_INTERPOLATEUPWIND_H

#include "InterpolationOperator.h"
class InterpolateUpwind : public InterpolationOperator {
public:
    void interpolateVelocities(CartesianGrid& grid);
    void clearData();
    void applyDirichletBCs(CartesianGrid& grid);

};


#endif //PROJECT_INTERPOLATEUPWIND_H
