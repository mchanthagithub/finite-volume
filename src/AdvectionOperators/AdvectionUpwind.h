//
// Created by maytee on 4/27/18.
//
// Applies upwind scheme to determine advection fluxes

#ifndef PROJECT_ADVECTIONUPWIND_H
#define PROJECT_ADVECTIONUPWIND_H


#include "AdvectionOperator.h"

class AdvectionUpwind : public AdvectionOperator{
public:
    void calculateAdvectionFluxesCartesian(CartesianGrid &grid);
    void calculateInterpolationValues(CartesianGrid& grid);
};


#endif //PROJECT_ADVECTIONUPWIND_H
