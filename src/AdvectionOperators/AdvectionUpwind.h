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

    // Calculate cell-averaged velocities at element centers
    void calculateInterpolationValues(CartesianGrid& grid);

    // Apply Dirichlet BCs to interpolated face values; should be called in calculateInterpolationValues
    void applyDirichletBCs(CartesianGrid &grid);

    void clearData();
    // Apply Neumann BCs
    //void applyNeumannBCs(CartesianGrid& grid);
};


#endif //PROJECT_ADVECTIONUPWIND_H
