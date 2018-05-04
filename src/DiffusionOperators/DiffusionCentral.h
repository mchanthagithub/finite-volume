//
// Created by maytee on 4/30/18.
//

#ifndef PROJECT_DIFFUSIONCENTRAL_H
#define PROJECT_DIFFUSIONCENTRAL_H

#include "DiffusionOperator.h"
#include "../InterpolationOperators/InterpolateUpwind.h"
class DiffusionCentral : public DiffusionOperator {
public:
    void calculateDiffusionFluxesCartesian(CartesianGrid &grid,InterpolateUpwind& interp);
    void calculateGradients(CartesianGrid& grid, InterpolateUpwind& interp);
    void clearData();

    void applyNeumannBCs(CartesianGrid& grid);

    void calculateGradients(CartesianGrid& grid);

    std::vector<Eigen::MatrixXd> diffusionCentralFluxValues;
    std::vector<Eigen::MatrixXd> centralGradients;

    // Is a vector numCells long with each element being another vector that is numFacesPerElement long
    // which then hold (nDim x nDim) matrix, which is the velocity gradient at each face
    std::vector<std::vector<Eigen::MatrixXd> >faceVelocityGradients;
};


#endif //PROJECT_DIFFUSIONCENTRAL_H
