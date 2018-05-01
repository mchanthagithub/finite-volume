//
// Created by maytee on 4/29/18.
//

#ifndef PROJECT_DIFFUSIONOPERATOR_H
#define PROJECT_DIFFUSIONOPERATOR_H

#include "../Grids/CartesianGrid.h"
class DiffusionOperator {
public:
  virtual void calculateDiffusionFluxesCartesian(CartesianGrid &grid) = 0;
  virtual void clearData()=0;

  // (numFacesPerCell x nDim) in size to account for x and y components
  // Note that flux is usually a scalar, but because momentum is a vector, there is a vector flux
  std::vector<Eigen::MatrixXd> diffusionFluxValues;

  // Total cell flux (numCells x nDim) in size
  Eigen::MatrixXd totalCellDiffusionFlux;

};


#endif //PROJECT_DIFFUSIONOPERATOR_H
