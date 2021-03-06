//
// Created by maytee on 4/27/18.
//

#ifndef PROJECT_OUTPUTUTILITIES_H
#define PROJECT_OUTPUTUTILITIES_H

#include "../Grids/Grid.h"
#include "../Grids/CartesianGrid.h"
#include <iostream>
#include <fstream>
class OutputUtilities {
public:
  void writeCartesianCellDataToVTU(CartesianGrid& grid,std::string fileName);
  void writeCartesianFaceDataToVTU(CartesianGrid& grid,std::string fileName);
  void writePlottingCartesianDataToVTU(CartesianGrid& grid, std::string fileName);
  std::string writeCellScalar(Eigen::VectorXd inputVector, int nDim, std::string inputString);
  std::string writeCellScalar(Eigen::VectorXi inputVector, int nDim, std::string inputString);

  std::string writeCellVector(Eigen::VectorXd inputVector, int nDim, std::string inputString);
  std::string writeCellVector(Eigen::VectorXi inputVector, int nDim, std::string inputString);
};


#endif //PROJECT_OUTPUTUTILITIES_H
