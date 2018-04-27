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
  void writeCartesianCellDataToVTU(Grid& grid,std::string fileName);
  void writeCartesianFaceDataToVTU(CartesianGrid& grid,std::string fileName);
};


#endif //PROJECT_OUTPUTUTILITIES_H
