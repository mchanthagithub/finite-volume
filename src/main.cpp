// Finite volume solver for 2D incompressible flows
#include <iostream>
#include "Grids/Grid.h"
#include "Grids/CartesianGrid.h"
#include "OutputUtilities/OutputUtilities.h"
#include <Eigen/Core>

int main() {
  std::cout << "Hello, World!" << std::endl;

  // Generate and initialize grid
  CartesianGrid newCartesianGrid(4,5,3.0,2.0);

  // Set BCs on grid
  newCartesianGrid.setBCs();

  // Initialize values on grid
  newCartesianGrid.setInitialValues();

  // Calculate advection terms


  OutputUtilities output;
  output.writeCartesianCellDataToVTU(newCartesianGrid,"/home/maytee/Documents/2.29/finiteVolumeSolver/cellData.vtu");
  output.writeCartesianFaceDataToVTU(newCartesianGrid,"/home/maytee/Documents/2.29/finiteVolumeSolver/faceData.vtu");
  return 0;
}