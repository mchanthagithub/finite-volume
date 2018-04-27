// Finite volume solver for 2D incompressible flows
#include <iostream>
#include "Grids/Grid.h"
#include "Grids/CartesianGrid.h"
#include <Eigen/Core>

int main() {
  std::cout << "Hello, World!" << std::endl;

  // Generate and initialize grid
  CartesianGrid newCartesianGrid(3,4,3.0,2.0);

  // Set BCs on grid
  newCartesianGrid.setBCs();

  return 0;
}