//
// Created by maytee on 4/26/18.
//
// Grid base class

#ifndef PROJECT_GRID_H
#define PROJECT_GRID_H
#include <Eigen/Core>
#include <vector>

class Grid {
public:
    virtual void setBCs() = 0;
    virtual void setInitialValues() = 0;

  // Cell center coordinates in x1,y1,x2,y2,...,xn,yn format
  Eigen::VectorXd centerXY;

  // Cell corner coordinates in x1,y1,...,xn,yn format
  Eigen::VectorXd cornerXY;

  // Values of cell-averaged velocity and pressure at cell centers
  Eigen::VectorXd uVel;
  Eigen::VectorXd vVel;
  Eigen::VectorXd pressure;

  // Hold the value of BCs
  // Are as long as the total number of different BCs for that field, ex. if there are BCs for uvel = 0 m/s, 2 m/s,
  // and 4 m/s, then uVelBCValue will be 3 elements long
  Eigen::VectorXd uVelBCValue;
  Eigen::VectorXd vVelBCValue;
  Eigen::VectorXd pressureBCValue;

  // Tells if cell has a BC and its type
  // Vector that is numCells in length, with each vector element being a matrix that is (1+numElementFaces)x3 in size
  // First element of each row and column is either 0 or 1: 0 for no BCs at all, 1 for BCs
  // 2nd through numElementFaces+1 element in 1st column is the type of u BC for that face: 0 = no BC, positive = Dirichlet, negative = Neumann
  //  "                       "            in 2nd column is type of v BC, and 3rd column is type of pressure BC\
  // For example a value of 2 in the first column second row would mean the south cell face has a dirichlet condition
  // listed second in the uVelBCValue vector
  std::vector<Eigen::MatrixXi> cellBCTypes;

  // Maps the local corner node IDs to global node IDs
  // Size of (number of nodes per cell) x number of cells
  Eigen::MatrixXi cornerMap;
  // Maps the local side ID to global side ID
  // Size of (number of faces per cell) x number of cells
  Eigen::MatrixXi faceMap;

  // Global number of cells, faces and corners
  int numCells;
  int numFaces;
  int numCorners;

  // Right now this assumes quads
  int numFacesPerElement = 4;
  int numCornersPerElement = 4;
  int nDim = 2;
};


#endif //PROJECT_GRID_H
