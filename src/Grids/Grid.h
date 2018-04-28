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


  // Tells if cell has a velocity BC
  Eigen::VectorXi cellHasVelocityBC;

  // Tells if a cell has a pressure BC
  Eigen::VectorXi cellHasPressureBC;

  // Matrix that is (numFaces) x (nDim+1) in size that tells if the face has a BC and what kind\
  // Column 1 is u-velocity BC, column 2 is v-velocity BC, column 3 is pressure BC here
  // 0 = no BC, positive = Dirichlet, negative = Neumann
  Eigen::MatrixXi BCTypes;

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
