//
// Created by maytee on 4/26/18.
//

#include "CartesianGrid.h"
#include<iostream>
// Nx = number of cells in X, Ny = number of cells in Y, Lx = length in X, Ly = length in Y
CartesianGrid::CartesianGrid(int in_Nx, int in_Ny, double in_Lx, double in_Ly)
{
  Nx = in_Nx;
  Ny = in_Ny;
  Lx = in_Lx;
  Ly = in_Ly;
  delX = Lx/Nx;
  delY = Ly/Ny;
  std::cout<<Nx<<" "<<Ny<<" "<<Lx<<" "<<Ly<<" "<<delX<<" "<<delY<<std::endl;

  int nFaces = Nx*(Ny+1) + Ny*(Nx+1);
  int nElem = Nx*Ny;
  int nCorners = (Nx+1)*(Ny+1);

  setVariableSizes(nElem,nFaces,nCorners);
  generateMesh();
  setBCs();
  std::cout<<"Centers:"<<std::endl;
  std::cout<<centerXY.transpose()<<std::endl;
  std::cout<<"Corners:"<<std::endl;
  std::cout<<cornerXY.transpose()<<std::endl;
  std::cout<<"BCs: "<<std::endl;
  for(int ii = 0; ii < numCells; ii++)
  {
    std::cout<<"cell "<<ii<<std::endl;
    std::cout<<cellBCTypes[ii]<<std::endl;
    std::cout<<std::endl;
  }
  std::cout<<"Faces: "<<std::endl;
  for(int ii = 0; ii < numCells; ii++)
  {
    std::cout<<"cell "<<ii<<std::endl;
    std::cout<<faceMap.col(ii)<<std::endl;
    std::cout<<std::endl;
  }
  std::cout<<"CornerMaps: "<<std::endl;
  for(int ii = 0; ii < numCells; ii++)
  {
    std::cout<<"cell "<<ii<<std::endl;
    std::cout<<cornerMap.col(ii)<<std::endl;
    std::cout<<std::endl;
  }


}

void CartesianGrid::setBCs()
{
  // Set top and bottom to be no slip and no penetration
  for(int cellNum = 0; cellNum < numCells; cellNum++)
  {
    if(cellNum < Nx)
    {
      // Bottom of cells
      cellBCTypes[cellNum](0,0) = 1;
      cellBCTypes[cellNum](1,0) = 2; // No slip in x
      cellBCTypes[cellNum](1,1) = 1; // No penetration in y
    }
    else if (cellNum >= numCells - Nx)
    {
      // Top of cells
      cellBCTypes[cellNum](0,0) = 1;
      cellBCTypes[cellNum](3,0) = 2; // No slip in x
      cellBCTypes[cellNum](3,1) = 1; // No penetration in y
    }

    if(cellNum % Nx == 0)
    {
      // Left side
      cellBCTypes[cellNum](0,0) = 1;
      cellBCTypes[cellNum](4,2) = 1; // Dirichlet pressure at inlet
    }
    else if((cellNum+1) % Nx == 0)
    {
      // Right side
      cellBCTypes[cellNum](0,0) = 1;
      cellBCTypes[cellNum](2,2) = 2; // Dirichlet pressure at outlet
    }
  }

  uVelBCValue.resize(2);
  uVelBCValue(0) = 1.0;
  uVelBCValue(1) = 0.0;

  vVelBCValue.resize(1);
  vVelBCValue(0) = 0.0;

  pressureBCValue.resize(2);
  pressureBCValue(0) = 1000;
  pressureBCValue(1) = 200;
}

void CartesianGrid::setInitialValues()
{
  uVel.setZero();
  vVel.setZero();
  pressure.setZero();
}

void CartesianGrid::setVariableSizes(int nCells, int nFaces, int nCorners)
{
  numCells = nCells;
  numFaces = nFaces;
  numCorners = nCorners;

  centerXY.setZero(numCells*nDim);
  cornerXY.setZero(numCorners*nDim);
  uVel.setZero(numCells);
  vVel.setZero(numCells);
  pressure.setZero(numCells);
  cornerMap.setZero(numCornersPerElement,numCells);
  faceMap.setZero(numFacesPerElement,numCells);
  cellBCTypes.resize(numCells);
  for(int ii = 0; ii < numCells; ii++)
    cellBCTypes[ii].setZero(1+numFacesPerElement,3);
}

void CartesianGrid::generateMesh()
{
  for(int kk = 0; kk < Ny+1; kk++)
  {
    for(int ii = 0; ii < Nx+1; ii++)
    {
      if(ii != Nx && kk != Ny)
      {
        int centerFlatIdx = (ii + kk*Nx)*nDim;
        centerXY(centerFlatIdx) = delX/2 + ii*delX;
        centerXY(centerFlatIdx+1) = delY/2 + kk*delY;
      }

      int cornerFlatIdx = (ii+kk*(Nx+1))*nDim;
      cornerXY(cornerFlatIdx) = delX*ii;
      cornerXY(cornerFlatIdx+1) = delY*kk;
    }
  }
  std::cout<<"Got here 2"<<std::endl;
  for(int kk = 0; kk < Ny; kk++)
  {
    for(int ii = 0; ii < Nx; ii++)
    {
      int cellFlatIdx = ii + kk*Nx;
      cornerMap(0,cellFlatIdx) = ii + kk*(Nx+1);
      cornerMap(1,cellFlatIdx) = ii+1 + kk*(Nx+1);
      cornerMap(2,cellFlatIdx) = ii+1 + (kk+1)*(Nx+1);
      cornerMap(3,cellFlatIdx) = ii + (kk+1)*(Nx+1);

      int offset = (2*Nx+1)*kk;
      faceMap(0,cellFlatIdx) = ii + offset;
      faceMap(1,cellFlatIdx) = ii+Nx+1 + offset;
      faceMap(2,cellFlatIdx) = ii + offset + 2*Nx+1;
      faceMap(3,cellFlatIdx) = ii+Nx + offset;
    }
  }
}