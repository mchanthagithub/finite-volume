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
  std::cout<<"Centers:"<<std::endl;
  std::cout<<centerXY.transpose()<<std::endl;
  std::cout<<"Corners:"<<std::endl;
  std::cout<<cornerXY.transpose()<<std::endl;
}

void CartesianGrid::setBCs()
{
  // Set top and bottom to be no slip and no penetration
  for(int cellNum = 0; cellNum < numCells; cellNum++)
  {
    if(cellNum < Nx)
    {
      cellBCTypes[cellNum](0,0) = 1;
      cellBCTypes[cellNum](1,0) = 1;

    }


  }
}

void CartesianGrid::setInitialValues()
{
  ;
}

void CartesianGrid::setVariableSizes(int nCells, int nFaces, int nCorners)
{
  numCells = nCells;
  numFaces = nFaces;
  numCorners = nCorners;

  centerXY.resize(numCells*nDim);
  cornerXY.resize(numCorners*nDim);
  uVel.resize(numCells);
  vVel.resize(numCells);
  pressure.resize(numCells);
  cornerMap.resize(nCorners,numCornersPerElement);
  faceMap.resize(numFaces,numFacesPerElement);
  cellBCTypes.resize(numCells);
  for(int ii = 0; ii < numCells; ii++)
    cellBCTypes[ii].resize(1+numFacesPerElement,2);
}

void CartesianGrid::generateMesh()
{
  for(int kk = 0; kk < Ny+1; kk++)
  {
    for(int ii = 0; ii < Nx+1; ii++)
    {
      if(ii != Nx && kk != Ny)
      {
        int centerFlatID = (ii + kk*Nx)*nDim;
        centerXY(centerFlatID) = delX/2 + ii*delX;
        centerXY(centerFlatID+1) = delY/2 + kk*delY;
      }

      int cornerFlatID = (ii+kk*(Nx+1))*nDim;
      cornerXY(cornerFlatID) = delX*ii;
      cornerXY(cornerFlatID+1) = delY*kk;
    }
  }

}