//
// Created by maytee on 4/26/18.
//

#include "CartesianGrid.h"
#include<iostream>
// Nx = number of cells in X, Ny = number of cells in Y, Lx = length in X, Ly = length in Y
CartesianGrid::CartesianGrid()
{
}
CartesianGrid::CartesianGrid(int in_Nx, int in_Ny, double in_Lx, double in_Ly, double smallX, double smallY)
{
  Nx = in_Nx;
  Ny = in_Ny;
  Lx = in_Lx;
  Ly = in_Ly;
  delX = Lx/Nx;
  delY = Ly/Ny;
  minX = smallX;
  minY = smallY;
  std::cout<<Nx<<" "<<Ny<<" "<<Lx<<" "<<Ly<<" "<<delX<<" "<<delY<<std::endl;

  int nFaces = Nx*(Ny+1) + Ny*(Nx+1);
  int nElem = Nx*Ny;
  int nCorners = (Nx+1)*(Ny+1);

  setVariableSizes(nElem,nFaces,nCorners);
  generateMesh();
  setBCs();

  std::cout<<"Centers:"<<std::endl;
  //std::cout<<centerXY.transpose()<<std::endl;
  std::cout<<"Corners:"<<std::endl;
  //std::cout<<cornerXY.transpose()<<std::endl;
  std::cout<<"BCs: "<<std::endl;
  for(int ii = 0; ii < numCells; ii++)
  {
    //std::cout<<"cell "<<ii<<std::endl;
    //std::cout<<BCTypes<<std::endl;
    //std::cout<<std::endl;
  }
  std::cout<<"Faces: "<<std::endl;
  for(int ii = 0; ii < numCells; ii++)
  {
    //std::cout<<"cell "<<ii<<std::endl;
    //std::cout<<faceMap.col(ii)<<std::endl;
    //std::cout<<std::endl;
  }
  std::cout<<"CornerMaps: "<<std::endl;
  for(int ii = 0; ii < numCells; ii++)
  {
    //std::cout<<"cell "<<ii<<std::endl;
    //std::cout<<cornerMap.col(ii)<<std::endl;
    //std::cout<<std::endl;
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
      cellHasVelocityBC(cellNum) = 1;
      BCTypes(faceMap(0,cellNum),0) = 2; // No slip in x
      BCTypes(faceMap(0,cellNum),1) = 1; // No penetration in y
    }
    else if (cellNum >= numCells - Nx)
    {
      // Top of cells
      cellHasVelocityBC(cellNum) = 1;
      BCTypes(faceMap(2,cellNum),0) = 2; // No slip in x
      BCTypes(faceMap(2,cellNum),1) = 1; // No penetration in y
    }

    if(cellNum % Nx == 0)
    {
      // Left side
      cellHasVelocityBC(cellNum) = 1;
      BCTypes(faceMap(3,cellNum),0) = 2; // No penetration in x
      BCTypes(faceMap(3,cellNum),1) = 1; // No slip in y

      //cellHasPressureBC(cellNum) = 1;
      //BCTypes(faceMap(3,cellNum),2) = 1; // Dirichlet presure at inlet

      //cellHasVelocityBC(cellNum) = 1;
      //BCTypes(faceMap(3,cellNum),0) = 1; // Entrance vel, note turn on either pressure BC or vel BC
    }
    else if((cellNum+1) % Nx == 0)
    {
      // Right side
      //cellHasPressureBC(cellNum) = 1;
      //BCTypes(faceMap(1,cellNum),2) = 1; // Dirichlet presure at otulet

      //cellHasVelocityBC(cellNum) = 1;
      //BCTypes(faceMap(1,cellNum),0) = 1; // Exit vel, note turn on either pressure BC or vel BC
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

  for(int kk = 0; kk < Ny; kk++)
  {
    for(int ii = 0; ii < Nx; ii++)
    {
      int flatIdx = ii + kk*Nx;
      double xCoord = centerXY(flatIdx*2);
      double yCoord = centerXY(flatIdx*2+1);
      double pi = 3.14159;
      if(xCoord <= 1.0 && yCoord <= 1.0)
      {
        //uVel(flatIdx) = 0.25*std::sin(pi*xCoord);
        //vVel(flatIdx) = 0.25*std::sin(pi*yCoord);
      }
      else if((xCoord <= 2.0 && xCoord >=1.0) && (yCoord <= 2.0 && yCoord >=1.0))
      {
        uVel(flatIdx) = 0.25*std::sin(pi*(xCoord-1.0));
        vVel(flatIdx) = 0.25*std::sin(pi*(yCoord-1.0));
      }
    }
  }
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

  cellHasPressureBC.setZero(numCells);
  cellHasVelocityBC.setZero(numCells);
  BCTypes.setZero(numFaces,nDim+1);
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
        centerXY(centerFlatIdx) = delX/2 + ii*delX + minX;
        centerXY(centerFlatIdx+1) = delY/2 + kk*delY + minY;
      }

      int cornerFlatIdx = (ii+kk*(Nx+1))*nDim;
      cornerXY(cornerFlatIdx) = delX*ii;
      cornerXY(cornerFlatIdx+1) = delY*kk;
    }
  }
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