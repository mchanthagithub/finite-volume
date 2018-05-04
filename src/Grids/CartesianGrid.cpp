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
  createMappings();
}

void CartesianGrid::setBCs() {
  // Set top and bottom to be no slip and no penetration

  for(int kk = 0; kk < Ny; kk++)
  {
    for(int ii = 0; ii < Nx; ii++)
    {
      int cellNum = ii + kk*Nx;
      if(kk == 0)
      {
        // Bottom of cells
        cellHasDirichletVelocityBC(cellNum) = 1;
        BCTypes(faceMap(0, cellNum), 0) = 1; // No slip in x
        BCTypes(faceMap(0, cellNum), 1) = 1; // No penetration in y

        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 1;
      }
      else if(kk == Ny-1)
      {
        // Top of cells
        cellHasDirichletVelocityBC(cellNum) = 1;
        BCTypes(faceMap(2, cellNum), 0) = 1; // No slip in x
        BCTypes(faceMap(2, cellNum), 1) = 1; // No penetration in y

        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 3;
      }

      if(ii == 0)
      {
        // Left side
        cellHasDirichletVelocityBC(cellNum) = 1;
        BCTypes(faceMap(3, cellNum), 0) = 2; // No penetration in x
        BCTypes(faceMap(3, cellNum), 1) = 1; // No slip in y

        // Dirichlet pressure BC at inlet
        //if(cellNum >= Nx && cellNum < numCells - Nx)
        {
          //cellHasPressureBC(cellNum, 0) = 1;
          //cellHasPressureBC(cellNum, 1) = -1;
        }

        if(kk != Ny-1)
        {
          cellHasPressureBC(cellNum, 0) = -1;
          cellHasPressureBC(cellNum, 1) = 4;
        }

      }
      else if(ii == Nx-1)
      {
        // Right side
        //cellHasPressureBC(cellNum) = 1;
        //BCTypes(faceMap(1,cellNum),2) = 1; // Dirichlet presure at otulet

        cellHasDirichletVelocityBC(cellNum) = -1 ; // Neumann
        BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

        faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        int temp = faceMap(1,cellNum);
        int sizetep = faceVelocityNeumannBCs[0].size();
        faceVelocityNeumannBCs[faceMap(1,cellNum)](0,0) = 1;
        faceVelocityNeumannBCs[faceMap(1,cellNum)](1,0) = 1;

        if(kk != Ny-1)
        {
          cellHasPressureBC(cellNum, 0) = -1;
          cellHasPressureBC(cellNum, 1) = 2;
        }

      }

      if (ii == Nx - 1 && kk == 0) {
        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 5; // Bottom right
      }
      else if (ii == Nx - 1 && kk == Ny - 1)
      {
        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 6; // Top right
      }
      else if (ii == 0 && kk == Ny - 1)
      {
        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 7; // Top left
      }

      else if (ii == 0 && kk == 0)
      {
        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 8; // Bottom Left
      }
    }
  }
  isAllNeumannPressureBCs = true;

  /*
  double height = 1.5;
  double length = 1.5;
  double minX = 1.25;
  double minY = 1.25;

  for(int kk = 0; kk < Ny; kk++)
  {
    for(int ii = 0; ii < Nx; ii++)
    {
      int flatElemIdx = ii + kk*Nx;
      double xCoord = centerXY(flatElemIdx*2);
      double yCoord = centerXY(flatElemIdx*2+1);
      if(xCoord >= minX && xCoord <= minX+length && yCoord >= minY && yCoord <= minY+height)
      {
        cellHasDirichletVelocityBC(flatElemIdx,0) = 1;
        cellHasPressureBC(flatElemIdx,0) = -1;
        for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
        {
          BCTypes(faceMap(faceNum,flatElemIdx),0) = 2;
          BCTypes(faceMap(faceNum,flatElemIdx),1) = 1;
        }
      }
    }
  }
  */

  int xBegin = (2*Nx)/6;
  int xEnd = (4*Nx)/6;

  int yBegin = (2*Ny)/6;
  int yEnd = (4*Ny)/6;

  bool useSquare = false;
  //bool useSquare = true;
  if (useSquare)
  {
    for(int kk = yBegin; kk <= yEnd; kk++)
    {
      for(int ii = xBegin; ii <= xEnd; ii++)
      {
        int flatElemIdx = ii + kk*Nx;
        int eastElemIdx =flatElemIdx + 1;
        int westElemIdx =flatElemIdx - 1;
        int southElemIdx =flatElemIdx - Nx;
        int northElemIdx = flatElemIdx + Nx;
        if(kk == yBegin)
        {
          cellHasPressureBC(southElemIdx,0) = -1;
          cellHasPressureBC(southElemIdx,1) = 3;
        }
        if(kk == yEnd)
        {
          cellHasPressureBC(northElemIdx,0) = -1;
          cellHasPressureBC(northElemIdx,1) = 1;
        }
        if(ii == xBegin)
        {
          cellHasPressureBC(westElemIdx,0) = -1;
          cellHasPressureBC(westElemIdx,1) = 2;
        }
        if(ii == xEnd)
        {
          cellHasPressureBC(eastElemIdx,0) = -1;
          cellHasPressureBC(eastElemIdx,1) = 4;
        }
      }
    }

    for(int kk = 0; kk < Ny; kk++)
    {
      for (int ii = 0; ii < Nx; ii++)
      {
        int flatElemIdx = ii + kk*Nx;
        int eastElemIdx =flatElemIdx + 1;
        int westElemIdx =flatElemIdx + 1;
        int southElemIdx =flatElemIdx - Nx;
        int northElemIdx = flatElemIdx + Nx;

        // Square
        if(ii >= xBegin && ii <= xEnd && kk >= yBegin && kk <= yEnd)
        {
          cellHasDirichletVelocityBC(flatElemIdx) = 1;
          cellHasPressureBC(flatElemIdx,0) = -1;
          for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
          {
            BCTypes(faceMap(faceNum,flatElemIdx),0) = 2;
            BCTypes(faceMap(faceNum,flatElemIdx),1) = 1;
          }
        }
      }
    }
  }


  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 10.0;
  uVelBCValue(2) = -1.0;

  uVelNeumannBCValue.resize(1);
  uVelNeumannBCValue(0) = 0.0;

  vVelNeumannBCValue.resize(1);
  vVelNeumannBCValue(0) = 0.0;

  vVelBCValue.resize(3);
  vVelBCValue(0) = 0.0;
  vVelBCValue(1) = 1.0;
  vVelBCValue(2) = -1.0;

  pressureDirichletBCValue.resize(2);
  pressureDirichletBCValue(0) = 1000;
  pressureDirichletBCValue(1) = 200;

  pressureNeumannBCValue.resize(1);
  pressureNeumannBCValue(0) = 0.0;
}

void CartesianGrid::setInitialValues()
{
  uVel.setZero();
  vVel.setZero();
  pressure = 1000.0*Eigen::VectorXd::Ones(numCells);


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
        //uVel(flatIdx) = 0.5*std::sin(pi*(xCoord-1.0));
        //vVel(flatIdx) = 0.5*std::sin(pi*(yCoord-1.0));
      }
    }
  }
}

void CartesianGrid::setVariableSizes(int nCells, int nFaces, int nCorners)
{
  numCells = nCells;
  numFaces = nFaces;
  numCorners = nCorners;
  nDim = 2;
  centerXY.setZero(numCells*nDim);
  cornerXY.setZero(numCorners*nDim);
  uVel.setZero(numCells);
  vVel.setZero(numCells);
  pressure.setZero(numCells);
  pressureGradient.setZero(numCells,nDim);
  cornerMap.setZero(numCornersPerElement,numCells);
  faceMap.setZero(numFacesPerElement,numCells);

  cellHasPressureBC.setZero(numCells,nDim);
  cellHasDirichletVelocityBC.setZero(numCells);
  BCTypes.setZero(numFaces,nDim+1);
  advectionFluxes.setZero(numCells,nDim);
  diffusionFluxes.setZero(numCells,nDim);
  faceVelocityNeumannBCs.resize(numFaces);
  Eigen::MatrixXi temp = Eigen::MatrixXi::Zero(nDim,nDim);
  for(int ii = 0; ii < numFaces; ii++)
    faceVelocityNeumannBCs[ii] = temp;

  faceHasNeumannVelocityBC.resize(numFaces);
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
      cornerXY(cornerFlatIdx) = delX*ii + minX;
      cornerXY(cornerFlatIdx+1) = delY*kk + minY;
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

void CartesianGrid::createMappings()
{
 // Construct mapping between active center nodes and global center nodes
  mappingGlobalToActive.setZero(numCells);
  std::vector<int> mappingActiveToGlobalVector;
  int dofCtr = 0;
  for(int cellNum = 0; cellNum < numCells; cellNum++)
  {
    bool isMasked = false;
    int faceBCCtr = 0;
    for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
    {
      if(BCTypes(faceMap(faceNum,cellNum),0) != 0)
        faceBCCtr++;
    }
    if (faceBCCtr == numFacesPerElement)
      isMasked = true;

    if(isMasked)
      mappingGlobalToActive(cellNum) = -2;
    else
    {
      // Is a dirichlet BC
      if(cellHasPressureBC(cellNum,0) > 0)
        mappingGlobalToActive(cellNum) = -1;
      // Is a Neumann BC
      else
      {
        mappingGlobalToActive(cellNum) = dofCtr;
        mappingActiveToGlobalVector.push_back(cellNum);
        dofCtr++;
      }
    }
  }
  totalDOF = dofCtr;
  mappingActiveToGlobal.setZero(dofCtr);
  for(int ii = 0; ii < dofCtr; ii++)
    mappingActiveToGlobal(ii) = mappingActiveToGlobalVector[ii];

}
