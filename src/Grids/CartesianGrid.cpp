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
  //setChannel();
  //setChannelWithObstacle();
  //setChannelWithCylinder();
  //setBurgers();
  setDiffusion();


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

void CartesianGrid::setChannel()
{
  //uVel.setZero();
  uVel.setOnes();
  uVel = 0.25*uVel;
  vVel.setZero();
  pressure = 1000.0*Eigen::VectorXd::Ones(numCells);
  isChannel = true;
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

        faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

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
  isAllNeumannPressureBCs = false;

  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 2.0;
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

void CartesianGrid::setChannelWithObstacle()
{
  uVel.setZero();
  uVel = 0.25*uVel;
  vVel.setZero();
  pressure = 1000.0*Eigen::VectorXd::Ones(numCells);
  isChannel = true;
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

        faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

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
  isAllNeumannPressureBCs = false;

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

  int xBegin = (2*Nx)/8;
  int xEnd = (4*Nx)/8;

  int yBegin = (2*Ny)/6;
  int yEnd = (4*Ny)/6;

  bool useSquare = true;
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

          cellHasDirichletVelocityBC(southElemIdx) = 1;
          BCTypes(faceMap(2, southElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(2, southElemIdx), 1) = 1; // No slip in y
        }
        if(kk == yEnd)
        {
          cellHasPressureBC(northElemIdx,0) = -1;
          cellHasPressureBC(northElemIdx,1) = 1;

          cellHasDirichletVelocityBC(northElemIdx) = 1;
          BCTypes(faceMap(0, northElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(0, northElemIdx), 1) = 1; // No slip in y
        }
        if(ii == xBegin)
        {
          cellHasPressureBC(westElemIdx,0) = -1;
          cellHasPressureBC(westElemIdx,1) = 2;

          cellHasDirichletVelocityBC(westElemIdx) = 1;
          BCTypes(faceMap(1, westElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(1, westElemIdx), 1) = 1; // No slip in y
        }
        if(ii == xEnd)
        {
          cellHasPressureBC(eastElemIdx,0) = -1;
          cellHasPressureBC(eastElemIdx,1) = 4;

          cellHasDirichletVelocityBC(eastElemIdx) = 1;
          BCTypes(faceMap(3, eastElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(3, eastElemIdx), 1) = 1; // No slip in y
        }
      }
    }

    for(int kk = 0; kk < Ny; kk++)
    {
      for (int ii = 0; ii < Nx; ii++)
      {
        int flatElemIdx = ii + kk*Nx;

        // Square
        if(ii >= xBegin && ii <= xEnd && kk >= yBegin && kk <= yEnd)
        {
          cellHasDirichletVelocityBC(flatElemIdx) = 1;
          cellHasPressureBC(flatElemIdx,0) = -1;
          for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
          {
            BCTypes(faceMap(faceNum, flatElemIdx), 0) = 1; // No slip in x
            BCTypes(faceMap(faceNum, flatElemIdx), 1) = 1; // No penetration in y
          }
        }
      }
    }
  }


  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 1.0;
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
  pressureDirichletBCValue(0) = 0.0;
  pressureDirichletBCValue(1) = 0.0;

  pressureNeumannBCValue.resize(1);
  pressureNeumannBCValue(0) = 0.0;
}

void CartesianGrid::setChannelWithCylinder()
{
  uVel.setZero();
  vVel.setZero();
  pressure = 1000.0*Eigen::VectorXd::Ones(numCells);
  isChannel = true;
  for(int kk = 0; kk < Ny; kk++)
  {
    for(int ii = 0; ii < Nx; ii++)
    {
      int globalIdx = ii + kk*Nx;
      if(kk < Ny/2)
        uVel(globalIdx) = 0.85;
      else
        uVel(globalIdx) = 0.80;
    }
  }
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

        faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

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
  isAllNeumannPressureBCs = false;

  double x0 = 0.80;
  double y0 = 0.5;
  double radius = 0.15;

  for(int ii = 0; ii < Nx; ii++)
  {
    for(int kk = 0; kk < Ny; kk++)
    {
      int globalIdx = ii+kk*Nx;
      double xCoord = centerXY(globalIdx*2);
      double yCoord = centerXY(globalIdx*2+1);
      double distance = sqrt((xCoord-x0)*(xCoord-x0) + (yCoord-y0)*(yCoord-y0) );
      if(distance < radius)
      {
        // Mask out cells
        uVel(globalIdx) = 0.0;
        vVel(globalIdx) = 0.0;

        cellHasDirichletVelocityBC(globalIdx) = 1;
        cellHasPressureBC(globalIdx,0) = -1;
        for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
        {
          BCTypes(faceMap(faceNum, globalIdx), 0) = 1; // No slip in x
          BCTypes(faceMap(faceNum, globalIdx), 1) = 1; // No penetration in y
        }
      }
    }
  }


  for(int ii = 0; ii < Nx; ii++)
  {
    for(int kk = 0; kk < Ny; kk++)
    {
      int globalIdx = ii+kk*Nx;
      int westIdx = globalIdx - 1;
      int eastIdx = globalIdx + 1;
      int northIdx = globalIdx + Nx;
      int southIdx = globalIdx - Nx;

      bool currentIsMasked = true;
      bool westIsMasked = true;
      bool eastIsMasked = true;
      bool northIsMasked = true;
      bool southIsMasked = true;
      for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
      {
        if( BCTypes(faceMap(faceNum, globalIdx), 0) == 0)
          currentIsMasked = false;

        if(ii != 0)
          if( BCTypes(faceMap(faceNum, westIdx), 0) == 0)
            westIsMasked = false;

        if(ii != Nx-1)
          if( BCTypes(faceMap(faceNum, eastIdx), 0) == 0)
            eastIsMasked = false;

        if(kk != 0)
          if( BCTypes(faceMap(faceNum, southIdx), 0) == 0)
            southIsMasked = false;

        if(kk != Ny-1)
          if( BCTypes(faceMap(faceNum, northIdx), 0) == 0)
            northIsMasked = false;
      }

      if(currentIsMasked && westIsMasked && eastIsMasked && northIsMasked && southIsMasked)
        continue;

      if(currentIsMasked)
      {
        if(!westIsMasked)
        {
          if(cellHasPressureBC(westIdx,1) == 0)
          {
            cellHasPressureBC(westIdx,0) = -1;
            cellHasPressureBC(westIdx,1) = 2;
          }
          else if (cellHasPressureBC(westIdx,1) == 1)
          {
            cellHasPressureBC(westIdx,0) = -1;
            cellHasPressureBC(westIdx,1) = 5;
          }
          else if (cellHasPressureBC(westIdx,1) == 3)
          {
            cellHasPressureBC(westIdx,0) = -1;
            cellHasPressureBC(westIdx,1) = 6;
          }
          cellHasDirichletVelocityBC(westIdx) = 1;
          BCTypes(faceMap(1, westIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(1, westIdx), 1) = 1; // No slip in y
        }
        if(!northIsMasked)
        {
          if(cellHasPressureBC(northIdx,1) == 0)
          {
            cellHasPressureBC(northIdx,0) = -1;
            cellHasPressureBC(northIdx,1) = 1;
          }
          else if (cellHasPressureBC(northIdx,1) == 2)
          {
            cellHasPressureBC(northIdx,0) = -1;
            cellHasPressureBC(northIdx,1) = 5;
          }
          else if (cellHasPressureBC(northIdx,1) == 4)
          {
            cellHasPressureBC(northIdx,0) = -1;
            cellHasPressureBC(northIdx,1) = 8;
          }
          cellHasDirichletVelocityBC(northIdx) = 1;
          BCTypes(faceMap(0, northIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(0, northIdx), 1) = 1; // No slip in y
        }
        if(!eastIsMasked)
        {
          if(cellHasPressureBC(eastIdx,1) == 0)
          {
            cellHasPressureBC(eastIdx,0) = -1;
            cellHasPressureBC(eastIdx,1) = 4;
          }
          else if (cellHasPressureBC(eastIdx,1) == 1)
          {
            cellHasPressureBC(eastIdx,0) = -1;
            cellHasPressureBC(eastIdx,1) = 8;
          }
          else if (cellHasPressureBC(eastIdx,1) == 3)
          {
            cellHasPressureBC(eastIdx,0) = -1;
            cellHasPressureBC(eastIdx,1) = 7;
          }
          cellHasDirichletVelocityBC(eastIdx) = 1;
          BCTypes(faceMap(3, eastIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(3, eastIdx), 1) = 1; // No slip in y
        }

        if(!southIsMasked)
        {
          if(cellHasPressureBC(southIdx,1) == 0)
          {
            cellHasPressureBC(southIdx,0) = -1;
            cellHasPressureBC(southIdx,1) = 3;
          }
          else if (cellHasPressureBC(southIdx,1) == 2)
          {
            cellHasPressureBC(southIdx,0) = -1;
            cellHasPressureBC(southIdx,1) = 6;
          }
          else if (cellHasPressureBC(southIdx,1) == 4)
          {
            cellHasPressureBC(southIdx,0) = -1;
            cellHasPressureBC(southIdx,1) = 7;
          }
          cellHasDirichletVelocityBC(southIdx) = 1;
          BCTypes(faceMap(2, southIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(2, southIdx), 1) = 1; // No slip in y
        }
      }
    }
  }

  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 1.0;
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
  pressureDirichletBCValue(0) = 0.0;
  pressureDirichletBCValue(1) = 0.0;

  pressureNeumannBCValue.resize(1);
  pressureNeumannBCValue(0) = 0.0;
}


void CartesianGrid::setBurgers()
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
        uVel(flatIdx) = 0.25*std::sin(pi*xCoord);
        vVel(flatIdx) = 0.25*std::sin(pi*yCoord);
      }
      else if((xCoord <= 2.0 && xCoord >=1.0) && (yCoord <= 2.0 && yCoord >=1.0))
      {
        //uVel(flatIdx) = 0.5*std::sin(pi*(xCoord-1.0));
        //vVel(flatIdx) = 0.5*std::sin(pi*(yCoord-1.0));
      }
    }
  }

  isBurgers = true;
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
        BCTypes(faceMap(3, cellNum), 0) = 1; // No penetration in x
        BCTypes(faceMap(3, cellNum), 1) = 1; // No slip in y

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

        faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

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
  isAllNeumannPressureBCs = false;

  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 1.0;
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

void CartesianGrid::setDiffusion()
{
  uVel.setZero();
  vVel.setZero();
  pressure = 1000.0*Eigen::VectorXd::Ones(numCells);

  isDiffusion = true;
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
        BCTypes(faceMap(2, cellNum), 0) = 2; // No slip in x
        BCTypes(faceMap(2, cellNum), 1) = 1; // No penetration in y

        cellHasPressureBC(cellNum, 0) = -1;
        cellHasPressureBC(cellNum, 1) = 3;
      }

      if(ii == 0)
      {
        // Left side
        cellHasDirichletVelocityBC(cellNum) = 1;
        BCTypes(faceMap(3, cellNum), 0) = 1; // No penetration in x
        BCTypes(faceMap(3, cellNum), 1) = 1; // No slip in y

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

        cellHasDirichletVelocityBC(cellNum) = 1;
        BCTypes(faceMap(1, cellNum), 0) = 1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = 1; // No slip in y

        //faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        //BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        //BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

        //int temp = faceMap(1,cellNum);
        //int sizetep = faceVelocityNeumannBCs[0].size();
        //faceVelocityNeumannBCs[faceMap(1,cellNum)](0,0) = 1;
        //faceVelocityNeumannBCs[faceMap(1,cellNum)](1,0) = 1;

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
  isAllNeumannPressureBCs = false;

  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 5.0;
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

        faceHasNeumannVelocityBC(faceMap(1,cellNum)) = 1;
        BCTypes(faceMap(1, cellNum), 0) = -1; // No penetration in x
        BCTypes(faceMap(1, cellNum), 1) = -1; // No slip in y

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
  isAllNeumannPressureBCs = false;

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

  int xBegin = (1*Nx)/8;
  int xEnd = (3*Nx)/8;

  int yBegin = (2*Ny)/6;
  int yEnd = (4*Ny)/6;

  bool useSquare = true;
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

          cellHasDirichletVelocityBC(southElemIdx) = 1;
          BCTypes(faceMap(2, southElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(2, southElemIdx), 1) = 1; // No slip in y
        }
        if(kk == yEnd)
        {
          cellHasPressureBC(northElemIdx,0) = -1;
          cellHasPressureBC(northElemIdx,1) = 1;

          cellHasDirichletVelocityBC(northElemIdx) = 1;
          BCTypes(faceMap(0, northElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(0, northElemIdx), 1) = 1; // No slip in y
        }
        if(ii == xBegin)
        {
          cellHasPressureBC(westElemIdx,0) = -1;
          cellHasPressureBC(westElemIdx,1) = 2;

          cellHasDirichletVelocityBC(westElemIdx) = 1;
          BCTypes(faceMap(1, westElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(1, westElemIdx), 1) = 1; // No slip in y
        }
        if(ii == xEnd)
        {
          cellHasPressureBC(eastElemIdx,0) = -1;
          cellHasPressureBC(eastElemIdx,1) = 4;

          cellHasDirichletVelocityBC(eastElemIdx) = 1;
          BCTypes(faceMap(3, eastElemIdx), 0) = 1; // No penetration in x
          BCTypes(faceMap(3, eastElemIdx), 1) = 1; // No slip in y
        }
      }
    }

    for(int kk = 0; kk < Ny; kk++)
    {
      for (int ii = 0; ii < Nx; ii++)
      {
        int flatElemIdx = ii + kk*Nx;

        // Square
        if(ii >= xBegin && ii <= xEnd && kk >= yBegin && kk <= yEnd)
        {
          cellHasDirichletVelocityBC(flatElemIdx) = 1;
          cellHasPressureBC(flatElemIdx,0) = -1;
          for(int faceNum = 0; faceNum < numFacesPerElement; faceNum++)
          {
            BCTypes(faceMap(faceNum, flatElemIdx), 0) = 1; // No slip in x
            BCTypes(faceMap(faceNum, flatElemIdx), 1) = 1; // No penetration in y
          }
        }
      }
    }
  }


  uVelBCValue.resize(3);
  uVelBCValue(0) = 0.0;
  uVelBCValue(1) = 1.0;
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

  oldUVel.setZero(numCells);
  oldVVel.setZero(numCells);

  pressure.setZero(numCells);
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
  faceVels.setZero(numFaces*nDim);
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
      //if(cellHasPressureBC(cellNum,0) > 0)
      //  mappingGlobalToActive(cellNum) = -1;
      // Is a Neumann BC
      //else
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

  //pressure.setZero(numCells);
}
