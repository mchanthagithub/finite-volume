//
// Created by maytee on 4/27/18.
//

#include "AdvectionUpwind.h"
#include <iostream>
#include "cstdlib"
void AdvectionUpwind::calculateAdvectionFluxesCartesian(CartesianGrid &grid, InterpolateUpwind& interp)
{
  uVelInterp = interp.uVelInterp;
  vVelInterp = interp.vVelInterp;
  // Calculate advection flux at all surfaces for all cells, then average to obtain final flux term
  // Here use midpoint rule (fluxMid*surfaceArea = totalFlux + O(delX^2))

  //calculateInterpolationValues(grid);
  //std::cout<<"uVelInterp: "<<std::endl;
  //std::cout<<uVelInterp.transpose()<<std::endl;

  //std::cout<<"vVelInterp: "<<std::endl;
  //std::cout<<vVelInterp.transpose()<<std::endl;

  advectionFluxValues.resize(grid.numCells);
  for(int ii = 0; ii < grid.numCells; ii++)
    advectionFluxValues[ii].setZero(grid.numFacesPerElement,grid.nDim);
  totalCellAdvectionFlux.setZero(grid.numCells, grid.nDim);

  Eigen::MatrixXd averagedFaceFluxes;
  averagedFaceFluxes.setZero(grid.numFaces,grid.nDim);
  Eigen::MatrixXd faceNormals;
  faceNormals.setZero(grid.nDim,grid.numFacesPerElement);
  faceNormals(0,0) = 0.0; faceNormals(1,0) = -1.0;
  faceNormals(0,1) = 1.0; faceNormals(1,1) = 0.0;
  faceNormals(0,2) = 0.0; faceNormals(1,2) = 1,0;
  faceNormals(0,3) = -1.0; faceNormals(1,3) = 0.0;

  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    // Masked out cell
    if(grid.mappingGlobalToActive(cellNum) == -2)
      continue;
    /*
    // Are general formulas -> should move to a function in grid to calculate distances between nodes
    double southLength = grid.cornerXY(grid.cornerMap(0,cellNum)*grid.nDim) - grid.cornerXY(grid.cornerMap(1,cellNum)*grid.nDim);
    double northLength = grid.cornerXY(grid.cornerMap(3,cellNum)*grid.nDim) - grid.cornerXY(grid.cornerMap(2,cellNum)*grid.nDim);
    double eastLength = grid.cornerXY(grid.cornerMap(2,cellNum)*grid.nDim+1) - grid.cornerXY(grid.cornerMap(1,cellNum)*grid.nDim+1);
    double westLength = grid.cornerXY(grid.cornerMap(3,cellNum)*grid.nDim+1) - grid.cornerXY(grid.cornerMap(0,cellNum)*grid.nDim+1);
    */

    double southLength = grid.delX;
    double northLength = grid.delX;
    double eastLength = grid.delY;
    double westLength = grid.delY;

    Eigen::VectorXd sideLengths;
    sideLengths.setZero(grid.numFacesPerElement);
    sideLengths(0) = southLength; sideLengths(1) = eastLength; sideLengths(2) = northLength; sideLengths(3) = westLength;
    double cellArea = southLength*eastLength;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      Eigen::VectorXd faceVelocity(grid.nDim);
      faceVelocity << uVelInterp(cellNum,faceNum), vVelInterp(cellNum,faceNum);

      advectionFluxValues[cellNum].row(faceNum) += (faceVelocity *
        (faceVelocity.dot(faceNormals.col(faceNum)))*sideLengths(faceNum)).transpose()*(1.0/cellArea);


    }
  }

  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + kk*grid.Nx;
      for(int faceNum = 0 ; faceNum < grid.numFacesPerElement; faceNum++)
      {
        averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),0) += std::abs(advectionFluxValues[flatElemIdx](faceNum,0)/2.0);
        averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),1) += std::abs(advectionFluxValues[flatElemIdx](faceNum,1)/2.0);
      }
    }
  }
  /*
  std::cout<<"faceFluxes: "<<std::endl;
  for(int ii = 0; ii < grid.numCells; ii++)
  {
    std::cout<<"cellnumflux: "<<ii<<std::endl;
    std::cout<<averagedFaceFluxes(grid.faceMap(0,ii),0)<<" "<<averagedFaceFluxes(grid.faceMap(1,ii),0)<<" "<<averagedFaceFluxes(grid.faceMap(2,ii),0)<<" "<<averagedFaceFluxes(grid.faceMap(3,ii),0)<<std::endl;
    std::cout<<averagedFaceFluxes(grid.faceMap(0,ii),1)<<" "<<averagedFaceFluxes(grid.faceMap(1,ii),1)<<" "<<averagedFaceFluxes(grid.faceMap(2,ii),1)<<" "<<averagedFaceFluxes(grid.faceMap(3,ii),1)<<std::endl;
  }
  */
  double signX, signY;
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    // Masked out cell
    if(grid.mappingGlobalToActive(cellNum) == -2)
      continue;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      //signX = (advectionFluxValues[cellNum](faceNum,0) > 0) - (advectionFluxValues[cellNum](faceNum,0) < 0);
      //signY = (advectionFluxValues[cellNum](faceNum,1) > 0) - (advectionFluxValues[cellNum](faceNum,1) < 0);
      //totalCellAdvectionFlux(cellNum,0) += averagedFaceFluxes(grid.faceMap(faceNum,cellNum),0)*signX;
      //totalCellAdvectionFlux(cellNum,1) += averagedFaceFluxes(grid.faceMap(faceNum,cellNum),1)*signY;
      totalCellAdvectionFlux(cellNum,0) += advectionFluxValues[cellNum](faceNum,0);
      totalCellAdvectionFlux(cellNum,1) += advectionFluxValues[cellNum](faceNum,1);
    }
  }
  //std::cout<<"Total cell fluxes: "<<std::endl;
  //std::cout<<totalCellAdvectionFlux<<std::endl;
}

void AdvectionUpwind::calculateAdvectionFluxesCartesian(CartesianGrid &grid, InterpolateQUICK& interp)
{
  uVelInterp = interp.uVelInterp;
  vVelInterp = interp.vVelInterp;
  // Calculate advection flux at all surfaces for all cells, then average to obtain final flux term
  // Here use midpoint rule (fluxMid*surfaceArea = totalFlux + O(delX^2))

  //calculateInterpolationValues(grid);
  //std::cout<<"uVelInterp: "<<std::endl;
  //std::cout<<uVelInterp.transpose()<<std::endl;

  //std::cout<<"vVelInterp: "<<std::endl;
  //std::cout<<vVelInterp.transpose()<<std::endl;

  advectionFluxValues.resize(grid.numCells);
  for(int ii = 0; ii < grid.numCells; ii++)
    advectionFluxValues[ii].setZero(grid.numFacesPerElement,grid.nDim);
  totalCellAdvectionFlux.setZero(grid.numCells, grid.nDim);

  Eigen::MatrixXd averagedFaceFluxes;
  averagedFaceFluxes.setZero(grid.numFaces,grid.nDim);

  Eigen::MatrixXd faceNormals;
  faceNormals.setZero(grid.nDim,grid.numFacesPerElement);
  faceNormals(0,0) = 0.0; faceNormals(1,0) = -1.0;
  faceNormals(0,1) = 1.0; faceNormals(1,1) = 0.0;
  faceNormals(0,2) = 0.0; faceNormals(1,2) = 1,0;
  faceNormals(0,3) = -1.0; faceNormals(1,3) = 0.0;

  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    /*
    // Are general formulas -> should move to a function in grid to calculate distances between nodes
    double southLength = grid.cornerXY(grid.cornerMap(0,cellNum)*grid.nDim) - grid.cornerXY(grid.cornerMap(1,cellNum)*grid.nDim);
    double northLength = grid.cornerXY(grid.cornerMap(3,cellNum)*grid.nDim) - grid.cornerXY(grid.cornerMap(2,cellNum)*grid.nDim);
    double eastLength = grid.cornerXY(grid.cornerMap(2,cellNum)*grid.nDim+1) - grid.cornerXY(grid.cornerMap(1,cellNum)*grid.nDim+1);
    double westLength = grid.cornerXY(grid.cornerMap(3,cellNum)*grid.nDim+1) - grid.cornerXY(grid.cornerMap(0,cellNum)*grid.nDim+1);
    */

    double southLength = grid.delX;
    double northLength = grid.delX;
    double eastLength = grid.delY;
    double westLength = grid.delY;

    Eigen::VectorXd sideLengths;
    sideLengths.setZero(grid.numFacesPerElement);
    sideLengths(0) = southLength; sideLengths(1) = eastLength; sideLengths(2) = northLength; sideLengths(3) = westLength;
    double cellArea = southLength*eastLength;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      Eigen::VectorXd faceVelocity(grid.nDim);
      faceVelocity << uVelInterp(cellNum,faceNum), vVelInterp(cellNum,faceNum);

      advectionFluxValues[cellNum].row(faceNum) += (faceVelocity *
        (faceVelocity.dot(faceNormals.col(faceNum)))*sideLengths(faceNum)).transpose()*(1.0/cellArea);

      Eigen::VectorXd temp = (faceVelocity *
        (faceVelocity.dot(faceNormals.col(faceNum)))*sideLengths(faceNum)).transpose();

    }
  }

  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + kk*grid.Nx;
      int northElemIdx = flatElemIdx + grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;
      int eastElemIdx = flatElemIdx + 1;
      int westElemIdx = flatElemIdx - 1;

      int nnElemIdx = flatElemIdx + 2*grid.Nx;
      int ssElemIdx = flatElemIdx - 2*grid.Nx;
      int eeElemIdx = flatElemIdx + 2;
      int wwElemIdx = flatElemIdx - 2;

      // Determine BCs in the current cell
      int globalSouthFaceIdx = grid.faceMap(0,flatElemIdx);
      int globalEastFaceIdx = grid.faceMap(1,flatElemIdx);
      int globalNorthFaceIdx = grid.faceMap(2,flatElemIdx);
      int globalWestFaceIdx = grid.faceMap(3,flatElemIdx);

      bool southBoundary = false;
      bool eastBoundary = false;
      bool northBoundary = false;
      bool westBoundary = false;

      bool ssBoundary = false;
      bool eeBoundary = false;
      bool nnBoundary = false;
      bool wwBoundary = false;
      // BC to the south
      if(grid.BCTypes(globalSouthFaceIdx,0) != 0 || grid.BCTypes(globalSouthFaceIdx,1) != 0)
        southBoundary = true;
      if(grid.BCTypes(globalEastFaceIdx,0) != 0 || grid.BCTypes(globalEastFaceIdx,1) != 0)
        eastBoundary = true;
      if(grid.BCTypes(globalNorthFaceIdx,0) != 0 || grid.BCTypes(globalNorthFaceIdx,1) != 0)
        northBoundary = true;
      if(grid.BCTypes(globalWestFaceIdx,0) != 0 || grid.BCTypes(globalWestFaceIdx,1) != 0)
        westBoundary = true;

      for(int faceNum = 0 ; faceNum < grid.numFacesPerElement; faceNum++)
      {
        // No BCs
        if(grid.cellHasDirichletVelocityBC(flatElemIdx) == 0)
        {
          averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),0) += std::abs(advectionFluxValues[flatElemIdx](faceNum,0)/2.0);
          averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),1) += std::abs(advectionFluxValues[flatElemIdx](faceNum,1)/2.0);
        }
        else
        {

        }
      }
    }
  }

  double signX, signY;
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      signX = (advectionFluxValues[cellNum](faceNum,0) > 0) - (advectionFluxValues[cellNum](faceNum,0) < 0);
      signY = (advectionFluxValues[cellNum](faceNum,1) > 0) - (advectionFluxValues[cellNum](faceNum,1) < 0);
      totalCellAdvectionFlux(cellNum,0) += averagedFaceFluxes(grid.faceMap(faceNum,cellNum),0)*signX;
      totalCellAdvectionFlux(cellNum,1) += averagedFaceFluxes(grid.faceMap(faceNum,cellNum),1)*signY;
    }
  }

  //std::cout<<"Total cell fluxes: "<<std::endl;
  //std::cout<<totalCellAdvectionFlux<<std::endl;
}


void AdvectionUpwind::applyDirichletBCs(CartesianGrid &grid)
{
  // Account for BCs
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    if(grid.cellHasDirichletVelocityBC(cellNum) == 0)
      continue;

    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      // Dirichlet on u-vel
      if(grid.BCTypes(grid.faceMap(faceNum,cellNum),0) > 0)
        uVelInterp(cellNum,faceNum) = grid.uVelBCValue(grid.BCTypes(grid.faceMap(faceNum,cellNum),0)-1);

      // Dirichlet on v-vel
      if(grid.BCTypes(grid.faceMap(faceNum,cellNum),1) > 0)
        vVelInterp(cellNum,faceNum) = grid.vVelBCValue(grid.BCTypes(grid.faceMap(faceNum,cellNum),1)-1);
    }
  }
}

void AdvectionUpwind::clearData()
{
  for(int ii = 0; ii < advectionFluxValues.size(); ii++)
    advectionFluxValues[ii].setZero();
  totalCellAdvectionFlux.setZero();
  uVelInterp.setZero();
  vVelInterp.setZero();
}

