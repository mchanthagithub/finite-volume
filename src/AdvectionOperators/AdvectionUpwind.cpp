//
// Created by maytee on 4/27/18.
//

#include "AdvectionUpwind.h"
#include <iostream>
#include "cstdlib"
void AdvectionUpwind::calculateAdvectionFluxesCartesian(CartesianGrid &grid, InterpolateUpwind& interp)
{
  double rho = 1000;
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
        averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),0) += std::abs(advectionFluxValues[flatElemIdx](faceNum,0)/2.0);
        averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),1) += std::abs(advectionFluxValues[flatElemIdx](faceNum,1)/2.0);
        // No BCs
        if(grid.cellHasVelocityBC(flatElemIdx) == 0)
        {
          //averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),0) += std::abs(advectionFluxValues[flatElemIdx](faceNum,0)/2.0);
          //averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),1) += std::abs(advectionFluxValues[flatElemIdx](faceNum,1)/2.0);
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

void AdvectionUpwind::calculateAdvectionFluxesCartesian(CartesianGrid &grid, InterpolateQUICK& interp)
{
  double rho = 1000;
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
        if(grid.cellHasVelocityBC(flatElemIdx) == 0)
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


void AdvectionUpwind::calculateInterpolationValues(CartesianGrid &grid)
{
  // Go cell by cell and figure out whether the velocity is upwind or not in both the x and y directions
  // For east face, check element east, for west face, check element west, for south face check element south, etc
  // At boundary of boubdary elements, eg the west face of a left-side cell, just take the current cell center value,
  // though still check east element for east face
  uVelInterp.setZero(grid.numCells,grid.numFacesPerElement);
  vVelInterp.setZero(grid.numCells,grid.numFacesPerElement);

  // Loop through all elements going left to right, bottom to top
  // Check all faces of each element and account for boundaries at the end
  // Boundary values will overwrite whatever is calculated here, but still calculate for simplicity
  bool isTopElement = false;
  bool isBotElement = false;
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    isTopElement = (kk == grid.Ny-1);
    isBotElement = (kk == 0);
    bool isRightElement = false;
    bool isLeftElement = false;
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      isLeftElement = (ii == 0);
      isRightElement = (ii == grid.Nx-1);
      int flatElemIdx = ii + kk*grid.Nx;
      int elemNorthIdx = flatElemIdx + grid.Nx;
      int elemSouthIdx = flatElemIdx - grid.Nx;
      int elemEastIdx = flatElemIdx + 1;
      int elemWestIdx = flatElemIdx - 1;

      // Bottom edge of grid or top edge of grid, just set the velocity vector equal to the current center node value
      if(isBotElement)
      {
        vVelInterp(flatElemIdx,0) = grid.vVel(flatElemIdx);
        uVelInterp(flatElemIdx,0) = grid.uVel(flatElemIdx);
      }
      else if(isTopElement)
      {
        vVelInterp(flatElemIdx,2) = grid.vVel(flatElemIdx);
        uVelInterp(flatElemIdx,2) = grid.uVel(flatElemIdx);
      }

      // Left edge of grid or right edge of grid, just set velocity vector equal to the current center node value
      if(isLeftElement)
      {
        uVelInterp(flatElemIdx,3) = grid.uVel(flatElemIdx);
        vVelInterp(flatElemIdx,3) = grid.vVel(flatElemIdx);
      }
      else if(isRightElement)
      {
        uVelInterp(flatElemIdx,1) = grid.uVel(flatElemIdx);
        vVelInterp(flatElemIdx,1) = grid.vVel(flatElemIdx);
      }

      if(!isBotElement)
      {
        if(!isTopElement)
        {
          // Check south face
          if(grid.vVel(elemSouthIdx) >= 0)
          {
            vVelInterp(flatElemIdx,0) = grid.vVel(elemSouthIdx);
            uVelInterp(flatElemIdx,0) = grid.uVel(elemSouthIdx);
          }
          else
          {
            vVelInterp(flatElemIdx,0) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,0) = grid.uVel(flatElemIdx);
          }

          // Check north face
          if(grid.vVel(elemNorthIdx) >= 0)
          {
            vVelInterp(flatElemIdx,2) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,2) = grid.uVel(flatElemIdx);
          }
          else
          {
            vVelInterp(flatElemIdx,2) = grid.vVel(elemNorthIdx);
            uVelInterp(flatElemIdx,2) = grid.uVel(elemNorthIdx);
          }
        }
        else
        {
          // Take care of south face of top element
          if(grid.vVel(elemSouthIdx) >= 0)
          {
            vVelInterp(flatElemIdx,0) = grid.vVel(elemSouthIdx);
            uVelInterp(flatElemIdx,0) = grid.uVel(elemSouthIdx);
          }
          else
          {
            vVelInterp(flatElemIdx,0) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,0) = grid.uVel(flatElemIdx);
          }
        }
      }
      else
      {
        // Take care of north face of bottom elements
        if(grid.vVel(elemNorthIdx) >= 0)
        {
          vVelInterp(flatElemIdx,2) = grid.vVel(flatElemIdx);
          uVelInterp(flatElemIdx,2) = grid.uVel(flatElemIdx);
        }
        else
        {
          vVelInterp(flatElemIdx,2) = grid.vVel(elemNorthIdx);
          uVelInterp(flatElemIdx,2) = grid.uVel(elemNorthIdx);
        }
      }

      if(!isLeftElement)
      {
        if(!isRightElement)
        {
          // Check west face
          if(grid.uVel(elemWestIdx) >= 0)
          {
            uVelInterp(flatElemIdx,3) = grid.uVel(elemWestIdx);
            vVelInterp(flatElemIdx,3) = grid.vVel(elemWestIdx);
          }
          else
          {
            uVelInterp(flatElemIdx,3) = grid.uVel(flatElemIdx);
            vVelInterp(flatElemIdx,3) = grid.vVel(flatElemIdx);
          }

          // Check east face
          if(grid.uVel(elemEastIdx) >= 0)
          {
            uVelInterp(flatElemIdx,1) = grid.uVel(flatElemIdx);
            vVelInterp(flatElemIdx,1) = grid.vVel(flatElemIdx);
          }
          else
          {
            uVelInterp(flatElemIdx,1) = grid.uVel(elemEastIdx);
            vVelInterp(flatElemIdx,1) = grid.vVel(elemEastIdx);
          }
        }
        else
        {
          // Take care of west face of right element
          if(grid.uVel(elemWestIdx) >= 0)
          {
            uVelInterp(flatElemIdx,3) = grid.uVel(elemWestIdx);
            vVelInterp(flatElemIdx,3) = grid.vVel(elemWestIdx);
          }
          else
          {
            uVelInterp(flatElemIdx,3) = grid.uVel(flatElemIdx);
            vVelInterp(flatElemIdx,3) = grid.vVel(flatElemIdx);
          }
        }
      }
      else
      {
        // Take care of east face of left element
        if(grid.uVel(elemEastIdx) >= 0)
        {
          uVelInterp(flatElemIdx,1) = grid.uVel(flatElemIdx);
          vVelInterp(flatElemIdx,1) = grid.vVel(flatElemIdx);
        }
        else
        {
          uVelInterp(flatElemIdx,1) = grid.uVel(elemEastIdx);
          vVelInterp(flatElemIdx,1) = grid.vVel(elemEastIdx);
        }
      }
    }
  }

  // Now apply BCs
  applyDirichletBCs(grid);

}


void AdvectionUpwind::applyDirichletBCs(CartesianGrid &grid)
{
  // Account for BCs
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    if(grid.cellHasVelocityBC(cellNum) == 0)
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

