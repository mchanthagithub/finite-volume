//
// Created by maytee on 4/30/18.
//

#include "DiffusionCentral.h"
#include "iostream"
void DiffusionCentral::calculateDiffusionFluxesCartesian(CartesianGrid &grid)
{
  double nu = 0.001;
  Eigen::MatrixXd faceNormals;
  faceNormals.setZero(grid.nDim,grid.numFacesPerElement);
  faceNormals(0,0) = 0.0; faceNormals(1,0) = -1.0;
  faceNormals(0,1) = 1.0; faceNormals(1,1) = 0.0;
  faceNormals(0,2) = 0.0; faceNormals(1,2) = 1,0;
  faceNormals(0,3) = -1.0; faceNormals(1,3) = 0.0;
  if(diffusionCentralFluxValues.size() == 0)
    diffusionCentralFluxValues.resize(grid.numCells);

  totalCellDiffusionFlux.setZero(grid.numCells,grid.nDim);

  calculateCentralGradients(grid);

  for(int ii = 0; ii < grid.numCells; ii++)
  {
    diffusionCentralFluxValues[ii].setZero(grid.numFacesPerElement,grid.nDim);
  }

  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    double southLength = grid.delX;
    double northLength = grid.delX;
    double eastLength = grid.delY;
    double westLength = grid.delY;

    Eigen::VectorXd sideLengths;
    sideLengths.setZero(grid.numFacesPerElement);
    sideLengths(0) = southLength; sideLengths(1) = eastLength; sideLengths(2) = northLength; sideLengths(3) = westLength;
    double cellArea = southLength*eastLength;
    //std::cout<<"cell gradient: "<<cellNum<<std::endl;
    //std::cout<<centralGradients[cellNum]<<std::endl;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      //std::cout<<"Multuplied by normal:"<<std::endl;
      //std::cout<<centralGradients[cellNum]*faceNormals.col(faceNum)<<std::endl;

      diffusionCentralFluxValues[cellNum].row(faceNum) += sideLengths(faceNum)*(nu/cellArea)*centralGradients[cellNum]*faceNormals.col(faceNum);
    }
  }


  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    std::cout<<"====cell: "<<cellNum<<std::endl;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      std::cout<<"facenum: "<<faceNum<<" "<<diffusionCentralFluxValues[cellNum].row(faceNum)<<std::endl;
      totalCellDiffusionFlux.row(cellNum) += diffusionCentralFluxValues[cellNum].row(faceNum);
    }
    std::cout<<"totalDiffusion: "<<" "<<totalCellDiffusionFlux.row(cellNum)<<std::endl;
  }

}

void DiffusionCentral::calculateGradients(CartesianGrid &grid, InterpolateUpwind& interp)
{

  Eigen::VectorXd uVelInterp = interp.uVelInterp;
  Eigen::VectorXd vVelInterp = interp.vVelInterp;

  // For each cell calculate gradient based on surounding cells. If a boundary cell then do a forward/backward difference
  // to calculate the gradient
  if(faceVelocityGradients.size() == 0)
    faceVelocityGradients.resize(grid.numCells);

  for(int ii = 0; ii < grid.numCells; ii++)
    faceVelocityGradients[ii].resize(grid.numFacesPerElement);

  for(int ii = 0 ; ii <grid.numCells; ii++)
    for(int kk = 0; kk < grid.numFacesPerElement; kk++)
      faceVelocityGradients[ii][kk].setZero(grid.nDim,grid.nDim);

  // Part of gradient that is derivatives with respect to y. partialY(0) is delu/delY, partialY(1) is delv/delY
  Eigen::VectorXd partialY;

  // Part of gradient that is derivatives with respect to x. partialX(0) is delu/delX, partialX(2) is delv/delX
  Eigen::VectorXd partialX;
  // Loop through all elements going left to right, bottom to top
  // Check all faces of each element and account for boundaries at the end
  // Boundary values will overwrite whatever is calculated here, but still calculate for simplicity
  double delY = grid.delY;
  double delX = grid.delX;
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + kk*grid.Nx;
      int northElemIdx = flatElemIdx + grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;
      int eastElemIdx = flatElemIdx + 1;
      int westElemIdx = flatElemIdx - 1;

      for(int faceNum = 0 ; faceNum < grid.numFacesPerElement; faceNum++)
      {
        // No BCs
        if(grid.cellHasVelocityBC(flatElemIdx) == 0)
        {
          if(faceNum == 0)
          {
            partialY(0) = (grid.uVel(flatElemIdx) - grid.uVel(southElemIdx))/delY;
            partialY(1) = (grid.vVel(flatElemIdx) - grid.vVel(southElemIdx))/delY;

            // Use interpolated values at center of cell faces
            partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/(2*delX);
            partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/(2*delX);
          }
          else if(faceNum == 1)
          {
            partialX(0) = (grid.uVel(eastElemIdx) - grid.uVel(flatElemIdx))/delX;
            partialX(1) = (grid.vVel(eastElemIdx) - grid.vVel(flatElemIdx))/delX;

            // Use interpolated values
            partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/(2*delY);
            partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/(2*delX);
          }
          else if(faceNum == 2)
          {
            partialY(0) = (grid.uVel(northElemIdx) - grid.uVel(flatElemIdx))/delY;
            partialY(1) = (grid.vVel(northElemIdx) - grid.vVel(flatElemIdx))/delY;

            // Use interpolated values at center of cell faces
            partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/(2*delX);
            partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/(2*delX);
          }
          else if(faceNum == 3)
          {
            partialX(0) = (grid.uVel(flatElemIdx) - grid.uVel(westElemIdx))/delX;
            partialX(1) = (grid.vVel(flatElemIdx) - grid.vVel(westElemIdx))/delX;

            // Use interpolated values
            partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/(2*delY);
            partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/(2*delX);
          }
        }
        faceVelocityGradients[flatElemIdx][faceNum].col(0) = partialX;
        faceVelocityGradients[flatElemIdx][faceNum].col(1) = partialY;
      }
      // Deal with BCs
      int globalSouthFaceIdx = grid.faceMap(flatElemIdx,0);
      int globalEastFaceIdx = grid.faceMap(flatElemIdx,1);
      int globalNorthFaceIdx = grid.faceMap(flatElemIdx,2);
      int globalWestFaceIdx = grid.faceMap(flatElemIdx,3);

      bool southBoundary = false;
      bool eastBoundary = false;
      bool northBoundary = false;
      bool westBoundary = false;
      // BC to the south
      if(grid.BCTypes(globalSouthFaceIdx,0) != 0 || grid.BCTypes(globalSouthFaceIdx,1) != 0)
        southBoundary = true;
      if(grid.BCTypes(globalEastFaceIdx,0) != 0 || grid.BCTypes(globalEastFaceIdx,1) != 0)
        eastBoundary = true;
      if(grid.BCTypes(globalNorthFaceIdx,0) != 0 || grid.BCTypes(globalNorthFaceIdx,1) != 0)
        northBoundary = true;
      if(grid.BCTypes(globalWestFaceIdx,0) != 0 || grid.BCTypes(globalWestFaceIdx,1) != 0)
        westBoundary = true;


    }
  }
}

void DiffusionCentral::calculateCentralGradients(CartesianGrid &grid)
{
  // For each cell calculate gradient based on surounding cells. If a boundary cell then do a forward/backward difference
  // to calculate the gradient
  if(centralGradients.size() == 0)
    centralGradients.resize(grid.numCells);

  for(int ii = 0; ii < grid.numCells; ii++)
    centralGradients[ii].setZero(grid.nDim,grid.nDim);

  // Part of gradient that is derivatives with respect to y. partialY(0) is delu/delY, partialY(1) is delv/delY
  Eigen::VectorXd partialY;

  // Part of gradient that is derivatives with respect to x. partialX(0) is delu/delX, partialX(2) is delv/delX
  Eigen::VectorXd partialX;
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + kk*grid.Nx;
      int eastElemIdx = flatElemIdx + 1;
      int westElemIdx = flatElemIdx - 1;
      int northElemIdx = flatElemIdx + grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;

      partialY.setZero(grid.nDim);
      partialX.setZero(grid.nDim);
      // Left boundary
      if(ii == 0)
      {
        partialX(0) = (grid.uVel(eastElemIdx) - grid.uVel(flatElemIdx))/grid.delX;
        partialX(1) = (grid.vVel(eastElemIdx) - grid.vVel(flatElemIdx))/grid.delX;
      }
      // Right boundary
      else if(ii == grid.Nx-1)
      {
        partialX(0) = (grid.uVel(flatElemIdx) - grid.uVel(westElemIdx))/grid.delX;
        partialX(1) = (grid.vVel(flatElemIdx) - grid.vVel(westElemIdx))/grid.delX;
      }
      // Interior points, just use central difference with elements to left and right
      else
      {
        partialX(0) = (grid.uVel(eastElemIdx) - grid.uVel(westElemIdx))/(2*grid.delX);
        partialX(1) = (grid.vVel(westElemIdx) - grid.vVel(westElemIdx))/(2*grid.delX);
      }

      // Bottom boundary
      if(kk == 0)
      {
        partialY(0) = (grid.uVel(northElemIdx) - grid.uVel(flatElemIdx))/grid.delY;
        partialY(1) = (grid.vVel(northElemIdx) - grid.vVel(flatElemIdx))/grid.delY;
      }
      // Top boundary
      else if(kk == grid.Ny-1)
      {
        partialY(0) = (grid.uVel(flatElemIdx) - grid.uVel(southElemIdx))/grid.delY;
        partialY(1) = (grid.vVel(flatElemIdx) - grid.vVel(southElemIdx))/grid.delY;
      }
      else
      {
        partialY(0) = (grid.uVel(northElemIdx) - grid.uVel(southElemIdx))/(2*grid.delY);
        partialY(1) = (grid.vVel(northElemIdx) - grid.vVel(southElemIdx))/(2*grid.delY);
      }

      centralGradients[flatElemIdx].col(0) = partialX;
      centralGradients[flatElemIdx].col(1) = partialY;
    }
  }
}

void DiffusionCentral::clearData()
{
  for(int ii = 0; ii < diffusionCentralFluxValues.size(); ii++)
  {
    diffusionCentralFluxValues[ii].setZero();
    centralGradients[ii].setZero();
  }
  totalCellDiffusionFlux.setZero();
}