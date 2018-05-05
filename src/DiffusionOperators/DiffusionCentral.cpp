//
// Created by maytee on 4/30/18.
//

#include "DiffusionCentral.h"
#include "iostream"
void DiffusionCentral::calculateDiffusionFluxesCartesian(CartesianGrid &grid, InterpolateUpwind& interp)
{
  Eigen::MatrixXd faceNormals;
  faceNormals.setZero(grid.nDim,grid.numFacesPerElement);
  faceNormals(0,0) = 0.0; faceNormals(1,0) = -1.0;
  faceNormals(0,1) = 1.0; faceNormals(1,1) = 0.0;
  faceNormals(0,2) = 0.0; faceNormals(1,2) = 1,0;
  faceNormals(0,3) = -1.0; faceNormals(1,3) = 0.0;
  if(diffusionCentralFluxValues.size() == 0)
    diffusionCentralFluxValues.resize(grid.numCells);

  totalCellDiffusionFlux.setZero(grid.numCells,grid.nDim);

  Eigen::MatrixXd averagedFaceFluxes;
  averagedFaceFluxes.setZero(grid.numFaces,grid.nDim);
  //calculateCentralGradients(grid);

  calculateGradients(grid,interp);
  applyNeumannBCs(grid);
  for(int ii = 0; ii < grid.numCells; ii++)
  {
    diffusionCentralFluxValues[ii].setZero(grid.numFacesPerElement,grid.nDim);
  }

  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    //std::cout<<"cell: "<<cellNum<<std::endl;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      //std::cout<<"Multuplied by normal:"<<std::endl;
      //std::cout<< faceVelocityGradients[cellNum][faceNum]*faceNormals.col(faceNum)<<std::endl;
      //std::cout<<"facenum: "<<faceNum<<std::endl;
      //std::cout<<faceVelocityGradients[cellNum][faceNum]<<std::endl;
    }
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
      //std::cout<< faceVelocityGradients[cellNum][faceNum]*faceNormals.col(faceNum)<<std::endl;

      diffusionCentralFluxValues[cellNum].row(faceNum) += sideLengths(faceNum)*(1.0/cellArea)*faceVelocityGradients[cellNum][faceNum]*faceNormals.col(faceNum);
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
        averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),0) += std::abs(diffusionCentralFluxValues[flatElemIdx](faceNum,0)/2.0);
        averagedFaceFluxes(grid.faceMap(faceNum,flatElemIdx),1) += std::abs(diffusionCentralFluxValues[flatElemIdx](faceNum,1)/2.0);
        // No BCs
        if(grid.cellHasDirichletVelocityBC(flatElemIdx) == 0)
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
      signX = (diffusionCentralFluxValues[cellNum](faceNum,0) > 0) - (diffusionCentralFluxValues[cellNum](faceNum,0) < 0);
      signY = (diffusionCentralFluxValues[cellNum](faceNum,1) > 0) - (diffusionCentralFluxValues[cellNum](faceNum,1) < 0);
      totalCellDiffusionFlux(cellNum,0) += averagedFaceFluxes(grid.faceMap(faceNum,cellNum),0)*signX;
      totalCellDiffusionFlux(cellNum,1) += averagedFaceFluxes(grid.faceMap(faceNum,cellNum),1)*signY;
      //totalCellDiffusionFlux(cellNum,0) += diffusionCentralFluxValues[cellNum](faceNum,0);
      //totalCellDiffusionFlux(cellNum,1) += diffusionCentralFluxValues[cellNum](faceNum,1);
    }
    //std::cout<<"Total x diffusion: "<<cellNum<<" "<<totalCellDiffusionFlux(cellNum,0)<<std::endl;
    //std::cout<<"Total y diffusion: "<<cellNum<<" "<<totalCellDiffusionFlux(cellNum,1)<<std::endl;
  }

}
/*
void DiffusionCentral::calculateGradients(CartesianGrid &grid, InterpolateUpwind& interp)
{

  Eigen::MatrixXd uVelInterp = interp.uVelInterp;
  Eigen::MatrixXd vVelInterp = interp.vVelInterp;

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

      // Determine BCs in the current cell
      int globalSouthFaceIdx = grid.faceMap(0,flatElemIdx);
      int globalEastFaceIdx = grid.faceMap(1,flatElemIdx);
      int globalNorthFaceIdx = grid.faceMap(2,flatElemIdx);
      int globalWestFaceIdx = grid.faceMap(3,flatElemIdx);

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

      for(int faceNum = 0 ; faceNum < grid.numFacesPerElement; faceNum++)
      {
        partialX.setZero(grid.nDim);
        partialY.setZero(grid.nDim);
        // No BCs
        if(grid.cellHasDirichletVelocityBC(flatElemIdx) == 0)
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
        else
        {
          int globalFaceIdx = grid.faceMap(faceNum,flatElemIdx);

          //partialY(0) = 0;
          //partialY(1) = 0;
          //partialX(0) = 0;
          //partialX(1) = 0;

          if(faceNum == 0 || faceNum == 2)
          {
            // Take care of partial X derivatives
            if(eastBoundary)
            {
              partialX(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/delX;
              partialX(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/delX;
            }
            else if(westBoundary)
            {
              partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/delX;
              partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/delX;
            }
            else
            {
              partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/(2*delX);
              partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/(2*delX);
            }

            // Take care of partial Y derivatives
            if(faceNum == 0)
            {
              if(southBoundary)
              {
                partialY(0) = (grid.uVel(flatElemIdx) - uVelInterp(flatElemIdx,faceNum))/(delY/2.0);
                partialY(1) = (grid.vVel(flatElemIdx) - vVelInterp(flatElemIdx,faceNum))/(delY/2.0);
              }
              else
              {
                partialY(0) = (grid.uVel(flatElemIdx) - grid.uVel(southElemIdx,faceNum))/(delY);
                partialY(1) = (grid.vVel(flatElemIdx) - grid.vVel(southElemIdx,faceNum))/(delY);
              }
            }
            else if(faceNum == 2)
            {
              if(northBoundary)
              {
                partialY(0) = (uVelInterp(flatElemIdx,faceNum) - grid.uVel(flatElemIdx))/(delY/2.0);
                partialY(1) = (vVelInterp(flatElemIdx,faceNum) - grid.vVel(flatElemIdx))/(delY/2.0);
              }
              else
              {
                partialY(0) = (grid.uVel(northElemIdx) - grid.uVel(flatElemIdx))/(delY);
                partialY(1) = (grid.vVel(northElemIdx) - grid.vVel(flatElemIdx))/(delY);
              }
            }
          }
          else if(faceNum == 1 || faceNum == 3)
          {
            // Take care of partial Y derivatives
            if(northBoundary)
            {
              partialY(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/delY;
              partialY(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/delY;
            }
            else if(southBoundary)
            {
              partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/delY;
              partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/delY;
            }
            else
            {
              partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/(2*delY);
              partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/(2*delY);
            }

            // Take care of partial X derivatives
            if(faceNum == 1)
            {
              if(eastBoundary)
              {
                partialX(0) = (uVelInterp(flatElemIdx,faceNum) - grid.uVel(flatElemIdx))/(delX/2.0);
                partialX(1) = (vVelInterp(flatElemIdx,faceNum) - grid.vVel(flatElemIdx))/(delX/2.0);
              }
              else
              {
                partialX(0) = (grid.uVel(eastElemIdx) - grid.uVel(flatElemIdx))/delX;
                partialX(1) = (grid.vVel(eastElemIdx) - grid.vVel(flatElemIdx))/delX;
              }
            }
            else if(faceNum == 3)
            {
              if(westBoundary)
              {
                partialX(0) = (grid.uVel(flatElemIdx) - uVelInterp(flatElemIdx,faceNum))/(delX/2.0);
                partialX(1) = (grid.vVel(flatElemIdx) - vVelInterp(flatElemIdx,faceNum))/(delX/2.0);
              }
              else
              {
                partialX(0) = (grid.uVel(flatElemIdx) - grid.uVel(westElemIdx))/delX;
                partialX(1) = (grid.vVel(flatElemIdx) - grid.vVel(westElemIdx))/delX;
              }
            }
          }
        }
        faceVelocityGradients[flatElemIdx][faceNum].col(0) = partialX;
        faceVelocityGradients[flatElemIdx][faceNum].col(1) = partialY;
      }
    }
  }
}
*/
void DiffusionCentral::calculateGradients(CartesianGrid &grid, InterpolateUpwind& interp)
{
  Eigen::MatrixXd uVelInterp = interp.uVelInterp;
  Eigen::MatrixXd vVelInterp = interp.vVelInterp;
  // For each cell calculate gradient based on surounding cells. Do a forward or backward difference for the cross derivatives
  // depending on upwinding
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
      partialX.setZero(grid.nDim);
      partialY.setZero(grid.nDim);
      int flatElemIdx = ii + kk*grid.Nx;
      int northElemIdx = flatElemIdx + grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;
      int eastElemIdx = flatElemIdx + 1;
      int westElemIdx = flatElemIdx - 1;

      // Determine BCs in the current cell
      int globalSouthFaceIdx = grid.faceMap(0,flatElemIdx);
      int globalEastFaceIdx = grid.faceMap(1,flatElemIdx);
      int globalNorthFaceIdx = grid.faceMap(2,flatElemIdx);
      int globalWestFaceIdx = grid.faceMap(3,flatElemIdx);

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

      for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
      {
        int globalFaceIdx = grid.faceMap(faceNum,flatElemIdx);

        //partialY(0) = 0;
        //partialY(1) = 0;
        //partialX(0) = 0;
        //partialX(1) = 0;

        if(faceNum == 0 || faceNum == 2)
        {
          // Take care of partial X derivatives
          if(eastBoundary)
          {
            partialX(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/delX;
            partialX(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/delX;
          }
          else if(westBoundary)
          {
            partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/delX;
            partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/delX;
          }
          else
          {
            if(uVelInterp(eastElemIdx,faceNum) >= uVelInterp(westElemIdx,faceNum))
            {
              partialX(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/(delX);
              partialX(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/(delX);
            }
            else
            {
              partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/(delX);
              partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/(delX);
            }
          }

          // Take care of partial Y derivatives
          if(faceNum == 0)
          {
            if(southBoundary)
            {
              //partialY(0) = (grid.uVel(flatElemIdx) - uVelInterp(flatElemIdx,faceNum))/(delY/2.0);
              //partialY(1) = (grid.vVel(flatElemIdx) - vVelInterp(flatElemIdx,faceNum))/(delY/2.0);

              partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/(delY);
              partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/(delY);
            }
            else
            {
              partialY(0) = (grid.uVel(flatElemIdx) - grid.uVel(southElemIdx,faceNum))/(delY);
              partialY(1) = (grid.vVel(flatElemIdx) - grid.vVel(southElemIdx,faceNum))/(delY);
            }
          }
          else if(faceNum == 2)
          {
            if(northBoundary)
            {
              //partialY(0) = (uVelInterp(flatElemIdx,faceNum) - grid.uVel(flatElemIdx))/(delY/2.0);
              //partialY(1) = (vVelInterp(flatElemIdx,faceNum) - grid.vVel(flatElemIdx))/(delY/2.0);
              partialY(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/(delY);
              partialY(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/(delY);
            }
            else
            {
              partialY(0) = (grid.uVel(northElemIdx) - grid.uVel(flatElemIdx))/(delY);
              partialY(1) = (grid.vVel(northElemIdx) - grid.vVel(flatElemIdx))/(delY);
            }
          }
        }
        else if(faceNum == 1 || faceNum == 3)
        {
          // Take care of partial Y derivatives
          if(northBoundary)
          {
            partialY(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/delY;
            partialY(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/delY;
          }
          else if(southBoundary)
          {
            partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/delY;
            partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/delY;
          }
          else
          {
            if(vVelInterp(northElemIdx,faceNum) >= vVelInterp(southElemIdx,faceNum))
            {
              partialY(0) = (uVelInterp(northElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/(delY);
              partialY(1) = (vVelInterp(northElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/(delY);
            }
            else
            {
              partialY(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(southElemIdx,faceNum))/(delY);
              partialY(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(southElemIdx,faceNum))/(delY);
            }
          }

          // Take care of partial X derivatives
          if(faceNum == 1)
          {
            if(eastBoundary)
            {
              //partialX(0) = (uVelInterp(flatElemIdx,faceNum) - grid.uVel(flatElemIdx))/(delX/2.0);
              //partialX(1) = (vVelInterp(flatElemIdx,faceNum) - grid.vVel(flatElemIdx))/(delX/2.0);
              partialX(0) = (uVelInterp(flatElemIdx,faceNum) - uVelInterp(westElemIdx,faceNum))/(delX);
              partialX(1) = (vVelInterp(flatElemIdx,faceNum) - vVelInterp(westElemIdx,faceNum))/(delX);
            }
            else
            {
              partialX(0) = (grid.uVel(eastElemIdx) - grid.uVel(flatElemIdx))/delX;
              partialX(1) = (grid.vVel(eastElemIdx) - grid.vVel(flatElemIdx))/delX;
            }
          }
          else if(faceNum == 3)
          {
            if(westBoundary)
            {
              //partialX(0) = (grid.uVel(flatElemIdx) - uVelInterp(flatElemIdx,faceNum))/(delX/2.0);
              //partialX(1) = (grid.vVel(flatElemIdx) - vVelInterp(flatElemIdx,faceNum))/(delX/2.0);

              partialX(0) = (uVelInterp(eastElemIdx,faceNum) - uVelInterp(flatElemIdx,faceNum))/(delX);
              partialX(1) = (vVelInterp(eastElemIdx,faceNum) - vVelInterp(flatElemIdx,faceNum))/(delX);
            }
            else
            {
              partialX(0) = (grid.uVel(flatElemIdx) - grid.uVel(westElemIdx))/delX;
              partialX(1) = (grid.vVel(flatElemIdx) - grid.vVel(westElemIdx))/delX;
            }
          }
        }
        faceVelocityGradients[flatElemIdx][faceNum].col(0) = partialX;
        faceVelocityGradients[flatElemIdx][faceNum].col(1) = partialY;
      }
    }
  }
}

void DiffusionCentral::clearData()
{
  for(int ii = 0; ii < diffusionCentralFluxValues.size(); ii++)
  {
    diffusionCentralFluxValues[ii].setZero();
  }
  for(int ii = 0; ii < centralGradients.size(); ii++)
    centralGradients[ii].setZero();

  for(int ii = 0; ii < faceVelocityGradients.size(); ii++)
    for(int kk = 0; kk <faceVelocityGradients[ii].size();kk++)
      faceVelocityGradients[ii][kk].setZero();

  totalCellDiffusionFlux.setZero();
}

void DiffusionCentral::applyNeumannBCs(CartesianGrid &grid)
{
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      int globalFaceIdx = grid.faceMap(faceNum,cellNum);
      if(grid.faceHasNeumannVelocityBC(globalFaceIdx) > 0)
      {
        for(int kk = 0; kk < grid.nDim; kk++)
        {
          for(int ii = 0; ii < grid.nDim; ii++)
          {
            if(grid.faceVelocityNeumannBCs[globalFaceIdx](ii,kk) > 0)
            {
              faceVelocityGradients[cellNum][faceNum](ii,kk) = grid.uVelNeumannBCValue(grid.faceVelocityNeumannBCs[globalFaceIdx](ii,kk)-1);
            }
          }
        }
      }
    }
  }
}