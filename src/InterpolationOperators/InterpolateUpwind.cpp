//
// Created by maytee on 5/1/18.
//

#include "InterpolateUpwind.h"
#include "iostream"
void InterpolateUpwind::interpolateVelocities(CartesianGrid &grid)
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
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + kk*grid.Nx;
      int elemNorthIdx = flatElemIdx + grid.Nx;
      int elemSouthIdx = flatElemIdx - grid.Nx;
      int elemEastIdx = flatElemIdx + 1;
      int elemWestIdx = flatElemIdx - 1;

      for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
      {
        if(faceNum == 0)
        {
          if(kk == 0)
          {
            // Bottom edge of grid or top edge of grid, just set the velocity vector equal to the current center node value
            vVelInterp(flatElemIdx,0) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,0) = grid.uVel(flatElemIdx);
          }
          else
          {
            // Check south face
            //if(grid.vVel(elemSouthIdx) >= grid.vVel(flatElemIdx))
            if(grid.vVel(flatElemIdx) >= 0.0)
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
        else if(faceNum == 1)
        {
          if(ii == grid.Nx-1)
          {
            // Right edge of grid just set the velocity vector equal to the current center node value
            vVelInterp(flatElemIdx,1) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,1) = grid.uVel(flatElemIdx);
          }
          else
          {
            if(flatElemIdx == 0)
            {
              //std::cout<<"current uvel: "<<grid.uVel(flatElemIdx);
              //std::cout<<" east uvel: "<<grid.uVel(elemEastIdx)<<std::endl;
            }
            // Check east face
            //if(grid.uVel(elemEastIdx) >= grid.uVel(flatElemIdx))
            if(grid.uVel(elemEastIdx) >= 0.0)
            {
              //if(flatElemIdx == 0)
              //  std::cout<<"take current"<<std::endl;
              //vVelInterp(flatElemIdx,1) = grid.vVel(elemEastIdx);
              //uVelInterp(flatElemIdx,1) = grid.uVel(elemEastIdx);

              vVelInterp(flatElemIdx,1) = grid.vVel(flatElemIdx);
              uVelInterp(flatElemIdx,1) = grid.uVel(flatElemIdx);
            }
            else
            {
              //if(flatElemIdx == 0)
              //  std::cout<<"take east"<<std::endl;
              //vVelInterp(flatElemIdx,1) = grid.vVel(flatElemIdx);
              //uVelInterp(flatElemIdx,1) = grid.uVel(flatElemIdx);
              vVelInterp(flatElemIdx,1) = grid.vVel(elemEastIdx);
              uVelInterp(flatElemIdx,1) = grid.uVel(elemEastIdx);
            }
          }
        }
        else if(faceNum == 2)
        {
          if(kk == grid.Ny-1)
          {
            // Top edge of grid, just set the velocity vector equal to the current center node value
            vVelInterp(flatElemIdx,2) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,2) = grid.uVel(flatElemIdx);
          }
          else
          {
            // Check north face
            //if(grid.vVel(elemNorthIdx) >= grid.vVel(flatElemIdx))
            if(grid.vVel(elemNorthIdx) >= 0.0)
            {
              //vVelInterp(flatElemIdx,2) = grid.vVel(elemNorthIdx);
              //uVelInterp(flatElemIdx,2) = grid.uVel(elemNorthIdx);

              vVelInterp(flatElemIdx,2) = grid.vVel(flatElemIdx);
              uVelInterp(flatElemIdx,2) = grid.uVel(flatElemIdx);
            }
            else
            {
              //vVelInterp(flatElemIdx,2) = grid.vVel(flatElemIdx);
              //uVelInterp(flatElemIdx,2) = grid.uVel(flatElemIdx);
              vVelInterp(flatElemIdx,2) = grid.vVel(elemNorthIdx);
              uVelInterp(flatElemIdx,2) = grid.uVel(elemNorthIdx);
            }
          }
        }
        else if(faceNum == 3)
        {
          if(ii == 0)
          {
            // Left of grid just set the velocity vector equal to the current center node value
            vVelInterp(flatElemIdx,3) = grid.vVel(flatElemIdx);
            uVelInterp(flatElemIdx,3) = grid.uVel(flatElemIdx);
          }
          else
          {
            // Check west face
            //if(grid.uVel(elemWestIdx) >= grid.uVel(flatElemIdx))
            if(grid.uVel(flatElemIdx) >= 0.0)
            {
              vVelInterp(flatElemIdx,3) = grid.vVel(elemWestIdx);
              uVelInterp(flatElemIdx,3) = grid.uVel(elemWestIdx);
            }
            else
            {
              vVelInterp(flatElemIdx,3) = grid.vVel(flatElemIdx);
              uVelInterp(flatElemIdx,3) = grid.uVel(flatElemIdx);
            }
          }
        }
      }
    }
  }

  applyDirichletBCs(grid);
}

void InterpolateUpwind::clearData()
{
  uVelInterp.setZero();
  vVelInterp.setZero();
}
void InterpolateUpwind::applyDirichletBCs(CartesianGrid &grid)
{
  // Account for BCs
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    if(grid.cellHasDirichletVelocityBC(cellNum) <= 0)
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
