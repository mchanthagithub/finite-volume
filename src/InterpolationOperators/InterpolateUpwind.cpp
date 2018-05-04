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
  applyDirichletBCs(grid);

  /*
  for(int ii = 0; ii < grid.numCells; ii++)
  {
    std::cout<<"cell: "<<ii<<std::endl;
    std::cout<<uVelInterp.row(ii)<<std::endl;
    std::cout<<vVelInterp.row(ii)<<std::endl;

  }
   */
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
