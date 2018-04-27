//
// Created by maytee on 4/27/18.
//

#include "AdvectionUpwind.h"
#include <iostream>
void AdvectionUpwind::calculateAdvectionFluxesCartesian(CartesianGrid &grid)
{
  double rho = 1000;
  // Calculate advection flux at all surfaces for all cells, then average to obtain final flux term
  // Here use midpoint rule (fluxMid*surfaceArea = totalFlux + O(delX^2))

  calculateInterpolationValues(grid);
  std::cout<<"uVelInterp: "<<std::endl;
  std::cout<<uVelInterp.transpose()<<std::endl;

  std::cout<<"vVelInterp: "<<std::endl;
  std::cout<<vVelInterp.transpose()<<std::endl;

  advectionFluxValues.setZero(grid.numCells);
  Eigen::MatrixXd localFluxes;
  localFluxes.setZero(grid.numCells,grid.numFacesPerElement);

  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    // Need to check for boundary conditions
    if(grid.cellBCTypes[cellNum](0,0) > 0)
    {

    }
  }


}

void AdvectionUpwind::calculateInterpolationValues(CartesianGrid &grid)
{
  // Go cell by cell and figure out whether the velocity is upwind or not in both the x and y directions
  // Take convention that will take cell center value to the right, and top. of each element to check
  // if the velocity is positive (away from current cell and so upwind is P) or if velocity is towards cell
  // (and so upwind is E or N) except for rightmost and topmost elements, then use element directly left and bottom
  uVelInterp.setZero(grid.numCells);
  vVelInterp.setZero(grid.numCells);

  bool isTop = false;
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    isTop = (kk == grid.Ny-1);
    bool isRight = false;
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      isRight = (ii == grid.Nx-1);
      int flatElemIdx = ii + kk*grid.Nx;

      if(isRight)
      {
        if(grid.uVel(flatElemIdx-1) > 0 )
          uVelInterp(flatElemIdx) = grid.uVel(flatElemIdx-1);
        else
          uVelInterp(flatElemIdx) = grid.uVel(flatElemIdx);
      }
      else
      {
        if(grid.uVel(flatElemIdx+1) > 0)
          uVelInterp(flatElemIdx) = grid.uVel(flatElemIdx);
        else
          uVelInterp(flatElemIdx) = grid.uVel(flatElemIdx+1);
      }

      if(isTop)
      {
        if(grid.vVel(flatElemIdx-grid.Nx) > 0)
          vVelInterp(flatElemIdx) = grid.vVel(flatElemIdx - grid.Nx);
        else
          vVelInterp(flatElemIdx) = grid.vVel(flatElemIdx);
      }
      else
      {
        if(grid.vVel(flatElemIdx+grid.Nx) > 0)
          vVelInterp(flatElemIdx) = grid.vVel(flatElemIdx);
        else
          vVelInterp(flatElemIdx) = grid.vVel(flatElemIdx + grid.Nx);
      }
    }
  }

}
