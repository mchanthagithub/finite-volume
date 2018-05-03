//
// Created by maytee on 4/30/18.
//

#include <iostream>
#include "ExplicitPressure.h"

void ExplicitPressure::calculatePressure(CartesianGrid& grid, Eigen::MatrixXd H)
{

  pressure.setZero(grid.numCells);

  double delX = grid.delX;
  double delY = grid.delY;
  double factor = delX*delX*delY*delY/(delX*delX+delY*delY);
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii+kk*grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;
      int northElemIdx = flatElemIdx + grid.Nx;
      int eastElemIdx = flatElemIdx + 1;
      int westElemIdx = flatElemIdx - 1;

      // No pressure BCs on the current cell;
      if(grid.cellHasPressureBC(flatElemIdx) == 0 )
      {
        pressure(flatElemIdx) = factor*( (grid.pressure(eastElemIdx)+grid.pressure(westElemIdx))/(delX*delX) +
                                         (grid.pressure(northElemIdx)+grid.pressure(southElemIdx))/(delY*delY) -
                                         (H(flatElemIdx,0)-H(westElemIdx,0))/delX -
                                         (H(flatElemIdx,1) - H(southElemIdx,1))/delY );
      }
      else
      {
        bool hasSouthBoundary = false;
        bool hasNorthBoundary = false;
        bool hasEastBoundary = false;
        bool hasWestBoundary = false;
        if(grid.BCTypes(grid.faceMap(0,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(0,flatElemIdx),1) != 0)
          hasSouthBoundary = true;
        if(grid.BCTypes(grid.faceMap(1,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(1,flatElemIdx),1) != 0)
          hasEastBoundary = true;
        if(grid.BCTypes(grid.faceMap(2,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(2,flatElemIdx),1) != 0)
          hasNorthBoundary = true;
        if(grid.BCTypes(grid.faceMap(3,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(3,flatElemIdx),1) != 0)
          hasWestBoundary = true;

        if(hasSouthBoundary)
        {
          // Corner
          if(hasEastBoundary || hasWestBoundary)
          {
            pressure(flatElemIdx) = grid.pressure(flatElemIdx);
          }
          else
          {
            pressure(flatElemIdx) = grid.pressure(eastElemIdx) + grid.pressure(westElemIdx) -
                                    delX*delX*( (H(flatElemIdx,0)-H(westElemIdx,0))/delX +
                                              (H(northElemIdx,1)-H(flatElemIdx,1))/delY );
          }
        }
        else if (hasNorthBoundary)
        {
          // Corner
          if(hasEastBoundary || hasWestBoundary)
          {
            pressure(flatElemIdx) = grid.pressure(flatElemIdx);
          }
          else
          {
            pressure(flatElemIdx) = grid.pressure(eastElemIdx) + grid.pressure(westElemIdx) -
                                    delX*delX*( (H(flatElemIdx,0)-H(westElemIdx,0))/delX +
                                    (H(flatElemIdx,1)-H(southElemIdx,1))/delY );
          }
        }

        if(hasEastBoundary)
        {
          // Corner
          if(hasNorthBoundary || hasSouthBoundary)
          {
            pressure(flatElemIdx) = grid.pressure(flatElemIdx);
          }
          else
          {
            pressure(flatElemIdx) = grid.pressure(northElemIdx) + grid.pressure(southElemIdx) -
                                    delY*delY*( (H(flatElemIdx,0)-H(westElemIdx,0))/delX +
                                    (H(flatElemIdx,1)-H(southElemIdx,1))/delY );
          }
        }
        else if(hasWestBoundary)
        {
          // Corner
          if (hasNorthBoundary || hasSouthBoundary)
          {
            pressure(flatElemIdx) = grid.pressure(flatElemIdx);
          }
          else
          {
            pressure(flatElemIdx) = grid.pressure(northElemIdx) + grid.pressure(southElemIdx) -
                                    delY * delY * ((H(eastElemIdx, 0) - H(flatElemIdx, 0)) / delX +
                                                   (H(flatElemIdx, 1) - H(southElemIdx, 1)) / delY);
          }
        }
      }
      std::cout<<"cell: "<<flatElemIdx<<" pressure: "<<pressure(flatElemIdx)<<std::endl;
    }
  }
}

void ExplicitPressure::calculatePressureGradient(CartesianGrid &grid, Eigen::MatrixXd H)
{
  calculatePressure(grid,H);
  pressureGradient.setZero(grid.numCells,grid.nDim);

  double delX = grid.delX;
  double delY = grid.delY;
  double factor = delX*delX*delY*delY/(delX*delX+delY*delY);
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii+kk*grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;
      int northElemIdx = flatElemIdx + grid.Nx;
      int eastElemIdx = flatElemIdx + 1;
      int westElemIdx = flatElemIdx - 1;


      // No pressure BCs on the current cell;
      if(grid.cellHasPressureBC(flatElemIdx) == 0 )
      {
        // Forward difference
        pressureGradient(flatElemIdx,0) = (pressure(eastElemIdx) - pressure(flatElemIdx)) / delX;
        pressureGradient(flatElemIdx,1) = (pressure(northElemIdx) - pressure(flatElemIdx)) / delY;
      }
      else
      {
        bool hasSouthBoundary = false;
        bool hasNorthBoundary = false;
        bool hasEastBoundary = false;
        bool hasWestBoundary = false;
        if(grid.BCTypes(grid.faceMap(0,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(0,flatElemIdx),1) != 0)
          hasSouthBoundary = true;
        if(grid.BCTypes(grid.faceMap(1,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(1,flatElemIdx),1) != 0)
          hasEastBoundary = true;
        if(grid.BCTypes(grid.faceMap(2,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(2,flatElemIdx),1) != 0)
          hasNorthBoundary = true;
        if(grid.BCTypes(grid.faceMap(3,flatElemIdx),0) != 0 || grid.BCTypes(grid.faceMap(3,flatElemIdx),1) != 0)
          hasWestBoundary = true;

        if(hasSouthBoundary)
        {
          // Corner
          if(hasEastBoundary || hasWestBoundary)
          {
            pressureGradient(flatElemIdx,0) = 0;
            pressureGradient(flatElemIdx,1) = 0;
          }
          else
          {
            pressureGradient(flatElemIdx,0) = (pressure(eastElemIdx) - pressure(flatElemIdx))/delX;
            pressureGradient(flatElemIdx,1) = 0;
          }
        }
        else if (hasNorthBoundary)
        {
          // Corner
          if(hasEastBoundary || hasWestBoundary)
          {
            pressureGradient(flatElemIdx,0) = 0;
            pressureGradient(flatElemIdx,1) = 0;
          }
          else
          {
            pressureGradient(flatElemIdx,0) = (pressure(eastElemIdx) - pressure(flatElemIdx))/delX;
            pressureGradient(flatElemIdx,1) = 0;
          }
        }
        if(hasEastBoundary)
        {
          // Corner
          if(hasNorthBoundary || hasSouthBoundary)
          {
            pressureGradient(flatElemIdx,0) = 0;
            pressureGradient(flatElemIdx,1) = 0;
          }
          else
          {
            pressureGradient(flatElemIdx,0) = 0;
            pressureGradient(flatElemIdx,1) = (pressure(northElemIdx) - pressure(flatElemIdx))/delY;
          }
        }
        else if(hasWestBoundary)
        {
          // Corner
          if(hasNorthBoundary || hasSouthBoundary)
          {
            pressureGradient(flatElemIdx,0) = 0;
            pressureGradient(flatElemIdx,1) = 0;
          }
          else
          {
            pressureGradient(flatElemIdx,0) = 0;
            pressureGradient(flatElemIdx,1) = (pressure(northElemIdx) - pressure(flatElemIdx))/delY;
          }
        }
      }
      //std::cout<<"cell: "<<flatElemIdx<<" pressure gradient: "<<pressureGradient.row(flatElemIdx)<<std::endl;
    }
  }

}

void ExplicitPressure::clearData()
{
  pressure.setZero();
  pressureGradient.setZero();

}