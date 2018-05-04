//
// Created by maytee on 5/1/18.
//

#include "InterpolateQUICK.h"
#include "iostream"
void InterpolateQUICK::interpolateVelocities(CartesianGrid& grid)
{
  uVelInterp.setZero(grid.numCells,grid.numFacesPerElement);
  vVelInterp.setZero(grid.numCells,grid.numFacesPerElement);
  double uCoeff = 6.0/8.0;
  double dCoeff = 3.0/8.0;
  double uuCoeff = 1.0/8.0;
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

      if(southElemIdx < 0 || grid.BCTypes(grid.faceMap(0,southElemIdx),0) != 0 ||  grid.BCTypes(grid.faceMap(0,southElemIdx),1) != 0 )
        ssBoundary = true;
      if(northElemIdx >= grid.Nx*grid.Ny || grid.BCTypes(grid.faceMap(2,northElemIdx),0) != 0 ||  grid.BCTypes(grid.faceMap(2,northElemIdx),1) != 0 )
        nnBoundary = true;
      if(ii == grid.Nx-1 || grid.BCTypes(grid.faceMap(1,eastElemIdx),0) != 0 ||  grid.BCTypes(grid.faceMap(1,eastElemIdx),1) != 0 )
        eeBoundary = true;
      if(ii == 0 || grid.BCTypes(grid.faceMap(3,westElemIdx),0) != 0 ||  grid.BCTypes(grid.faceMap(3,westElemIdx),1) != 0 )
        wwBoundary = true;

      bool surroundingHasBCs = (ssBoundary || nnBoundary || eeBoundary || wwBoundary);


      for(int faceNum = 0 ; faceNum < grid.numFacesPerElement; faceNum++)
      {
        // No BCs
        if(grid.cellHasDirichletVelocityBC(flatElemIdx) == 0 && !surroundingHasBCs)
        {
          if(faceNum == 0)
          {
            if(grid.vVel(southElemIdx) >= grid.vVel(flatElemIdx))
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(southElemIdx)+dCoeff*grid.uVel(flatElemIdx)-uuCoeff*grid.uVel(ssElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(southElemIdx)+dCoeff*grid.vVel(flatElemIdx)-uuCoeff*grid.vVel(ssElemIdx);
            }
            else
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(flatElemIdx)+dCoeff*grid.uVel(southElemIdx)-uuCoeff*grid.uVel(northElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(flatElemIdx)+dCoeff*grid.vVel(southElemIdx)-uuCoeff*grid.vVel(northElemIdx);
            }
          }
          else if(faceNum == 1)
          {
            if(grid.uVel(eastElemIdx) >= grid.uVel(flatElemIdx))
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(flatElemIdx)+dCoeff*grid.uVel(eastElemIdx)-uuCoeff*grid.uVel(westElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(flatElemIdx)+dCoeff*grid.vVel(eastElemIdx)-uuCoeff*grid.vVel(westElemIdx);
            }
            else
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(eastElemIdx)+dCoeff*grid.uVel(flatElemIdx)-uuCoeff*grid.uVel(eeElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(eastElemIdx)+dCoeff*grid.vVel(flatElemIdx)-uuCoeff*grid.vVel(eeElemIdx);
            }
          }
          else if(faceNum == 2)
          {
            if(grid.vVel(northElemIdx) >= grid.vVel(flatElemIdx))
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(flatElemIdx)+dCoeff*grid.uVel(northElemIdx)-uuCoeff*grid.uVel(southElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(flatElemIdx)+dCoeff*grid.vVel(northElemIdx)-uuCoeff*grid.vVel(southElemIdx);
            }
            else
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(northElemIdx)+dCoeff*grid.uVel(flatElemIdx)-uuCoeff*grid.uVel(nnElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(northElemIdx)+dCoeff*grid.vVel(flatElemIdx)-uuCoeff*grid.vVel(nnElemIdx);
            }
          }
          else if(faceNum == 3)
          {
            if(grid.uVel(westElemIdx) >= grid.uVel(flatElemIdx))
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(westElemIdx)+dCoeff*grid.uVel(flatElemIdx)-uuCoeff*grid.uVel(wwElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(westElemIdx)+dCoeff*grid.vVel(flatElemIdx)-uuCoeff*grid.vVel(wwElemIdx);
            }
            else
            {
              uVelInterp(flatElemIdx,faceNum) = uCoeff*grid.uVel(flatElemIdx)+dCoeff*grid.uVel(westElemIdx)-uuCoeff*grid.uVel(eastElemIdx);
              vVelInterp(flatElemIdx,faceNum) = uCoeff*grid.vVel(flatElemIdx)+dCoeff*grid.vVel(westElemIdx)-uuCoeff*grid.vVel(eastElemIdx);
            }
          }
          //std::cout<<"u: "<<uVelInterp(flatElemIdx,faceNum)<<" v: "<<vVelInterp(flatElemIdx,faceNum)<<std::endl;
        }
        else
        {
          /*
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
           */
        }
      }
    }
  }
}

void InterpolateQUICK::clearData()
{
  uVelInterp.setZero();
  vVelInterp.setZero();
}

void InterpolateQUICK::applyDirichletBCs(CartesianGrid& grid)
{
;
}
