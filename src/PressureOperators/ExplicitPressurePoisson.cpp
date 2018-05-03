//
// Created by maytee on 5/2/18.
//

#include "ExplicitPressurePoisson.h"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "Eigen/SparseCholesky"
#include "iostream"
void ExplicitPressurePoisson::calculatePressure(CartesianGrid &grid, Eigen::MatrixXd H)
{
  createMappings(grid);
  // Create matrix to hold central difference operator
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(totalDOF,totalDOF);
  Eigen::VectorXd RHS = Eigen::VectorXd::Zero(totalDOF);
  double delX = grid.delX;
  double delY = grid.delY;
  double xFactor = 1.0/(grid.delX*grid.delX);
  double yFactor = 1.0/(grid.delY*grid.delY);

  double diagonalFactor = (grid.delY*grid.delY+grid.delX*grid.delX)/(grid.delX*grid.delX*grid.delY*grid.delY);
  for(int kk = 0; kk < grid.Ny; kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + kk*grid.Nx;
      int northElemIdx = flatElemIdx + grid.Nx;
      int southElemIdx = flatElemIdx - grid.Nx;
      int eastElemIdx =flatElemIdx + 1;
      int westElemIdx =flatElemIdx - 1;


      // Is masked out cell
      if(mappingGlobalToActive(flatElemIdx) == -2)
        continue;
      else if(mappingGlobalToActive(flatElemIdx) >= 0)
      {
        int activeIdx = mappingGlobalToActive(flatElemIdx);
        int eastActiveIdx, westActiveIdx, southActiveIdx, northActiveIdx;
        if(ii == grid.Nx-1)
          eastActiveIdx = -3;
        else
          eastActiveIdx = mappingGlobalToActive(eastElemIdx);

        if(ii == 0)
          westActiveIdx = -3;
        else
          westActiveIdx = mappingGlobalToActive(westElemIdx);

        if(kk == 0)
          southActiveIdx = -3;
        else
          southActiveIdx = mappingGlobalToActive(southElemIdx);

        if(kk == grid.Ny-1)
          northActiveIdx = -3;
        else
          northActiveIdx = mappingGlobalToActive(northElemIdx);

        // Cell has neumann condition
        if(grid.cellHasPressureBC(flatElemIdx,0) < 0)
        {
          // Condition on the y direction, no cell/masked out cell bellow current one
          if(grid.cellHasPressureBC(flatElemIdx,1) == 1)
          {
            RHS(activeIdx) += 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            A(activeIdx,northActiveIdx) += 2.0*yFactor;
            if(grid.cellHasPressureBC(westElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(westElemIdx,0)-1)*xFactor;
              A(activeIdx,eastActiveIdx) += xFactor;
            }
            else if(grid.cellHasPressureBC(eastElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(eastElemIdx,0)-1)*xFactor;
              A(activeIdx,westActiveIdx) += xFactor;
            }
            else
            {
              A(activeIdx,eastActiveIdx) += xFactor;
              A(activeIdx,westActiveIdx) += xFactor;
            }
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Condition on the x direction, no cell/masked out cell to right of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 2)
          {
            RHS(activeIdx) -= 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            A(activeIdx,westActiveIdx) += 2.0*xFactor;
            if(grid.cellHasPressureBC(northElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(northElemIdx,0)-1)*yFactor;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            else if(grid.cellHasPressureBC(southElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(southElemIdx,0)-1)*yFactor;
              A(activeIdx,northActiveIdx) += yFactor;
            }
            else
            {
              A(activeIdx,northActiveIdx) += yFactor;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            // Need to do a backward difference in x here
            RHS(activeIdx) += (H(flatElemIdx) - H(westElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Condition on the y direction, no cell/masked out cell to top of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 3)
          {
            RHS(activeIdx) -= 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            A(activeIdx,southActiveIdx) += 2.0*yFactor;
            if(grid.cellHasPressureBC(westElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(westElemIdx,0)-1)*xFactor;
              A(activeIdx,eastActiveIdx) += xFactor;
            }
            else if(grid.cellHasPressureBC(eastElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(eastElemIdx,0)-1)*xFactor;
              A(activeIdx,westActiveIdx) += xFactor;
            }
            else
            {
              A(activeIdx,eastActiveIdx) += xFactor;
              A(activeIdx,westActiveIdx) += xFactor;
            }
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            // Need to do a backward difference in y here
            RHS(activeIdx) += (H(flatElemIdx) - H(southElemIdx))/delY;
          }
          // Condition on the x direction, no cell/masked out cell to left of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 4)
          {
            RHS(activeIdx) += 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            A(activeIdx,eastActiveIdx) += 2.0*xFactor;
            if(grid.cellHasPressureBC(northElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(northElemIdx,0)-1)*yFactor;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            else if(grid.cellHasPressureBC(southElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(southElemIdx,0)-1)*yFactor;
              A(activeIdx,northActiveIdx) += yFactor;
            }
            else
            {
              A(activeIdx,northActiveIdx) += yFactor;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 5)
          {
            RHS(activeIdx) -= 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) += 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            A(activeIdx,westActiveIdx) = 2.0*xFactor;
            A(activeIdx,northActiveIdx) = 2.0*yFactor;
            // Need to do backward difference on x
            RHS(activeIdx) += (H(flatElemIdx) - H(westElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and top of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 6)
          {
            RHS(activeIdx) -= 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) -= 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            A(activeIdx,westActiveIdx) = 2.0*xFactor;
            A(activeIdx,southActiveIdx) = 2.0*yFactor;
            // Need to do backward difference on x and on y
            RHS(activeIdx) += (H(flatElemIdx) - H(westElemIdx))/delX;
            RHS(activeIdx) += (H(flatElemIdx) - H(southElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to left and top of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 7)
          {
            RHS(activeIdx) += 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) -= 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            A(activeIdx,eastActiveIdx) = 2.0*xFactor;
            A(activeIdx,southActiveIdx) = 2.0*yFactor;
            // Need to do backward difference on on y
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(flatElemIdx) - H(southElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to left and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 8)
          {
            RHS(activeIdx) += 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) += 2.0*grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            A(activeIdx,eastActiveIdx) = 2.0*xFactor;
            A(activeIdx,northActiveIdx) = 2.0*yFactor;
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
        }
        // No neumann boundary conditions on current cell
        else
        {
          bool westActive = true;
          bool eastActive = true;
          bool northActive = true;
          bool southActive = true;

          if(westActiveIdx < 0)
          {
            RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(westElemIdx,0)-1)*xFactor;
            westActive = false;
          }
          if(eastActiveIdx < 0)
          {
            RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(eastElemIdx,0)-1)*xFactor;
            eastActive = false;
          }
          if(northActiveIdx < 0)
          {
            RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(northElemIdx,0)-1)*yFactor;
            northActive = false;
          }
          if(southActiveIdx < 0)
          {
            RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(southElemIdx,0)-1)*yFactor;
            southActive = false;
          }

          if(eastActive)
            A(activeIdx,eastActiveIdx) += xFactor;
          if(westActive)
            A(activeIdx,westActiveIdx) += xFactor;
          if(northActive)
            A(activeIdx,northActiveIdx) += yFactor;
          if(southActive)
            A(activeIdx,southActiveIdx) += yFactor;
        }
        A(activeIdx,activeIdx) = -2*diagonalFactor;
      }
    }
  }
  /*
  std::cout<<"A MATRIX: "<<std::endl;
  std::cout<<A<<std::endl;
  std::cout<<"B VECTOR: "<<std::endl;
  std::cout<<RHS<<std::endl;
  */
  Eigen::SparseMatrix<double> sparseA = A.sparseView();


}

void ExplicitPressurePoisson::createMappings(CartesianGrid &grid)
{
 // Construct mapping between active center nodes and global center nodes
  mappingGlobalToActive.setZero(grid.numCells);
  std::vector<int> mappingActiveToGlobalVector;
  int dofCtr = 0;
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    bool isMasked = false;
    int faceBCCtr = 0;
    for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
    {
      if(grid.BCTypes(grid.faceMap(faceNum,cellNum),0) != 0)
        faceBCCtr++;
    }
    if (faceBCCtr == grid.numFacesPerElement)
      isMasked = true;

    if(isMasked)
      mappingGlobalToActive(cellNum) = -2;
    else
    {
      // Is a dirichlet BC
      if(grid.cellHasPressureBC(cellNum,0) > 0)
        mappingGlobalToActive(cellNum) = -1;
      // Is a Neumann BC
      else
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
}

void ExplicitPressurePoisson::clearData()
{

;
}

void ExplicitPressurePoisson::calculatePressureGradient(CartesianGrid &grid, Eigen::MatrixXd H)
{
  ;
}