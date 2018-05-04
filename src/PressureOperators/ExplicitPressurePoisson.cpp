//
// Created by maytee on 5/2/18.
//

#include "ExplicitPressurePoisson.h"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "Eigen/SparseCholesky"
#include "Eigen/IterativeLinearSolvers"
#include "iostream"
void ExplicitPressurePoisson::calculatePressure(CartesianGrid &grid, Eigen::MatrixXd H)
{
  if(!haveCreatedMappings)
    createMappings(grid);
  if(!haveConstructedAMatrix)
    constructAMatrix(grid);

  constructRHSVector(grid,H);

  //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

  solver.compute(sparseA);

  //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
  //solver.compute(sparseA);
  activePressure = solver.solve(RHS);
  //std::cout<<"active pressure: "<<std::endl;
  //std::cout<<activePressure<<std::endl;



  int dofCtr = 0;
  for(int cellNum = 0; cellNum < grid.numCells; cellNum++)
  {
    if(mappingGlobalToActive(cellNum) >= 0)
    {
      pressure(cellNum) = activePressure(dofCtr);
      dofCtr++;
      // Re-adjust the pressure to get rid of the 0 pressure term
      //if(dofCtr == activePressure.size())
        //pressure(cellNum) = activePressure(dofCtr-2);
    }
    else if(mappingGlobalToActive(cellNum) == -2)
    {
      pressure(cellNum) = 0;
    }
    else
    {
      int mapping = mappingGlobalToActive(cellNum);
      int blah = grid.cellHasPressureBC(cellNum,0)-1;
      pressure(cellNum) = grid.pressureDirichletBCValue(grid.cellHasPressureBC(cellNum,0)-1);
      //std::cout<<pressure(cellNum)<<std::endl;
    }
  }
  //std::cout<<"total pressure: "<<std::endl;
  //std::cout<<pressure<<std::endl;
}

void ExplicitPressurePoisson::constructAMatrix(CartesianGrid& grid)
{
  if(!haveCreatedMappings)
    createMappings(grid);
  // Create matrix to hold central difference operator
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(totalDOF,totalDOF);
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
            A(activeIdx,northActiveIdx) += yFactor;
            if(grid.cellHasPressureBC(westElemIdx,0) > 0)
            {
              A(activeIdx,eastActiveIdx) += xFactor;
            }
            else if(grid.cellHasPressureBC(eastElemIdx,0) > 0)
            {
              A(activeIdx,westActiveIdx) += xFactor;
            }
            else
            {
              A(activeIdx,eastActiveIdx) += xFactor;
              A(activeIdx,westActiveIdx) += xFactor;
            }
            A(activeIdx,activeIdx) -= (2.0*xFactor + yFactor);
          }
          // Condition on the x direction, no cell/masked out cell to right of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 2)
          {
            A(activeIdx,westActiveIdx) += xFactor;
            if(grid.cellHasPressureBC(northElemIdx,0) > 0)
            {
              std::cout<<"HERE1"<<std::endl;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            else if(grid.cellHasPressureBC(southElemIdx,0) > 0)
            {
              std::cout<<"HERE2"<<std::endl;
              A(activeIdx,northActiveIdx) += yFactor;
            }
            else
            {
              A(activeIdx,northActiveIdx) += yFactor;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            A(activeIdx,activeIdx) -= (xFactor + 2.0*yFactor);
          }
          // Condition on the y direction, no cell/masked out cell to top of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 3)
          {
            A(activeIdx,southActiveIdx) += yFactor;
            if(grid.cellHasPressureBC(westElemIdx,0) > 0)
            {
              A(activeIdx,eastActiveIdx) += xFactor;
            }
            else if(grid.cellHasPressureBC(eastElemIdx,0) > 0)
            {
              A(activeIdx,westActiveIdx) += xFactor;
            }
            else
            {
              A(activeIdx,eastActiveIdx) += xFactor;
              A(activeIdx,westActiveIdx) += xFactor;
            }
            A(activeIdx,activeIdx) -= (2.0*xFactor + yFactor);
          }
          // Condition on the x direction, no cell/masked out cell to left of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 4)
          {
            A(activeIdx,eastActiveIdx) += xFactor;
            if(grid.cellHasPressureBC(northElemIdx,0) > 0)
            {
              A(activeIdx,southActiveIdx) += yFactor;
            }
            else if(grid.cellHasPressureBC(southElemIdx,0) > 0)
            {
              A(activeIdx,northActiveIdx) += yFactor;
            }
            else
            {
              A(activeIdx,northActiveIdx) += yFactor;
              A(activeIdx,southActiveIdx) += yFactor;
            }
            A(activeIdx,activeIdx) -= (xFactor + 2.0*yFactor);
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 5)
          {
            A(activeIdx,westActiveIdx) = xFactor;
            A(activeIdx,northActiveIdx) = yFactor;
            A(activeIdx,activeIdx) -= (xFactor + yFactor);
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and top of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 6)
          {
            A(activeIdx,westActiveIdx) = xFactor;
            A(activeIdx,southActiveIdx) = yFactor;
            A(activeIdx,activeIdx) -= (xFactor + yFactor);
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to left and top of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 7)
          {
            A(activeIdx,eastActiveIdx) = xFactor;
            A(activeIdx,southActiveIdx) = yFactor;
            A(activeIdx,activeIdx) -= (xFactor + yFactor);
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to left and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 8)
          {
            A(activeIdx,eastActiveIdx) = xFactor;
            A(activeIdx,northActiveIdx) = yFactor;
            A(activeIdx,activeIdx) -= (xFactor + yFactor);
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
            westActive = false;
          }
          if(eastActiveIdx < 0)
          {
            eastActive = false;
          }
          if(northActiveIdx < 0)
          {
            northActive = false;
          }
          if(southActiveIdx < 0)
          {
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

          A(activeIdx,activeIdx) -= 2.0*xFactor + 2.0*yFactor;
        }
      }
    }
  }
  //std::cout<<"A: "<<std::endl;
  //std::cout<<A<<std::endl;

  // If all neumann pressure conditions then need to deal with singular matrix
  if(grid.isAllNeumannPressureBCs)
  {
    Eigen::VectorXd temp = Eigen::VectorXd::Ones(totalDOF);
    A.conservativeResize(A.rows()+1,A.cols());
    //temp(totalDOF-1) = 1.0;
    //A.row(totalDOF-1) = temp;
    A.row(A.rows()-1) = temp;
  }

  //std::cout<<"A after singular treatment: "<<std::endl;
  //std::cout<<A<<std::endl;

  sparseA = A.sparseView();

  haveConstructedAMatrix = true;
}

void ExplicitPressurePoisson::constructRHSVector(CartesianGrid &grid, Eigen::MatrixXd H)
{
  RHS.setZero(totalDOF);
  // Create matrix to hold central difference operator
  double delX = grid.delX;
  double delY = grid.delY;
  double xFactor = 1.0/(grid.delX*grid.delX);
  double yFactor = 1.0/(grid.delY*grid.delY);

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
            RHS(activeIdx) += grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            if(grid.cellHasPressureBC(westElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(westElemIdx,0)-1)*xFactor;
            }
            else if(grid.cellHasPressureBC(eastElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(eastElemIdx,0)-1)*xFactor;
            }
            else
            {
            }
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Condition on the x direction, no cell/masked out cell to right of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 2)
          {
            RHS(activeIdx) -= grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            if(grid.cellHasPressureBC(northElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(northElemIdx,0)-1)*yFactor;
            }
            else if(grid.cellHasPressureBC(southElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(southElemIdx,0)-1)*yFactor;
            }
            else
            {
            }
            // Need to do a backward difference in x here
            RHS(activeIdx) += (H(flatElemIdx) - H(westElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Condition on the y direction, no cell/masked out cell to top of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 3)
          {
            RHS(activeIdx) -= grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            if(grid.cellHasPressureBC(westElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(westElemIdx,0)-1)*xFactor;
            }
            else if(grid.cellHasPressureBC(eastElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(eastElemIdx,0)-1)*xFactor;
            }
            else
            {
            }
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            // Need to do a backward difference in y here
            RHS(activeIdx) += (H(flatElemIdx) - H(southElemIdx))/delY;
          }
          // Condition on the x direction, no cell/masked out cell to left of current one
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 4)
          {
            RHS(activeIdx) += grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            if(grid.cellHasPressureBC(northElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(northElemIdx,0)-1)*yFactor;
            }
            else if(grid.cellHasPressureBC(southElemIdx,0) > 0)
            {
              RHS(activeIdx) -= grid.pressureDirichletBCValue(grid.cellHasPressureBC(southElemIdx,0)-1)*yFactor;
            }
            else
            {
            }
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 5)
          {
            RHS(activeIdx) -= grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) += grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            // Need to do backward difference on x
            RHS(activeIdx) += (H(flatElemIdx) - H(westElemIdx))/delX;
            RHS(activeIdx) += (H(northElemIdx) - H(flatElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and top of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 6)
          {
            RHS(activeIdx) -= grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) -= grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            // Need to do backward difference on x and on y
            RHS(activeIdx) += (H(flatElemIdx) - H(westElemIdx))/delX;
            RHS(activeIdx) += (H(flatElemIdx) - H(southElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to left and top of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 7)
          {
            RHS(activeIdx) += grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) -= grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
            // Need to do backward difference on on y
            RHS(activeIdx) += (H(eastElemIdx) - H(flatElemIdx))/delX;
            RHS(activeIdx) += (H(flatElemIdx) - H(southElemIdx))/delY;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to left and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 8)
          {
            RHS(activeIdx) += grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delX;
            RHS(activeIdx) += grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)/delY;
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

        }
      }
    }
  }
  /*
  std::cout<<"A MATRIX: "<<std::endl;
  std::cout<<A<<std::endl;
  std::cout<<"B VECTOR: "<<std::endl;
  std::cout<<RHS<<std::endl;
  */
  //std::cout<<"RHS: "<<std::endl;
  //std::cout<<RHS<<std::endl;
  // If all neumann pressure conditions then need to treat for singularity
  if(grid.isAllNeumannPressureBCs)
  {
    //RHS(totalDOF-1) = 0.0;
    RHS.conservativeResize(RHS.size()+1);
    RHS(RHS.size()-1) = 0.0;
  }

  //std::cout<<"RHS after singular treatment: "<<std::endl;
  //`std::cout<<RHS<<std::endl;

}

void ExplicitPressurePoisson::calculatePressureGradient(CartesianGrid &grid, Eigen::MatrixXd H)
{
  calculatePressure(grid,H);
  double delX = grid.delX;
  double delY = grid.delY;
  pressureGradient.setZero(grid.numCells,grid.nDim);
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
          if(grid.cellHasPressureBC(flatElemIdx,1) == 1 || grid.cellHasPressureBC(flatElemIdx,1) == 3)
          {
            //std::cout<<"cellnum: "<<flatElemIdx<<std::endl;
            //std::cout<<grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1)<<std::endl;
            // 1: Condition on the y direction, no cell/masked out cell bellow current one
            // 3: Condition on the y direction, no cell/masked out cell to top of current one
            //pressureGradient(flatElemIdx,0) = (pressure(flatElemIdx) - pressure(westElemIdx))/delX;
            pressureGradient(flatElemIdx,0) = (pressure(eastElemIdx) - pressure(westElemIdx))/(2.0*delX);
            pressureGradient(flatElemIdx,1) = grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1);
            //std::cout<<"13cellNum: "<<flatElemIdx<<" "<<pressureGradient(flatElemIdx,1)<<std::endl;
          }
          else if(grid.cellHasPressureBC(flatElemIdx,1) == 2 || grid.cellHasPressureBC(flatElemIdx,1) == 4)
          {
            // 2: Condition on the x direction, no cell/masked out cell to right of current one
            // 4: Condition on the x direction, no cell/masked out cell to left of current one
            pressureGradient(flatElemIdx,0) = grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1);
            //pressureGradient(flatElemIdx,1) = (pressure(flatElemIdx) - pressure(southElemIdx))/delY;
            pressureGradient(flatElemIdx,1) = (pressure(northElemIdx) - pressure(southElemIdx))/(2.0*delY);
            //std::cout<<"24cellNum: "<<flatElemIdx<<" "<<pressureGradient(flatElemIdx,1)<<std::endl;
          }
          // Neumann condition on both x and y direction, no cell/masked out cell to right and bottom of current
          else if(grid.cellHasPressureBC(flatElemIdx,1) >= 5 || grid.cellHasPressureBC(flatElemIdx,1) <= 8)
          {
            pressureGradient(flatElemIdx,0) = grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1);
            pressureGradient(flatElemIdx,0) = grid.pressureNeumannBCValue(-1*grid.cellHasPressureBC(flatElemIdx,0)-1);
            //std::cout<<"58cellNum: "<<flatElemIdx<<" "<<pressureGradient(flatElemIdx,1)<<std::endl;
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
            westActive = false;
          if(eastActiveIdx < 0)
            eastActive = false;
          if(northActiveIdx < 0)
            northActive = false;
          if(southActiveIdx < 0)
            southActive = false;


          // Need to do a forward difference in x
          if(!westActive)
            pressureGradient(flatElemIdx,0) = (pressure(eastElemIdx) - pressure(flatElemIdx))/delX;
          else
          {
            //pressureGradient(flatElemIdx,0) = (pressure(flatElemIdx) - pressure(westElemIdx))/delX;
            pressureGradient(flatElemIdx,0) = (pressure(eastElemIdx) - pressure(westElemIdx))/(2.0*delX);
          }

          // Need to do a forward difference in y
          if(!southActive)
            pressureGradient(flatElemIdx,1) = (pressure(northElemIdx) - pressure(flatElemIdx))/delY;
          else
          {
            //pressureGradient(flatElemIdx,1) = (pressure(flatElemIdx) - pressure(southElemIdx))/delY;
            pressureGradient(flatElemIdx,1) = (pressure(northElemIdx) - pressure(southElemIdx))/(2.0*delY);
          }
          //std::cout<<"NoneumcellNum: "<<flatElemIdx<<" "<<pressureGradient(flatElemIdx,1)<<std::endl;


        }
      }
    }
  }

  //std::cout<<"pressure gradient:"<<std::endl;
  //std::cout<<pressureGradient<<std::endl;

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

  activePressure.setZero(totalDOF);
  pressure.setZero(grid.numCells);
  haveCreatedMappings = true;
}

void ExplicitPressurePoisson::clearData()
{
  activePressure.setZero();
  pressureGradient.setZero();
  pressure.setZero();
  //totalDOF = 0;
  //mappingGlobalToActive.setZero();
  //mappingActiveToGlobal.setZero();
;
}

