// Finite volume solver for 2D incompressible flows
#include <iostream>
#include "Grids/Grid.h"
#include "Grids/CartesianGrid.h"
#include "OutputUtilities/OutputUtilities.h"
#include "AdvectionOperators/AdvectionUpwind.h"
#include "DiffusionOperators/DiffusionCentral.h"
#include "TimeIntegrators/TimeIntegrator.h"
#include "TimeIntegrators/ForwardEuler.h"
#include <Eigen/Core>

int main() {
  std::cout << "Hello, World!" << std::endl;

  // Generate and initialize grid
  int numX = 20;
  int numY = 20;
  double minX = 0.0;
  double minY = 0.0;
  double Lx = 3.0;
  double Ly = 3.0;
  CartesianGrid grid(numX,numY,Lx,Ly, minX,minY);
  CartesianGrid plottingGrid(numX-1,numY-1,Lx-grid.delX,Ly-grid.delY,grid.delX,grid.delY);
  // Set BCs on grid
  grid.setBCs();

  // Initialize values on grid
  grid.setInitialValues();

  // Interpolate velocity values
  InterpolateUpwind interpObj;

  // Calculate advection terms
  AdvectionUpwind advectionObj;

  // Calculate diffusion terms
  DiffusionCentral diffusionObj;

  double finalTime = 4.0;
  double delT = 0.001;
  std::cout<<"CFL CHECK: "<<delT/(grid.delX*grid.delY)<<std::endl;
  assert(delT/(grid.delX*grid.delY) < 1.0);
  double totalTime = delT;
  OutputUtilities output;
  output.writeCartesianCellDataToVTU(grid,"/home/maytee/Documents/2.29/finiteVolumeSolver/cellData_0.vtu");
  output.writeCartesianFaceDataToVTU(grid,"/home/maytee/Documents/2.29/finiteVolumeSolver/faceData_0.vtu");
  int timeStep = 0;
  double totalNumSteps = 4.0/delT;
  int outputInterval = (int)totalNumSteps/100;
  ForwardEuler integrator;
  while(totalTime < finalTime)
  {
    interpObj.clearData();
    interpObj.interpolateVelocities(grid);

    advectionObj.clearData();
    advectionObj.calculateAdvectionFluxesCartesian(grid,interpObj);

    diffusionObj.clearData();
    //diffusionObj.calculateDiffusionFluxesCartesian(grid);
    Eigen::VectorXd oldVelocities;
    Eigen::VectorXd newVelocities;
    oldVelocities.setZero(grid.numCells*grid.nDim);
    newVelocities.setZero(grid.numCells*grid.nDim);
    for(int ii = 0; ii < grid.numCells; ii++)
    {
      oldVelocities(ii*2) = grid.uVel(ii);
      oldVelocities(ii*2+1) = grid.vVel(ii);
    }
    Eigen::VectorXd newPressure;
    newPressure.setZero(grid.numCells*grid.nDim);
    Eigen::VectorXd fluxTerms;
    fluxTerms.setZero(grid.numCells*grid.nDim);
    for(int ii = 0; ii < grid.numCells; ii++)
    {
      fluxTerms(ii*2) = advectionObj.totalCellAdvectionFlux(ii,0);
      fluxTerms(ii*2+1) = advectionObj.totalCellAdvectionFlux(ii,1);
      //fluxTerms(ii*2) += diffusionObj.totalCellDiffusionFlux(ii,0);
      //fluxTerms(ii*2+1) += diffusionObj.totalCellDiffusionFlux(ii,1);
    }
    integrator.integrate(delT,grid,oldVelocities,newVelocities,newPressure,newPressure,fluxTerms);
    /*
    std::cout<<"Old velocities: "<<std::endl;
    std::cout<<oldVelocities<<std::endl;
    std::cout<<"New velocities: "<<std::endl;
    std::cout<<newVelocities<<std::endl;
    */
    // Update grid variables
    for(int ii = 0; ii < grid.numCells; ii++)
    {
      grid.uVel(ii) = newVelocities(ii*2);
      grid.vVel(ii) = newVelocities(ii*2+1);
    }
    grid.pressure = newPressure;

    totalTime += delT;
    timeStep++;
    if(timeStep % outputInterval == 0)
    {
      std::string fileNumber = std::to_string(timeStep/outputInterval)+".vtu";
      std::cout<<"Saving output at "<<std::to_string(totalTime)<<" seconds to cellData_"+fileNumber<<std::endl;
      output.writeCartesianCellDataToVTU(grid,"/home/maytee/Documents/2.29/finiteVolumeSolver/cellData_"+fileNumber);
      output.writeCartesianFaceDataToVTU(grid,"/home/maytee/Documents/2.29/finiteVolumeSolver/faceData_"+fileNumber);
      std::cout<<"CFL CHECK: "<<delT/(grid.delX*grid.delY)<<std::endl;
    }
  }
  return 0;
}