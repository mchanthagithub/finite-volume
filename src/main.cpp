// Finite volume solver for 2D incompressible flows
#include <iostream>
#include "Grids/Grid.h"
#include "Grids/CartesianGrid.h"
#include "OutputUtilities/OutputUtilities.h"
#include "AdvectionOperators/AdvectionUpwind.h"
#include "DiffusionOperators/DiffusionCentral.h"
#include "PressureOperators/ExplicitPressure.h"
#include "PressureOperators/ExplicitPressurePoisson.h"
#include "TimeIntegrators/TimeIntegrator.h"
#include "TimeIntegrators/ForwardEuler.h"
#include <Eigen/Core>

int main() {
  std::cout << "Hello, World!" << std::endl;
  //double rho = 1.225; // density of air
  //double mu = 18.0*1e-6; // dynamic viscosity of air
  //double mu = 1.0*1e-3; // dynamic viscosity of air
  double rho = 1.0;
  double mu = 1.0;

  double Nu = mu/rho; // kinematic viscosity
  // Generate and initialize grid
  int numX = 40;
  int numY = 40;
  double minX = 0.0;
  double minY = 0.0;
  double Lx = 2.0;
  double Ly = 2.0;
  CartesianGrid grid(numX,numY,Lx,Ly, minX,minY);
  CartesianGrid plottingGrid(numX-1,numY-1,Lx-grid.delX,Ly-grid.delY,grid.delX/2.0,grid.delY/2.0);
  plottingGrid.nodeVelocities.setZero(grid.numCells*grid.nDim);

  // Interpolate velocity values
  InterpolateUpwind interpObj;

  // Calculate advection terms
  AdvectionUpwind advectionObj;

  // Calculate diffusion terms
  DiffusionCentral diffusionObj;

  // Pressure
  ExplicitPressurePoisson poissonPressObj;

  bool hasAdvection = true;
  bool hasDiffusion = true;
  bool hasPressure = true;

  if(grid.isBurgers)
  {
    hasPressure = false;
    hasDiffusion = false;
  }
  if(grid.isDiffusion)
  {
    hasPressure = false;
    hasAdvection = false;
  }

  //hasPressure = false;
  //hasDiffusion = false;
  //hasAdvection = false;
  std::string writeDir = "/home/maytee/Documents/2.29/outputs/";
  double finalTime = 6.0;
  double delT = 0.0002;
  std::cout<<"CFL CHECK: "<<delT/(grid.delX*grid.delY)<<std::endl;
  assert(delT/(grid.delX*grid.delY) < 1.0);
  double totalTime = delT;
  for(int ii = 0; ii < grid.numCells; ii++)
  {
    plottingGrid.nodeVelocities(ii*2) = grid.uVel(ii);
    plottingGrid.nodeVelocities(ii*2+1) = grid.vVel(ii);
  }
  OutputUtilities output;
  output.writeCartesianCellDataToVTU(grid,writeDir+"cellData_0.vtu");
  output.writeCartesianFaceDataToVTU(grid,writeDir+"faceData_0.vtu");
  output.writePlottingCartesianDataToVTU(plottingGrid,writeDir+"plottingData_0.vtu");
  int timeStep = 0;
  double numStepsPerSec = 1.0/delT;
  int numOutputsPerSec = 50;
  int outputInterval = (int)(numStepsPerSec/(numOutputsPerSec));
  std::cout<<"output interval: "<<outputInterval<<std::endl;
  //outputInterval = 1;
  ForwardEuler integrator;
  Eigen::VectorXd fluxTermsVector;
  Eigen::VectorXd newPressureGradient;
  Eigen::VectorXd oldVelocities;
  Eigen::VectorXd newVelocities;
  Eigen::VectorXd faceVels;
  while(totalTime < finalTime)
  {
    std::cout<<"Time is: "<<totalTime<<std::endl;
    fluxTermsVector.setZero(grid.numCells*grid.nDim);
    newPressureGradient.setZero(grid.numCells*grid.nDim);
    oldVelocities.setZero(grid.numCells*grid.nDim);
    newVelocities.setZero(grid.numCells*grid.nDim);
    faceVels.setZero(grid.numFaces*grid.nDim);

    interpObj.clearData();
    interpObj.interpolateVelocities(grid);

    advectionObj.clearData();
    advectionObj.calculateAdvectionFluxesCartesian(grid,interpObj);

    diffusionObj.clearData();
    diffusionObj.calculateDiffusionFluxesCartesian(grid,interpObj);

    Eigen::MatrixXd totalFluxesMatrix;

    if(hasAdvection && hasDiffusion)
      totalFluxesMatrix = mu*diffusionObj.totalCellDiffusionFlux - rho*advectionObj.totalCellAdvectionFlux;
    else if(hasAdvection)
      totalFluxesMatrix = - rho*advectionObj.totalCellAdvectionFlux;
    else if(hasDiffusion)
      totalFluxesMatrix = mu*diffusionObj.totalCellDiffusionFlux;

    poissonPressObj.clearData(grid);
    if(hasPressure)
      poissonPressObj.calculatePressureGradient(grid,totalFluxesMatrix);

    for(int ii = 0; ii < grid.numCells; ii++)
    {
      oldVelocities(ii*2) = grid.uVel(ii);
      oldVelocities(ii*2+1) = grid.vVel(ii);
    }
    double tol = 1e-6;
    for(int ii = 0; ii < grid.numCells; ii++)
    {
      if(hasAdvection)
      {
        fluxTermsVector(ii*2) -= rho*advectionObj.totalCellAdvectionFlux(ii,0);
        fluxTermsVector(ii*2+1) -= rho*advectionObj.totalCellAdvectionFlux(ii,1);
      }
      if(hasDiffusion)
      {
        fluxTermsVector(ii*2) += mu*diffusionObj.totalCellDiffusionFlux(ii,0);
        fluxTermsVector(ii*2+1) += mu*diffusionObj.totalCellDiffusionFlux(ii,1);
      }
      if(hasPressure)
      {
        newPressureGradient(ii*2) += poissonPressObj.pressureGradient(ii,0);
        newPressureGradient(ii*2+1) += poissonPressObj.pressureGradient(ii,1);
      }

      for(int faceNum = 0; faceNum < grid.numFacesPerElement; faceNum++)
      {
        int globalFaceIdx = grid.faceMap(faceNum,ii);
        /*
        if(std::abs(faceVels(globalFaceIdx*2)) > 0.0)
        {
          if(std::abs(faceVels(globalFaceIdx*2) - interpObj.uVelInterp(ii,faceNum)) > tol)
          {
            std::cout<<"cell: "<<ii<<" face: "<<faceNum<<" ERRORX==========================="<<std::endl;
            std::cout<<faceVels(globalFaceIdx*2)<<" "<<interpObj.uVelInterp(ii,faceNum)<<std::endl;
          }
        }
        if(std::abs(faceVels(globalFaceIdx*2+1)) > 0.0)
        {
          if(std::abs(faceVels(globalFaceIdx*2+1) - interpObj.vVelInterp(ii,faceNum)) > tol)
          {
            std::cout<<"cell: "<<ii<<" face: "<<faceNum<<" ERRORY==========================="<<std::endl;
            std::cout<<faceVels(globalFaceIdx*2+1)<<" "<<interpObj.vVelInterp(ii,faceNum)<<std::endl;
          }
        }
        */
        faceVels(globalFaceIdx*2) = interpObj.uVelInterp(ii,faceNum);
        faceVels(globalFaceIdx*2+1) = interpObj.vVelInterp(ii,faceNum);
        //std::cout<<"face: "<<faceNum<<" "<<interpObj.uVelInterp(ii,faceNum)<<" "<<interpObj.vVelInterp(ii,faceNum)<<std::endl;
      }
    }
    grid.faceVels = faceVels;
    grid.oldUVel = grid.uVel;
    grid.oldVVel = grid.vVel;

    integrator.integrate(delT,grid,oldVelocities,newVelocities,newPressureGradient,fluxTermsVector,rho);
    grid.pressure = poissonPressObj.pressure;
    grid.pressureGradient = poissonPressObj.pressureGradient;
    //std::cout<<"grid press grad"<<std::endl;
    //std::cout<<grid.pressureGradient<<std::endl;
    // Print out the max and min velocity components
    std::cout<<"max-v-comp: "<<newVelocities.maxCoeff()<<" min-v-comp: "<<newVelocities.minCoeff()<<std::endl;

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

    totalTime += delT;
    timeStep++;
    if(timeStep % outputInterval == 0)
    {
      for(int ii = 0; ii < grid.numCells; ii++)
      {
        plottingGrid.nodeVelocities(ii*2) = grid.uVel(ii);
        plottingGrid.nodeVelocities(ii*2+1) = grid.vVel(ii);
      }
      grid.advectionFluxes = advectionObj.totalCellAdvectionFlux;
      grid.diffusionFluxes = diffusionObj.totalCellDiffusionFlux;
      std::string fileNumber = std::to_string(timeStep/outputInterval)+".vtu";
      std::cout<<"Saving output at "<<std::to_string(totalTime)<<" seconds to cellData_"+fileNumber<<std::endl;
      output.writeCartesianCellDataToVTU(grid,writeDir+"cellData_"+fileNumber);
      output.writeCartesianFaceDataToVTU(grid,writeDir+"faceData_"+fileNumber);
      output.writePlottingCartesianDataToVTU(plottingGrid,writeDir+"plottingData_"+fileNumber);
      std::cout<<"CFL CHECK: "<<delT/(grid.delX*grid.delY)<<std::endl;
    }
  }
  return 0;
}