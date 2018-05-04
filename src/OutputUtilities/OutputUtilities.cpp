//
// Created by maytee on 4/27/18.
//

#include "OutputUtilities.h"
void OutputUtilities::writeCartesianCellDataToVTU(CartesianGrid &grid, std::string fileName)
{
  std::ofstream outputFp;
  outputFp.open(fileName);
  int numNodes = grid.numCorners;
  int numCells = grid.numCells;

  outputFp<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  outputFp<<" <UnstructuredGrid>\n";
  outputFp<<"  <Piece NumberOfPoints=\""+std::to_string(numNodes)+"\" NumberOfCells=\""+std::to_string(numCells)+"\">\n";
  outputFp<<"   <Points>\n";
  outputFp<<"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int ii = 0; ii < numNodes;ii++)
           outputFp<<"    "+std::to_string(grid.cornerXY(ii*2))+" "+std::to_string(grid.cornerXY(ii*2+1))+" 0.0\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"   </Points>\n";
  outputFp<<"   <PointData Vectors=\"vectors\">\n";
  outputFp<<"   </PointData>\n";
  outputFp<<"   <Cells>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for(int ii = 0; ii < numCells;ii++)
          outputFp<<"    "+std::to_string(grid.cornerMap(0,ii))+" "+std::to_string(grid.cornerMap(1,ii))+" "+
                                  std::to_string(grid.cornerMap(2,ii))+" "+std::to_string(grid.cornerMap(3,ii))+"\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for(int ii = 0; ii < numCells;ii++)
      outputFp<<"    "+std::to_string(4*(ii+1))+"\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for(int ii = 0; ii < numCells;ii++)
      outputFp<<"    9\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"   </Cells>\n";
  outputFp<<"   <CellData Vectors=\"vectors\">\n";
  Eigen::VectorXd cellVelocities;
  Eigen::VectorXd cellPressureGradients;
  Eigen::VectorXd cellAdvectionFlux;
  Eigen::VectorXd cellDiffusionFlux;
  cellVelocities.setZero(grid.numCells*2);
  cellPressureGradients.setZero(grid.numCells*2);
  cellAdvectionFlux.setZero(grid.numCells*2);
  cellDiffusionFlux.setZero(grid.numCells*2);

  for(int ii = 0; ii < grid.numCells; ii++)
  {
    cellVelocities(ii*2) = grid.uVel(ii);
    cellVelocities(ii*2+1) = grid.vVel(ii);

    cellPressureGradients(ii*2) = grid.pressureGradient(ii,0);
    cellPressureGradients(ii*2+1) = grid.pressureGradient(ii,1);

    cellAdvectionFlux(ii*2) = grid.advectionFluxes(ii,0);
    cellAdvectionFlux(ii*2+1) = grid.advectionFluxes(ii,1);

    cellDiffusionFlux(ii*2) = grid.diffusionFluxes(ii,0);
    cellDiffusionFlux(ii*2+1) = grid.diffusionFluxes(ii,1);

  }
  outputFp<<writeCellVector(cellVelocities, grid.nDim,"cellCenteredVelocity");
  Eigen::VectorXi pressureBCs = grid.cellHasPressureBC.col(1);
  outputFp<<writeCellScalar(pressureBCs,grid.nDim,"pressureBC");
  outputFp<<writeCellScalar(grid.mappingGlobalToActive,grid.nDim,"dofNumber");
  outputFp<<writeCellScalar(grid.pressure,grid.nDim,"pressure");
  outputFp<<writeCellVector(cellPressureGradients,grid.nDim,"pressureGradients");
  outputFp<<writeCellVector(cellAdvectionFlux,grid.nDim,"advectionFlux");
  outputFp<<writeCellVector(cellDiffusionFlux,grid.nDim,"diffusionFlux");
  outputFp<<"   </CellData>\n";
  outputFp<<"  </Piece>\n";
  outputFp<<" </UnstructuredGrid>\n";
  outputFp<<"</VTKFile>\n";

  outputFp.close();
}

// Paraview can't visualize face data so need to work around, and create an unstructured mesh
// of line elements. The CellData function above can then be superimposed ontop of this FaceData
void OutputUtilities::writeCartesianFaceDataToVTU(CartesianGrid &grid, std::string fileName)
{
  std::ofstream outputFp;
  outputFp.open(fileName);
  int numNodes = grid.numCorners;
  int numCells = grid.numCells;
  int numFaces = grid.numFaces;
  outputFp<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  outputFp<<" <UnstructuredGrid>\n";
  outputFp<<"  <Piece NumberOfPoints=\""+std::to_string(numNodes)+"\" NumberOfCells=\""+std::to_string(numFaces)+"\">\n";
  outputFp<<"   <Points>\n";
  outputFp<<"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int ii = 0; ii < numNodes;ii++)
           outputFp<<"    "+std::to_string(grid.cornerXY(ii*2))+" "+std::to_string(grid.cornerXY(ii*2+1))+" 0.0\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"   </Points>\n";
  outputFp<<"   <PointData Vectors=\"vectors\">\n";
  outputFp<<"   </PointData>\n";
  outputFp<<"   <Cells>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  int faceCtr = 0;
  for(int kk = 0; kk < grid.Ny*grid.nDim;kk++)
  {
    for(int ii = 0; ii < grid.Nx; ii++)
    {
      int flatElemIdx = ii + (kk/2)*grid.Nx;
      // Bottom faces
      if(kk%grid.nDim == 0)
      {
        outputFp<<"    "+std::to_string(grid.cornerMap(0,flatElemIdx))+" "+std::to_string(grid.cornerMap(1,flatElemIdx))+"\n";
        faceCtr++;
      }
      // Left faces
      if(kk%grid.nDim != 0)
      {
        outputFp<<"    "+std::to_string(grid.cornerMap(0,flatElemIdx))+" "+std::to_string(grid.cornerMap(3,flatElemIdx))+"\n";
        faceCtr++;
      }
    }
    // Right faces
    int flatElemIdx = grid.Nx-1 + (kk/2)*grid.Nx;
    if(kk%grid.nDim != 0)
    {
      outputFp<<"    "+std::to_string(grid.cornerMap(1,flatElemIdx))+" "+std::to_string(grid.cornerMap(2,flatElemIdx))+"\n";
      faceCtr++;
    }
  }

  for(int ii = 0; ii < grid.Nx; ii++)
  {
    int flatElemIdx = ii + (grid.Ny-1)*grid.Nx;
    // Top faces
    outputFp<<"    "+std::to_string(grid.cornerMap(2,flatElemIdx))+" "+std::to_string(grid.cornerMap(3,flatElemIdx))+"\n";
    faceCtr++;
  }

  outputFp<<"    </DataArray>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for(int ii = 0; ii < numFaces;ii++)
      outputFp<<"    "+std::to_string(2*(ii+1))+"\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for(int ii = 0; ii < numFaces;ii++)
      outputFp<<"    3\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"   </Cells>\n";
  outputFp<<"   <CellData Vectors=\"vectors\">\n";
  Eigen::VectorXi faceBCs;
  faceBCs.setZero(grid.numFaces*2);
  for(int ii = 0; ii < grid.numFaces; ii++)
  {
    faceBCs(ii*2) = grid.BCTypes(ii,0);
    faceBCs(ii*2+1) = grid.BCTypes(ii,1);
  }
  outputFp<<writeCellVector(faceBCs,grid.nDim,"faceVelocityBCs");
  outputFp<<"   </CellData>\n";
  outputFp<<"  </Piece>\n";
  outputFp<<" </UnstructuredGrid>\n";
  outputFp<<"</VTKFile>\n";

  outputFp.close();
}

std::string OutputUtilities::writeCellScalar(Eigen::VectorXd inputVector, int nDim, std::string inputString)
{
  std::string outputStr;
  outputStr = "<DataArray type=\"Float32\" Name=\""+inputString+"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(int ii = 0; ii < inputVector.size(); ii++)
  {
    outputStr += std::to_string(inputVector(ii))+"\n";
  }
  outputStr += "</DataArray>\n";
  return outputStr;
}

std::string OutputUtilities::writeCellScalar(Eigen::VectorXi inputVector, int nDim, std::string inputString)
{
  std::string outputStr;
  outputStr = "<DataArray type=\"Int32\" Name=\""+inputString+"\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(int ii = 0; ii < inputVector.size(); ii++)
  {
    outputStr += std::to_string(inputVector(ii))+"\n";
  }
  outputStr += "</DataArray>\n";
  return outputStr;
}

std::string OutputUtilities::writeCellVector(Eigen::VectorXd inputVector, int nDim, std::string inputString)
{
  std::string outputStr;
  outputStr = "<DataArray type=\"Float32\" Name=\""+inputString+"\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int ii = 0; ii < inputVector.size()/nDim; ii++)
  {
    outputStr += std::to_string(inputVector(ii*nDim))+" "+std::to_string(inputVector(ii*nDim+1))+" 0.0\n";
  }
  outputStr += "</DataArray>\n";
  return outputStr;
}

std::string OutputUtilities::writeCellVector(Eigen::VectorXi inputVector, int nDim, std::string inputString)
{
  std::string outputStr;
  outputStr = "<DataArray type=\"Int32\" Name=\""+inputString+"\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int ii = 0; ii < inputVector.size()/nDim; ii++)
  {
    outputStr += std::to_string(inputVector(ii*nDim))+" "+std::to_string(inputVector(ii*nDim+1))+" 0\n";
  }
  outputStr += "</DataArray>\n";
  return outputStr;
}

void OutputUtilities::writePlottingCartesianDataToVTU(CartesianGrid &grid, std::string fileName)
{
  std::ofstream outputFp;
  outputFp.open(fileName);
  int numNodes = grid.numCorners;
  int numCells = grid.numCells;

  Eigen::VectorXd velocitySums;
  velocitySums.setZero(grid.nodeVelocities.size()/2);
  for(int ii = 0 ; ii < grid.nodeVelocities.size()/2;ii++)
    velocitySums(ii) = grid.nodeVelocities(ii*2) + grid.nodeVelocities(ii*2+1);

  outputFp<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  outputFp<<" <UnstructuredGrid>\n";
  outputFp<<"  <Piece NumberOfPoints=\""+std::to_string(numNodes)+"\" NumberOfCells=\""+std::to_string(numCells)+"\">\n";
  outputFp<<"   <Points>\n";
  outputFp<<"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int ii = 0; ii < numNodes;ii++)
           outputFp<<"    "+std::to_string(grid.cornerXY(ii*2))+" "+std::to_string(grid.cornerXY(ii*2+1))+" "+std::to_string(velocitySums(ii))+"\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"   </Points>\n";
  outputFp<<"   <PointData Vectors=\"vectors\">\n";
  std::string nodeVelocitiesOutputStr = writeCellVector(grid.nodeVelocities, grid.nDim,"cellCenteredVelocity");
  outputFp<<nodeVelocitiesOutputStr;

  std::string nodeVelocitySumsOutputStr = writeCellScalar(velocitySums, grid.nDim,"velocitySum");
  outputFp<<nodeVelocitySumsOutputStr;
  outputFp<<"   </PointData>\n";
  outputFp<<"   <Cells>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for(int ii = 0; ii < numCells;ii++)
          outputFp<<"    "+std::to_string(grid.cornerMap(0,ii))+" "+std::to_string(grid.cornerMap(1,ii))+" "+
                                  std::to_string(grid.cornerMap(2,ii))+" "+std::to_string(grid.cornerMap(3,ii))+"\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for(int ii = 0; ii < numCells;ii++)
      outputFp<<"    "+std::to_string(4*(ii+1))+"\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for(int ii = 0; ii < numCells;ii++)
      outputFp<<"    9\n";
  outputFp<<"    </DataArray>\n";
  outputFp<<"   </Cells>\n";
  outputFp<<"   <CellData Vectors=\"vectors\">\n";
  outputFp<<"   </CellData>\n";
  outputFp<<"  </Piece>\n";
  outputFp<<" </UnstructuredGrid>\n";
  outputFp<<"</VTKFile>\n";

  outputFp.close();
}
