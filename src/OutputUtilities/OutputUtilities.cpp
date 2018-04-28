//
// Created by maytee on 4/27/18.
//

#include "OutputUtilities.h"
void OutputUtilities::writeCartesianCellDataToVTU(Grid &grid, std::string fileName)
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
  cellVelocities.setZero(grid.numCells*2);
  for(int ii = 0; ii < grid.numCells; ii++)
  {
    cellVelocities(ii*2) = grid.uVel(ii);
    cellVelocities(ii*2+1) = grid.vVel(ii);
  }
  std::string cellVelocitiesOutputStr = writeCellVector(cellVelocities, grid.nDim,"cellCenteredVelocity");
  outputFp<<cellVelocitiesOutputStr;
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
  outputFp<<"  </Piece>\n";
  outputFp<<" </UnstructuredGrid>\n";
  outputFp<<"</VTKFile>\n";

  outputFp.close();
}

std::string OutputUtilities::writeCellScalar(Eigen::VectorXd inputVector)
{

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
