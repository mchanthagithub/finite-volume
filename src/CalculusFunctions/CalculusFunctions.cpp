//
// Created by maytee on 4/30/18.
//

#include "CalculusFunctions.h"

// Assume input is in x =[x1,y1,x2,y2,x3,y3,...,xn,yn]
// Use backward difference to approximate differences
Eigen::VectorXd divergenceVector(Eigen::MatrixXd x, int Nx, int Ny, int delX, int delY)
{
  assert(x.rows() == Nx*Ny);
  Eigen::VectorXd result;
  result.setZero(x.rows());

  double firstTerm;
  double secondTerm;
  for(int kk = 0; kk < Ny; kk++)
  {
    for(int ii = 0; ii < Nx; ii++)
    {
      int flatIdx = ii + kk*Nx;
      if(kk == 0)
        secondTerm = (x(ii,kk+1) - x(ii,kk))/delY;
      else
        secondTerm = (x(ii,kk) - x(ii,kk-1))/delY;

      if(ii ==0)
        firstTerm = (x(ii+1,kk) - x(ii,kk))/delX;
      else
        firstTerm = (x(ii,kk) - x(ii-1,kk))/delX;

      result(flatIdx) = firstTerm + secondTerm;
    }
  }
  return result;
}

Eigen::MatrixXd gradientScalar(Eigen::VectorXd x, int Nx, int Ny, int delX, int delY)
{

}
