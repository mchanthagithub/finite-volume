//
// Created by maytee on 4/28/18.
//

#include "ForwardEuler.h"
#include "iostream"
void ForwardEuler::integrate(double delT, Grid &grid, Eigen::VectorXd& oldVelocity, Eigen::VectorXd& newVelocity,
                             Eigen::VectorXd& pressureGradient, Eigen::VectorXd& fluxTerms)
{
  //std::cout<<"fluxTerms: "<<std::endl;
  //std::cout<<fluxTerms<<std::endl;
  newVelocity = oldVelocity + delT*(fluxTerms-pressureGradient);
}
