//
// Created by maytee on 4/28/18.
//

#include "ForwardEuler.h"

void ForwardEuler::integrate(double delT, Grid &grid, Eigen::VectorXd& oldVelocity, Eigen::VectorXd& newVelocity,
                             Eigen::VectorXd& oldPressureGradient, Eigen::VectorXd& newPressureGradient, Eigen::VectorXd& fluxTerms)
{
  newVelocity = oldVelocity - delT*(fluxTerms+oldPressureGradient);
}
