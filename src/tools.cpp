#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse.fill(0.0);

  // Check the validity of the following inputs:
  // * the estimation vector size should not be zero
  int est_size = estimations.size();
  if(est_size == 0){
    cout << "Error! The estimation vector size should not be zero.\n";
    return rmse;
  }
  // * the estimation vector size should equal ground truth vector size
  if(est_size != ground_truth.size()){
    cout << "Error! The estimation vector size should equal ground truth vector size.\n";
    return rmse;
  }

  // accumulate squared residuals
  VectorXd residual(4);
  for(int i=0;i<est_size;++i){
    residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse /= est_size;

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

void Tools::NormalizeAngle(double &angle){
  if((angle > M_PI) || (angle <= -M_PI)){
    angle -= (2*M_PI) * floor((angle + M_PI) / (2*M_PI));
  }
}
