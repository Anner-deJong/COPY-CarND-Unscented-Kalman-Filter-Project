#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

//ctor
UKF::UKF():
  n_x(5) // const state size
  {
  _initialized = false;
  P = Eigen::MatrixXd(n_x, n_x); // state covariance matrix
  
  }

//dtor
UKF::~UKF() {
}

void UKF::ProcessMeasurement(const MeasurementPackage &meas_pack) {
  if (!_initialized) {
    std::cout << "initializing!" << std::endl;
    x_ = Eigen::VectorXd(n_x);
    P  = Eigen::MatrixXd(n_x, n_x);
    return;
  }

  // else a regular update step
  Eigen::MatrixXd m(1,2);
  m << 3, 4;
  printA((Eigen::Matrix2d() << 1, 2, 3, 4).finished());
}

// generate sigma points function
Eigen::MatrixXd UKF::GenerateSigmaPoints() {
  // initialize Xsig with each column as the state vector x
  Eigen::MatrixXd Xsig(n_x, 2*n_x+1);
  Xsig.colwise() = x_;
  
  // multiply the square root matrix of P with the square root constant of lambda - n_x
  Eigen::MatrixXd A = P.llt().matrixL();
  double lambda = 3 - n_x; // spreading parameter (seems like this could be reduced to single line)
  A *= sqrt(lambda + n_x);

  // add and subtract the A matrix at the correct slices within Xsig
  Xsig.block(0,     1, n_x, n_x) += A;
  Xsig.block(0, n_x+1, n_x, n_x) -= A;

  return Xsig;
}

// print function
void UKF::printA(Eigen::Matrix2d m) {
  a = 3;
  std::cout << "m = " << m << std::endl << std::endl;
}


