#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

//ctor
UKF::UKF():
  n_x(5), // const state size
  std_a(-999), // std for longitudinal acceleration
  std_yawdd(-999) // std for yaw acceleration
  {
  _initialized = false;
  P = Eigen::MatrixXd(n_x, n_x); // state covariance matrix
  n_aug = n_x + P.rows();
  Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(n_aug, 2 * n_aug + 1);
  // PROPERLY DECLARE ALL STD (std_a + std_yawdd)
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

// augment sigma points
void UKF::AugmentSigmaPoints(Eigen::MatrixXd *Xsig_aug) {
  
  //create augmented mean state
  Eigen::VectorXd x_aug = Eigen::VectorXd(7);
  x_aug.head(n_x) = x;
  x_aug.tail(n_aug - n_x).setZero();
  
  //create augmented covariance matrix
  Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug, n_aug);
  P_aug.setZero();
  P_aug.topLeftCorner(n_x, n_x) = P;
  P_aug(n_aug-2, n_aug-2) = std_a*std_a;
  P_aug(n_aug-1, n_aug-1) = std_yawdd*std_yawdd;
  
  //create square root matrix
  // include the constant term sqrt(lambda + n_aug)
  Eigen::MatrixXd P_sqrt = P_aug.llt().matrixL();
  double lambda = 3 - n_x; // again, spreading parameter (seems like this could be reduced to single line)
  P_sqrt *= sqrt(lambda + n_aug);
  
  //create augmented sigma points
  Xsig_aug->colwise() = x_aug;
  Xsig_aug->block(0,       1, n_aug, n_aug) += P_sqrt;
  Xsig_aug->block(0, n_aug+1, n_aug, n_aug) -= P_sqrt;  
  
}


// print function
void UKF::printA(Eigen::Matrix2d m) {
  a = 3;
  std::cout << "m = " << m << std::endl << std::endl;
}


