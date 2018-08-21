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
  n_aug = n_x + P.rows();
  dt = 0;
  // PROPERLY DECLARE ALL STD (std_a + std_yawdd)
  }

//dtor
UKF::~UKF() {
}

void UKF::ProcessMeasurement(const MeasurementPackage &meas_pack) {
  if (!_initialized) {
    std::cout << "initializing!" << std::endl;
    x_ = Eigen::VectorXd(n_x);
    P  = Eigen::MatrixXd(n_x, n_x); // state covariance matrix
    return;
  }

  // else regular update step
  
  // TO BE DONE YET: SET dt CORRECTLY! 
  dt = meas_pack.timestamp_;
  
  ////// Part A: PREDICTION STEP //////
  // Step 1: Generating (augmented) sigma points
  Eigen::MatrixXd Xsig_aug(n_aug, 2 * n_aug + 1);
  AugmentSigmaPoints(Xsig_aug);
  
  // Step 2: Sigma point prediction
  Eigen::MatrixXd Xsig_pred(n_x, 2 * n_aug + 1);
  PredictSigmaPoints(Xsig_pred, Xsig_aug);
  
  // Step 3: Predict state and covariance
  
  
  ////// Part B: PREDICTION STEP //////
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
void UKF::AugmentSigmaPoints(Eigen::MatrixXd &Xsig_aug) {
  
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
  Xsig_aug.colwise() = x_aug;
  Xsig_aug.block(0,       1, n_aug, n_aug) += P_sqrt;
  Xsig_aug.block(0, n_aug+1, n_aug, n_aug) -= P_sqrt;
}

void UKF::PredictSigmaPoints(Eigen::MatrixXd &Xsig_pred, Eigen::MatrixXd &Xsig_aug){
  // double d_t    = delta_t;
  double dt_sq = dt * dt;
  
  for (int i = 0 ; i < (1 + n_aug*2) ; i++) {
      // get augmented state parameters
      double p_x       = Xsig_aug(0, i);
      double p_y       = Xsig_aug(1, i);
      double v         = Xsig_aug(2, i);
      double ups       = Xsig_aug(3, i);
      double ups_d     = Xsig_aug(4, i);
      double nu_a      = Xsig_aug(5, i);
      double nu_ups_dd = Xsig_aug(6, i);
      
      // regular addition
      // use '> 0.001' instead of '==0' to prevent division by zero stability?
      if (ups_d > 0.001) {
          Xsig_pred(0, i) = p_x + v/ups_d*(sin(ups + ups_d*dt) - sin(ups));
          Xsig_pred(1, i) = p_y + v/ups_d*(-cos(ups + ups_d*dt) + cos(ups));
      }
      else {
          Xsig_pred(0, i) = p_x + v*cos(ups)*dt;
          Xsig_pred(1, i) = p_y + v*sin(ups)*dt;
      }
      Xsig_pred(2, i) = v;
      Xsig_pred(3, i) = ups + ups_d * dt;
      Xsig_pred(4, i) = ups_d;
      
      // nu addition
      Xsig_pred(0, i) += 0.5 * dt_sq * cos(ups) * nu_a;
      Xsig_pred(1, i) += 0.5 * dt_sq * sin(ups) * nu_a;
      Xsig_pred(2, i) +=       dt    * nu_a;
      Xsig_pred(3, i) += 0.5 * dt_sq * nu_ups_dd;
      Xsig_pred(4, i) +=       dt    * nu_ups_dd;
  }
}


// print function
void UKF::printA(Eigen::Matrix2d m) {
  a = 3;
  std::cout << "m = " << m << std::endl << std::endl;
}


