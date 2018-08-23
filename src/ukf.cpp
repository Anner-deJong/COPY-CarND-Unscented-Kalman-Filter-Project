#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

//ctor
UKF::UKF():
  n_x(5) // const state size
  {
  _initialized = false;
  n_aug = n_x + P_.rows();
  // double lambda = 3 - n_x; // spreading parameter
  prev_timestamp = 0; // last measurment's timestamp for dt calculation
  _initialize_weights(weights_);
  
  /// Standard deviations and covariance matrices ///
  // PROPERLY DECLARE ALL STD (std_a + std_yawdd)
  // process noise nu
  std_a = -999;     // std for longitudinal acceleration
  std_yawdd = -999; // std for yaw acceleration
  
  // Radar
  double rdr_std_r   = 0.3;    // STD - radius in m
  double rdr_std_phi = 0.0175; // STD - angle in rad
  double rdr_std_rd  = 0.1;    // STD - radius change in m/s
  R_radar << rdr_std_r  *  rdr_std_r, 0, 0,
             0, rdr_std_phi*rdr_std_phi, 0,
             0, 0, rdr_std_rd * rdr_std_rd;
  
  // Laser
  double lsr_std_1 = -999;
  double lsr_std_2 = -999;
  R_laser << lsr_std_1*lsr_std_1, 0,
             0, lsr_std_2*lsr_std_2;
  }

//dtor
UKF::~UKF() {
}

void UKF::ProcessMeasurement(const MeasurementPackage &meas_pack) {
  if (!_initialized) {
    std::cout << "initializing!" << std::endl;
    prev_timestamp = meas_pack.timestamp_;
    x_ = VectorXd(n_x);
    P_ = MatrixXd(n_x, n_x); // state covariance matrix 
    
    if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
    }
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
      // x_ << meas_pack.raw_measurements_[0], meas_pack.raw_measurements_[1], 0, 0, 0;
    }
    _initialized = true;
    return;
  }

  // else: regular update step
  // get time elapsed between last two measurements
  dt = (meas_pack.timestamp_ - prev_timestamp) / 1000000.0; //dt expressed in seconds
  prev_timestamp = meas_pack.timestamp_;
  
  ////// Part A: PREDICTION STEP //////
  // Step 1: Generating (augmented) sigma points
  MatrixXd Xsig_aug(n_aug, 2 * n_aug + 1);
  AugmentSigmaPoints(Xsig_aug);
  
  // Step 2: Sigma point prediction
  MatrixXd Xsig_pred(n_x, 2 * n_aug + 1);
  PredictSigmaPoints(Xsig_pred, Xsig_aug);
  
  // Step 3: Predict state mean and covariance
  VectorXd x_pred(n_x);
  MatrixXd P_pred(n_x, n_x);
  PredictMeanAndCovariance(Xsig_pred);
  
  
  ////// Part B: UPDATE STEP //////
  // Step 1: Predict measurements
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar measurement prediction
    int n_z = 3; // radar measurement space size
    VectorXd z_pred = VectorXd(n_z); // predicted measurement
    MatrixXd S_pred = MatrixXd(n_z, n_z); // predicted measurement covariance matrix  
    PredictRadarMeasurement(z_pred, S_pred, Xsig_pred);
  }
  else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Lidar measurement prediction
  }
  
  // Step 2: Update state
  
  printA((Eigen::Matrix2d() << 1, 2, 3, 4).finished());
}





/// FUNCTION IMPLEMENTATIONS ///

// generate sigma points function
Eigen::MatrixXd UKF::GenerateSigmaPoints() {
  // initialize Xsig with each column as the state vector x
  MatrixXd Xsig(n_x, 2*n_x+1);
  Xsig.colwise() = x_;
  
  // multiply the square root matrix of P with the square root constant of lambda - n_x
  MatrixXd A = P_.llt().matrixL();
  A *= sqrt(lambda + n_x);

  // add and subtract the A matrix at the correct slices within Xsig
  Xsig.block(0,     1, n_x, n_x) += A;
  Xsig.block(0, n_x+1, n_x, n_x) -= A;

  return Xsig;
}

// augment sigma points
void UKF::AugmentSigmaPoints(Eigen::MatrixXd &Xsig_aug) {
  
  //create augmented mean state
  VectorXd x_aug(n_aug); 
  x_aug.head(n_x) = x_;
  x_aug.tail(n_aug - n_x).setZero();
  
  //create augmented covariance matrix
  MatrixXd P_aug(n_aug, n_aug);
  P_aug.setZero();
  P_aug.topLeftCorner(n_x, n_x) = P_;
  P_aug(n_aug-2, n_aug-2) = std_a*std_a;
  P_aug(n_aug-1, n_aug-1) = std_yawdd*std_yawdd;
  
  //create square root matrix
  // include the constant term sqrt(lambda + n_aug)
  MatrixXd P_sqrt = P_aug.llt().matrixL();
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
      double phi       = Xsig_aug(3, i);
      double phi_d     = Xsig_aug(4, i);
      double nu_a      = Xsig_aug(5, i);
      double nu_phi_dd = Xsig_aug(6, i);
      
      // regular addition
      // use '> 0.001' instead of '==0' to prevent division by zero stability?
      if (phi_d > 0.001) {
          Xsig_pred(0, i) = p_x + v/phi_d*(sin(phi + phi_d*dt) - sin(phi));
          Xsig_pred(1, i) = p_y + v/phi_d*(-cos(phi + phi_d*dt) + cos(phi));
      }
      else {
          Xsig_pred(0, i) = p_x + v*cos(phi)*dt;
          Xsig_pred(1, i) = p_y + v*sin(phi)*dt;
      }
      Xsig_pred(2, i) = v;
      Xsig_pred(3, i) = phi + phi_d * dt;
      Xsig_pred(4, i) = phi_d;
      
      // nu addition
      Xsig_pred(0, i) += 0.5 * dt_sq * cos(phi) * nu_a;
      Xsig_pred(1, i) += 0.5 * dt_sq * sin(phi) * nu_a;
      Xsig_pred(2, i) +=       dt    * nu_a;
      Xsig_pred(3, i) += 0.5 * dt_sq * nu_phi_dd;
      Xsig_pred(4, i) +=       dt    * nu_phi_dd;
  }
}

// Predict state mean and covariance
void UKF::PredictMeanAndCovariance(MatrixXd &Xsig_pred) {
  
  //predict state mean
  MatrixXd tiled_weights      = weights_.transpose().replicate(n_x, 1);
  MatrixXd weighted_Xsig_pred = Xsig_pred.array() * tiled_weights.array();
  x_ = weighted_Xsig_pred.rowwise().sum();
  _normalize_angle(x_(3));
  
  //predict state covariance matrix
  MatrixXd norm_Xsig_pred = Xsig_pred.colwise() - x_;
  VectorXd tmp_vec = norm_Xsig_pred.row(3);
  // This way of normalizing seems memory inefficient, isnt there a way to pass the row directly to a function??
  _normalize_angle(tmp_vec);
  norm_Xsig_pred.row(3) = tmp_vec;
  MatrixXd weighted_norm_Xsig_pred = norm_Xsig_pred.array() * tiled_weights.array();
  P_ = weighted_norm_Xsig_pred * norm_Xsig_pred.transpose();
  
}

// predict radar measurement
void UKF::PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S_pred, MatrixXd &Xsig_pred) {
  // initialize both predicted measurement and covariance matrix with zeros
  z_pred.setZero();
  S_pred.setZero();
  MatrixXd Zsig(z_pred.size(), 2*n_aug + 1);
  
  for (int i=0 ; i < 2*n_aug + 1 ; i++ ) {
    // store sigma points into seperate dimensions for readability
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v   = Xsig_pred(2, i);
    double phi = Xsig_pred(3, i);
    
    //transform sigma points into measurement space
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);               // rho
    Zsig(1, i) = atan2(p_y, p_x);                       // phi
    Zsig(2, i) = v*(p_x*cos(phi) + p_y*sin(phi))
                       / sqrt(p_x*p_x + p_y*p_y);       // phi_dot
    
    //calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate innovation(?) covariance matrix S
  VectorXd z_norm;
  for (int i=0 ; i < 2*n_aug + 1 ; i++ ) {
    z_norm = Zsig.col(i) - z_pred;
    _normalize_angle(z_norm(1));
    S_pred += weights_(i) * z_norm * z_norm.transpose();
  }
  // Adding the measurement covariance noise matrix to state covariance matrix
  S_pred += R_radar;
}

// print function
void UKF::printA(Eigen::Matrix2d m) {
  a = 3;
  std::cout << "m = " << m << std::endl << std::endl;
}

// initialize weight vector
void UKF::_initialize_weights(VectorXd &weights) {
  weights.setZero();
  weights(0) = lambda / (lambda + n_aug);
  weights.tail(2*n_aug).array() += 1 / (2*(lambda + n_aug));
}

// angle normalization
void UKF::_normalize_angle(double &angle){
  angle = atan2(sin(angle), cos(angle));
}
void UKF::_normalize_angle(VectorXd &angles){
  for (int i=0; i<angles.size(); i++) {
    _normalize_angle(angles[i]);
  }
  
}
