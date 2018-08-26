#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"



 /**************
 * constructor *
 **************/

UKF::UKF():
  n_x(5),     // state size
  n_aug(7),   // augmented state size
  n_z_rdr(3), // radar measurement space size
  n_z_lsr(2)  // laser measurement space size
  
  {
  
  // miscellaneous variables
  _initialized = false;
  lambda = 3 - n_x; // spreading parameter
  prev_timestamp = 0; // last measurment's timestamp for dt calculation
  _initialize_weights(weights_);
  
  // Sigma calculations related variables size initialization
  Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  
  /// Standard deviations and covariance matrices ///
  
  // process noise nu. TWEAK THESE TWO PARAMS PROPERLY
  std_a = M_PI/4;  // STD - longitudinal acceleration in m/s^2
  std_yawdd = 0.5; // STD - yaw acceleration in rad/s^2
  
  // Process state (optional, tweak these here or in the initialization step)
  double std_px    = 1;
  double std_py    = 1;
  double std_v     = 1;
  double std_psi   = 1;
  double std_psi_d = 1;
  P_ = MatrixXd(n_x, n_x);
  P_ << std_px   *   std_px, 0, 0, 0, 0,
        0, std_py   *   std_py, 0, 0, 0,
        0, 0, std_v    *    std_v, 0, 0,
        0, 0, 0, std_psi  *  std_psi, 0,
        0, 0, 0, 0, std_psi_d*std_psi_d;
  
  // Radar
  double rdr_std_r   = 0.3;  // STD - radius in m
  double rdr_std_phi = 0.03; // STD - angle in rad
  double rdr_std_r_d = 0.3;  // STD - radius change in m/s
  R_radar = MatrixXd(n_z_rdr, n_z_rdr);
  R_radar << rdr_std_r  *  rdr_std_r, 0, 0,
             0, rdr_std_phi*rdr_std_phi, 0,
             0, 0, rdr_std_r_d*rdr_std_r_d;
  
  // Laser
  double lsr_std_px = 0.15; // STD - position x in m
  double lsr_std_py = 0.15; // STD - position y in m
  R_laser = MatrixXd(n_z_lsr, n_z_lsr);
  R_laser << lsr_std_px*lsr_std_px, 0,
             0, lsr_std_py*lsr_std_py;  
  }

//dtor
UKF::~UKF() {
}


 /***********************************************
 * Main interaction API for usage of this class *
 ***********************************************/

void UKF::ProcessMeasurement(const MeasurementPackage &meas_pack) {
  if (!_initialized) {
    prev_timestamp = meas_pack.timestamp_;
    x_ = VectorXd(n_x);
    
    if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_pack.raw_measurements_[0];
      double phi = meas_pack.raw_measurements_[1];
      x_ << rho*cos(phi), rho*sin(phi), 0, 0, 0;
    }
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_pack.raw_measurements_[0], meas_pack.raw_measurements_[1], 0, 0, 0;
    }
    _initialized = true;
    return;
  }
  
  /// else: regular update step ///
  
  // get time elapsed between last two measurements
  dt = (meas_pack.timestamp_ - prev_timestamp) / 1000000.0; //dt expressed in seconds
  prev_timestamp = meas_pack.timestamp_;
  
  ////// Part A: PREDICTION STEP //////
  // Step 1: Generate (augmented) sigma points
  AugmentSigmaPoints();
  
  // Step 2: Predict sigma points
  PredictSigmaPoints();
  
  // Step 3: Predict state mean and covariance
  PredictMeanAndCovariance();
  
  
  ////// Part B: UPDATE STEP //////
  // Step 1: Predict (radar/laser) measurements
   
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar measurement prediction
    n_z = n_z_rdr;
    PredictRadarMeasurement();
    
    // filling z with the correct raw measurement data
    z = VectorXd(n_z_rdr);
    z << meas_pack.raw_measurements_[0], meas_pack.raw_measurements_[1], meas_pack.raw_measurements_[2];
    
  }
  else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Lidar measurement prediction
    // z_pred = VectorXd(n_z_lsr);
    // S_pred = MatrixXd(n_z_lsr, n_z_lsr);
    // Zsig   = MatrixXd(n_z_lsr, 2*n_aug + 1);
    n_z = n_z_lsr;
    PredictLaserMeasurement();

    // filling z with the correct raw measurement data
    z = VectorXd(n_z_lsr);
    z << meas_pack.raw_measurements_[0], meas_pack.raw_measurements_[1];
  }
  
  // Step 2: Update state
  UpdateState(Xsig_pred, z_pred, z, Zsig, S_pred);
}




 /***************************
 * FUNCTION IMPLEMENTATIONS *
 ***************************/

// Part A Step 1: Generate (augmented) sigma points
void UKF::AugmentSigmaPoints() {
  
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

// Part A Step 2: Predict sigma points
void UKF::PredictSigmaPoints() {
  // double d_t    = delta_t;
  double dt_sq = dt * dt;
  
  for (int i = 0 ; i < (1 + n_aug*2) ; i++) {
      // get augmented state parameters
      double p_x       = Xsig_aug(0, i);
      double p_y       = Xsig_aug(1, i);
      double v         = Xsig_aug(2, i);
      double psi       = Xsig_aug(3, i);
      double psi_d     = Xsig_aug(4, i);
      double nu_a      = Xsig_aug(5, i);
      double nu_psi_dd = Xsig_aug(6, i);
      
      // regular addition
      // use '> 0.001' instead of '==0' to prevent division by zero stability?
      if (psi_d > 0.001) {
          Xsig_pred(0, i) = p_x + v/psi_d*(sin(psi + psi_d*dt) - sin(psi));
          Xsig_pred(1, i) = p_y + v/psi_d*(-cos(psi + psi_d*dt) + cos(psi));
      }
      else {
          Xsig_pred(0, i) = p_x + v*cos(psi)*dt;
          Xsig_pred(1, i) = p_y + v*sin(psi)*dt;
      }
      Xsig_pred(2, i) = v;
      Xsig_pred(3, i) = psi + psi_d * dt;
      Xsig_pred(4, i) = psi_d;
      
      // nu addition
      Xsig_pred(0, i) += 0.5 * dt_sq * cos(psi) * nu_a;
      Xsig_pred(1, i) += 0.5 * dt_sq * sin(psi) * nu_a;
      Xsig_pred(2, i) +=       dt    * nu_a;
      Xsig_pred(3, i) += 0.5 * dt_sq * nu_psi_dd;
      Xsig_pred(4, i) +=       dt    * nu_psi_dd;
  }
}

// Part A Step 3: Predict state mean and covariance
void UKF::PredictMeanAndCovariance() {
  
  //predict state mean
  MatrixXd tiled_weights      = weights_.transpose().replicate(n_x, 1);
  MatrixXd weighted_Xsig_pred = Xsig_pred.array() * tiled_weights.array();
  x_ = weighted_Xsig_pred.rowwise().sum();
  _normalize_angle(x_(3));
  
  //predict state covariance matrix
  MatrixXd norm_Xsig_pred = Xsig_pred.colwise() - x_;
  _normalize_angle_row(norm_Xsig_pred, 3);
  MatrixXd weighted_norm_Xsig_pred = norm_Xsig_pred.array() * tiled_weights.array();
  P_ = weighted_norm_Xsig_pred * norm_Xsig_pred.transpose();
  
}

// Part B Step 1a: Predict radar measurement
void UKF::PredictRadarMeasurement() {
  // initialize both predicted measurement and covariance matrix and fill with zeros
  z_pred = VectorXd(n_z);
  S_pred = MatrixXd(n_z, n_z);
  Zsig   = MatrixXd(n_z, 2*n_aug + 1);
  z_pred.setZero();
  S_pred.setZero();
  
  for (int i=0 ; i < 2*n_aug + 1 ; i++ ) {
    // store sigma points into seperate dimensions for readability
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v   = Xsig_pred(2, i);
    double psi = Xsig_pred(3, i);
    
    //transform sigma points into measurement space
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);               // rho
    Zsig(1, i) = atan2(p_y, p_x);                       // phi
    Zsig(2, i) = v*(p_x*cos(psi) + p_y*sin(psi))
                       / sqrt(p_x*p_x + p_y*p_y);       // rho_dot
    
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
  // Adding the measurement covariance noise matrix to the predicted measurement covariance matrix
  S_pred += R_radar;
}

// Part B step 1b: Predict laser measurements
void UKF::PredictLaserMeasurement() {
  // initialize both predicted measurement and covariance matrix and fill with zeros
  z_pred = VectorXd(n_z);
  S_pred = MatrixXd(n_z, n_z);
  Zsig   = MatrixXd(n_z, 2*n_aug + 1);
  z_pred.setZero();
  S_pred.setZero();
  
  // calculate mean predicted measurement
  // better for readability, but for loop will likely be more computationally efficient
  Zsig = Xsig_pred.block(0, 0, n_z_lsr, 2*n_aug+1).array();
  MatrixXd Zsig_weighted = Zsig.array().rowwise() * weights_.transpose().array();
  z_pred = Zsig_weighted.rowwise().sum();
  
  //calculate innovation(?) covariance matrix S
  VectorXd z_norm;
  for (int i=0 ; i < 2*n_aug + 1 ; i++ ) {
    z_norm = Zsig.col(i) - z_pred;
    S_pred += weights_(i) * z_norm * z_norm.transpose();
  }
  // Adding the measurement covariance noise matrix to state covariance matrix
  S_pred += R_laser;
}

// Part B Step 2: Update state mean and covariance
void UKF::UpdateState(MatrixXd &Xsig_pred, VectorXd &z_pred, VectorXd &z, MatrixXd &Zsig, MatrixXd &S_pred) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, z_pred.size());
  Tc.setZero();
  for (int i=0; i < 2*n_aug + 1; i++ ) {
    
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    _normalize_angle(x_diff(1));
    _normalize_angle(z_diff(1));
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S_pred.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = (z - z_pred);
  
  // normalizing angle
  _normalize_angle(z_diff(1));
  
  x_ += K * z_diff;
  P_ -= K * S_pred * K.transpose();
}

// initialize weight vector
void UKF::_initialize_weights(VectorXd &weights) {
  weights = VectorXd(2*n_aug + 1);
  weights(0) = lambda / (lambda + n_aug);
  weights.tail(2*n_aug).array() += 1 / (2*(lambda + n_aug));
}

// angle normalization
void UKF::_normalize_angle(double &angle) {
  angle = atan2(sin(angle), cos(angle));
}
void UKF::_normalize_angle(VectorXd &angles) {
  for (int i=0; i<angles.size(); i++) {
    _normalize_angle(angles[i]);
  }
}
void UKF::_normalize_angle_row(MatrixXd &angles, int row) {
  for (int i=0; i < angles.cols(); i++) {
    _normalize_angle(angles(row, i));
  }
}


