#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

public:
  
  UKF();
  
  ~UKF();
  
  // state vector
  // preferably x_ is a private member variable, only readable through some get_x()
  // current main.cpp implementation however requires direct reading through object.x_ so it is public
  Eigen::VectorXd x_;

  // public function to process a measurment. This is all the outside world needs to call
  void ProcessMeasurement(const MeasurementPackage &meas_pack);

private:
  
  /*******************
  * Member Variables *
  *******************/
  
  // state sizes
  const int n_x;     // int to keep track of state size
  const int n_aug;   // size of augmentated state space
  int n_z;           // variable measurement space size
  const int n_z_rdr; // radar measurement space size
  const int n_z_lsr; // laser measurement space size

  // several miscellaneous variables
  bool _initialized;     // bool for initialization
  double prev_timestamp; // last measurment's timestamp for dt calculation
  double dt;             // time difference between measurements in sec
  
  // sigma calculations related variables (Part A)
  double lambda;     // spreading parameter for sigma point calculation
  VectorXd weights_; // weights for part of the prediction and update steps
  MatrixXd Xsig_aug;
  MatrixXd Xsig_pred;
  VectorXd x_pred;
  
  // Update step related variables (Part B)
  VectorXd z;      // actual raw measurements
  VectorXd z_pred; // predicted measurements
  MatrixXd S_pred; // predicted measurement covariance matrix
  MatrixXd Zsig;   // predicted sigma points in measurement space
  
  /// Standard deviations and covariance matrices ///
  // process noise nu
  double std_a;     // STD - longitudinal acceleration in m/s^2
  double std_yawdd; // STD - yaw acceleration in rad/s^2

  MatrixXd P_; // state covariance matrix
  MatrixXd R_radar; // Radar covariance matrix
  MatrixXd R_laser; // Laser covatiance matrix
  
  /*******************
  * Member functions *
  *******************/
  
  ////// Part A: PREDICTION STEP //////
  // Step 1: Generate (augmented) sigma points
  void AugmentSigmaPoints();

  // Step 2: Predict sigma points
  void PredictSigmaPoints();
  
  // Step 3: Predict state mean and covariance
  void PredictMeanAndCovariance();

  ////// Part B: UPDATE STEP //////
  // Step 1: Predict (radar/laser) measurement
  void PredictRadarMeasurement();
  void PredictLaserMeasurement();
  
  // Step 2: Update state mean and covariance
  void UpdateState();
  
  /// Helper functions ///
  void _normalize_angle(double &angle);
  void _normalize_angle(VectorXd &angles);
  void _normalize_angle_row(MatrixXd &angles, int row);
  void _initialize_weights(VectorXd &weights);

};


#endif
