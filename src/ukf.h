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
  
  // preferably x_ is a private member variable, only readable through some get_x()
  // current main.cpp implementation however requires direct reading through object.x_ so it is public
  // state matrix
  Eigen::VectorXd x_;

  int a;
  
  // functions
  void ProcessMeasurement(const MeasurementPackage &meas_pack);
  
  void printA(Eigen::Matrix2d m);

private:

	// include all other required matrices
  bool _initialized; // bool for initialization
  const int n_x; // int to keep track of state size
  const int n_aug; // size of augmentated state space
  const int n_z_rdr; // radar measurement space size
  const int n_z_lsr; // laser measurement space size

  double lambda; // spreading parameter for sigma point calculation
  double prev_timestamp; // last measurment's timestamp for dt calculation
  double dt; // time difference between measurements in sec
  VectorXd weights_; // weights for part of the prediction and update steps
  
  /// Standard deviations and covariance matrices ///
  // process noise nu
  double std_a;
  double std_yawdd;

  MatrixXd P_; // state covariance matrix
  MatrixXd R_radar; // Radar covariance matrix
  MatrixXd R_laser; // Laser covatiance matrix

  ////// Part A: PREDICTION STEP //////
  // Step 1: Generate (augmented) sigma points
  MatrixXd GenerateSigmaPoints();
  void AugmentSigmaPoints(MatrixXd &Xsig_aug);

  // Step 2: Predict sigma points
  void PredictSigmaPoints(MatrixXd &Xsig_pred, MatrixXd &Xsig_aug);
  
  // Step 3: Predict state mean and covariance
  void PredictMeanAndCovariance(MatrixXd &Xsig_pred);

  ////// Part B: UPDATE STEP //////
  // Step 1: Predict (radar/laser) measurement
  void PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &Zsig, MatrixXd &S_pred, MatrixXd &Xsig_pred);
  void PredictLaserMeasurement(VectorXd &z_pred, MatrixXd &Zsig, MatrixXd &S_pred, MatrixXd &Xsig_pred);
  
  // Step 2: Update state mean and covariance
  void UpdateState(MatrixXd &Xsig_pred, VectorXd &z_pred, VectorXd &z, MatrixXd &Zsig, MatrixXd &S_pred);
  
  /// Helper functions ///
  void _normalize_angle(double &angle);
  void _normalize_angle(VectorXd &angles);
  void _initialize_weights(VectorXd &weights);
};


#endif
