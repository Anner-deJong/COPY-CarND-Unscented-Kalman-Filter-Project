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
  int n_aug; // size of augmentated state space
  double lambda; // spreading parameter for sigma point calculation
  double prev_timestamp; // last measurment's timestamp for dt calculation
  double dt; // time difference between measurements in sec

  // measurement standart deviations
  const double std_a;
  const double std_yawdd;

  // state covariance matrix
  MatrixXd P_;


  ////// Part A: PREDICTION STEP //////
  // Step 1: Generating (augmented) sigma points
  MatrixXd GenerateSigmaPoints();
  void AugmentSigmaPoints(MatrixXd &Xsig_aug);

  // Step 2: Sigma point predictio
  void PredictSigmaPoints(MatrixXd &Xsig_pred, MatrixXd &Xsig_aug);
  
  // Step 3: Predict state mean and covariance
  void PredictMeanAndCovariance(MatrixXd &Xsig_pred);

  ////// Part B: UPDATE STEP //////
  // predict measurement
  // update state
  
  
  /// Helper functions ///
  void _normalize_angle(double &angle);
  void _normalize_angle(VectorXd &angles);
};


#endif
