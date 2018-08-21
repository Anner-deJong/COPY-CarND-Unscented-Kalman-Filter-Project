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
  double dt; // time difference in sec  
  // measurement standart deviations
  const double std_a;
  const double std_yawdd;

  // state and state covariance matrix
  Eigen::VectorXd x;
  Eigen::MatrixXd P;

  // PREDICTION STEP
  // generate sigma points
  Eigen::MatrixXd GenerateSigmaPoints();
  
  // augment sigma points
  void AugmentSigmaPoints(Eigen::MatrixXd &Xsig_aug);

  // predict sigma points
  void PredictSigmaPoints(MatrixXd &Xsig_pred, MatrixXd &Xsig_aug);
  
  // predict mean and covariance
	
	// UPDATE STEP
	// predict measurement
	// update state


};


#endif
