#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"

class UKF {

public:
  
  UKF();
  
  ~UKF();
  
  Eigen::MatrixXd x_;
  
  void ProcessMeasurement(MeasurementPackage meas_pack);
  
  int a;
  
  void printA(Eigen::Matrix2d m);

};


#endif
