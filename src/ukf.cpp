#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"

//ctor
UKF::UKF() {
  
}

//dtor
UKF::~UKF() {

}

void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  Eigen::MatrixXd m(1,2);
  m << 3, 4;
  printA((Eigen::Matrix2d() << 1, 2, 3, 4).finished());
}
// print function
void UKF::printA(Eigen::Matrix2d m) {
  a = 3;
  std::cout << "m = " << m << std::endl << std::endl;
}


