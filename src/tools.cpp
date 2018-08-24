#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size

  if (estimations.size() != ground_truth.size()) {
    cout << "CalculateRMSE() - Error - estimation size != ground truth size" << endl;
    return rmse;
  }
  // assuming all entries in either array have the same number of elements
  // thus we only need to check if the first entries of both arrays have equal size:
  float en = estimations[0].size(); // number of elements (n) in current (i) estimation
  float gn = ground_truth[0].size(); // number of elements (n) in current (i) ground truth

  if (en != gn) {
    cout << "CalculateRMSE() - Error - estimation[0] size != ground truth[0] size" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){

    // from Udacity answer, one less for-loop
    VectorXd difference = estimations[i] - ground_truth[i];
    VectorXd residual   = difference.array().square();
    rmse += residual;
  }

  //calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse  = rmse.array().sqrt();

  return rmse;

}
