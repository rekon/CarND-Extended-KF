#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);
  /**
  * A helper method to convert cartesian to polar.
  */
  VectorXd Cartesian2Polar(const VectorXd& x_state, float threshold = 0.0001);
  /**
  * A helper method to initially convert polar to cartesian.
  */
  VectorXd InitPolar2Cartesian(const VectorXd& x_state);
  /**
  * A helper method to get process covariance matrix Q.
  */
  MatrixXd GetProcessCovariance(const float noise_ax, const float noise_ay, const float dt = 1);

};

#endif /* TOOLS_H_ */
