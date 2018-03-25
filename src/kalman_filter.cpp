#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &R_lidar_in, MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;

  P_ = MatrixXd(4,4);
  P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  F_ = MatrixXd(4,4);
  F_ << 1, 0, 0.1, 0,
			  0, 1, 0, 0.1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

  H_ = MatrixXd(2,4);
  H_ << 1, 0, 0, 0,
			  0, 1, 0, 0;
  R_lidar = R_lidar_in;
  R_radar = R_radar_in;
  
  Q_ = Q_in;
  cout<<"EKF Initialized\n";

}

void KalmanFilter::Predict() {
	cout<<" In Kalman Predict()\n";
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  cout<<"Update KF\n";
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_lidar;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  cout<<"Update EKF\n";
  //Convert radar from cartesian to polar coordinates and initialize state.
  VectorXd x_polar = tools.Cartesian2Polar(x_);
  MatrixXd Hj_ = tools.CalculateJacobian( x_ );

  VectorXd z_pred = x_polar;
	VectorXd y = z - z_pred;
  //Normalize phi
  y(1) = atan2(sin(y(1)), cos(y(1)));

	MatrixXd Hjt = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Hjt + R_radar;
	MatrixXd Si = S.inverse();
	MatrixXd PHjt = P_ * Hjt;
	MatrixXd K = PHjt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - (K * Hj_)) * P_;
}
