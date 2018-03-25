#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0,0,0,0;

	if( estimations.size() == 0 || estimations.size() != ground_truth.size()){
    std::cout<<"Zero Estimations or ground_truth size and estimations size not equal";
    return rmse;
	}

    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array()*diff.array();
		rmse += diff;
	}
	
	//calculate the mean	
	rmse /= estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  cout<<"In CalculateJacobian()\n";
  MatrixXd Hj = MatrixXd::Zero(3, 4);
	//recover state parameters
	double px = x_state(0),
				 py = x_state(1),
				 vx = x_state(2),
				 vy = x_state(3),

  			 threshold = 0.0001,
  			 sq_sum = (px*px + py*py),
				 sqrt_sq_sum = sqrt(sq_sum),
				 cube_sqrt_sum = sqrt_sq_sum * sq_sum;

	//check division by zero
	if(( px == 0 && py == 0 ) || (sqrt_sq_sum < threshold)){
    cout<<"Division by zero error";
	} else {
    //compute the Jacobian matrix
    Hj << px/sqrt_sq_sum, py/sqrt_sq_sum, 0, 0,
        -py/sq_sum, px/sq_sum, 0, 0,
        (py*(vx*py - vy*px))/cube_sqrt_sum, (px*(vy*px -vx*py))/cube_sqrt_sum, px/sqrt_sq_sum, py/sqrt_sq_sum;
	}
	

	return Hj;
}

VectorXd Tools::Cartesian2Polar(const VectorXd& x_state, float threshold){
  cout<<"In Cartesian2Polar\n";
  VectorXd x_polar(3);
	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

  double rho = sqrt((px*px + py*py));
  double phi = (py == 0.0 && px == 0.0) ? 0.0 : atan2( py, px );
  double drho = (rho > threshold) ? ( px * vx + py * vy ) / rho : 0.0;
	
  x_polar << rho, phi, drho;
	return x_polar;
}

VectorXd Tools::InitPolar2Cartesian(const VectorXd& x_state){
  cout<<"In InitPolar2Cartesian\n";
  VectorXd x_cartesian(4);
	//recover state parameters
	double rho = x_state(0);
	double phi = x_state(1);
  double drho = x_state(2);

  double px = rho*cos( phi );
  double py = rho*sin( phi );

	
  x_cartesian << px, py, 0, 0;
	return x_cartesian;
}

MatrixXd Tools::GetProcessCovariance(const float noise_ax, const float noise_ay, const float dt){
  cout<<"In GetProcessCovariance()\n";
  double 
			dt2 = dt*dt,
			dt3 = dt2*dt,
			dt4 = dt3*dt,
			dt4ax = (dt4*noise_ax)/4,
	    dt4ay = (dt4*noise_ay)/4,
	    dt3ax = (dt3*noise_ax)/2,
	    dt3ay = (dt3*noise_ay)/2,
	    dt2ax = (dt2*noise_ax),
	    dt2ay = (dt2*noise_ay);
  
	MatrixXd Q_ = MatrixXd(4, 4);
	Q_ << dt4ax, 0, dt3ax, 0,
				0, dt4ay, 0, dt3ay,
				dt3ax, 0, dt2ax, 0,
				0, dt3ay, 0, dt2ay;
  return Q_;
}
