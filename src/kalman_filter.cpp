#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  MatrixXd Ft;
  
  Ft= F_.transpose();
  
  x_= F_ * x_;
  P_= F_ * P_ * Ft ;
  P_=P_+Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred(3,1);
  /*converting cartesian coordinates to polar for update step*/
  double rho=sqrt((x_(0)*x_(0))+(x_(1)*x_(1)));
  double phi=atan2(x_(1),x_(0));
  double rho_dot= ((x_(0)*x_(2))+(x_(1)*x_(3)))/rho;
  z_pred << rho,
  			phi,
  			rho_dot;
  
  VectorXd y = z - z_pred;
  
  /*Ensuring angle falls within range of -PI to PI*/
  while (y(1) > M_PI) y(1)-=2.*M_PI;
  while (y(1) < -M_PI) y(1)+=2.*M_PI;
  
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;
}
