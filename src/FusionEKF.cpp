#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;


  // init state vector
  VectorXd x_0(4,1);
  
  // init state covariance matrix
  MatrixXd P_0(4,4);

  // init state transition matrix
  MatrixXd F_0(4,4);

  //init process covariance matrix
  MatrixXd Q_0(4,4);

  //init  measurement matrix
  MatrixXd H_0(2,4);
  
  /*Initialise states & Matrices*/
  x_0 << 0,0,0,0;
  
  P_0 << 1,0,0,0,
         0,1,0,0,
         0,0,100,0,
         0,0,0,100;
  
  F_0 << 1,0,1,0,
         0,1,0,1,
         0,0,1,0,
         0,0,0,1;
  
  Q_0 << 1,0,1,0,
         0,1,0,1,
         1,0,1,0,
         0,1,0,1;  
 
  H_0 << 1,0,0,0,
         0,1,0,0;
 
  H_laser_ = H_0;
  
ekf_.Init(x_0, P_0, F_0, H_0,R_laser_,Q_0);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
	  /*RADAR measurement state update*/
      double r=measurement_pack.raw_measurements_(0);
	  double phi=measurement_pack.raw_measurements_(1);
	  double r_dot=measurement_pack.raw_measurements_(2);
       
	  ekf_.x_(0)=r*cos(phi);
	  ekf_.x_(1)=r*sin(phi);
	  ekf_.x_(2)=r_dot;
	  ekf_.x_(3)=0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /*Laser measurement state update*/
      ekf_.x_(0)=measurement_pack.raw_measurements_(0);
      ekf_.x_(1)=measurement_pack.raw_measurements_(1);
      ekf_.x_(2)=0;
      ekf_.x_(3)=0;
      }
	
    /*initise filter considering t=t0*/
    previous_timestamp_=measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  /*Converting us into s and determining deltaT*/
  double delta_T =(double)(measurement_pack.timestamp_- previous_timestamp_)/1000000;
  previous_timestamp_=measurement_pack.timestamp_;
 
  /*F matrix update*/
  ekf_.F_ << 1, 0, delta_T, 0,
  		0, 1, 0, delta_T,
  		0, 0, 1, 0,
  		0, 0, 0, 1;
 
  /*Q matrix update using equation Q= G*Q*Gt*/
  MatrixXd G(4,2);
  MatrixXd Qv(2,2);
  
  double deltaT_sqr2= delta_T*delta_T/2;
  
  G << deltaT_sqr2,   0,
  		0,            deltaT_sqr2,
  		delta_T,      0,
  		0,            delta_T;
  
  Qv << 9, 0,
  		0, 9;
  
  MatrixXd Gt=G.transpose();
  
  ekf_.Q_ = G*Qv*Gt;
  
  /*Predict step*/
  ekf_.Predict();
  
  
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /*Update Jocobian using state vector*/
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    
    ekf_.R_= R_radar_;
    
    /*EKF update for RADAR sensor*/
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);  

  } else {
    
    /*Kalman filter update for Laser sensor*/
    ekf_.H_ = H_laser_;
    ekf_.R_= R_laser_;
    
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
