#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
	
  cout << "Initialization of FusionEKF" << endl;
	
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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1; 

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  cout << "Initialization of FusionEKF complete" << endl;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  cout << "Processing measurement starts..." << endl;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  float rho     = measurement_pack.raw_measurements_(0);
      float theta    = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = rho     * cos(theta);
      ekf_.x_(1) = rho     * sin(theta);      
      ekf_.x_(2) = rho_dot * cos(theta);
      ekf_.x_(3) = rho_dot * sin(theta);
	  
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }
	
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
	cout << "Initialization with first reading complete" << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   cout << "1" << endl;
   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;
   
   float dt2 = dt*dt;
   float dt3 = dt2*dt;
   float dt4 = dt3*dt;
   cout << "2" << endl;
   ekf_.F_(0,2) = dt;
   ekf_.F_(1,3) = dt;
   
   float noise_ax = 9;
   float noise_ay = 9;
   cout << "3" << endl;
   ekf_.Q_ = MatrixXd(4, 4);
   ekf_.Q_ << dt4*noise_ax/4, 0, dt3*noise_ax/2, 0,
             0, dt4*noise_ay/4, 0, dt3*noise_ay/2,
			 dt3*noise_ax/2, 0, dt2*noise_ax, 0,
			 0, dt3*noise_ay/2, 0, dt2*noise_ay;
   cout << "4" << endl;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  cout << "5" << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	cout << "6" << endl;
    // Radar updates
	Tools t;
    Hj_ = t.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
	cout << "7" << endl;
	ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  cout << "8" << endl;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
