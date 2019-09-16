#include "kalman_filter.h"
#include <iostream>

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
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
//   std::cout << "KalmanFilter Update Enter:" << std::endl;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

//   std::cout << "KalmanFilter Update Tack0:" << std::endl;     
  //new estimate
  x_ = x_ + (K * y);
//   std::cout << "KalmanFilter Update Tack1:" << std::endl;     
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
//   std::cout << "KalmanFilter Update Tack2:" << std::endl;     
  P_ = (I - K * H_) * P_;
//   std::cout << "KalmanFilter Update Exit:" << std::endl;     
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
//   std::cout << "KalmanFilter UpdateEKF Enter:" << std::endl;
  //the function h(x) maps values from Cartesian coordinates to polar coordinates
  VectorXd z_pred = VectorXd(3);
   // float c = sqrt(x_(0)*x(0) + x(1)*x(1))
   // z_pred(0) = c;
   // z_pred(1) = atan(x(1)/x(0));
   // z_pred(1) = (x(0)*x(2) + x(1)*x(3))/c;
   // z_pred << c,atan2(x(1),x(0)),(x(0)*x(2) + x(1)*x(3))/c
   float px = x_[0];
   float py = x_[1];
   float vx = x_[2];
   float vy = x_[3];
   float rho;
   float phi;
   float rho_dot;
  
//   std::cout << "KalmanFilter UpdateEKF Tack0:" << std::endl;
  // if px and py are small at any point, set phi and rho_dot to zero
  if(fabs(px) < 0.0001 || fabs(py) < 0.0001)
  {
    if(fabs(px) < 0.0001){
      px = 0.0001; //set px to 0.0001
    }

    if(fabs(py) < 0.0001){
      py = 0.0001;   //set py to 0.0001
    }

    rho = sqrt(px*px + py*py);
    phi = 0;  //set phi to zero
    rho_dot = 0; //set rho_dot to zero
  }
  else
  {
    rho = sqrt(px*px + py*py);
    phi = atan2(py,px); 
    rho_dot = (px*vx + py*vy) /rho;
  } 
  z_pred << rho, phi, rho_dot;
  
//   std::cout << "KalmanFilter UpdateEKF Tack1:" << std::endl; 
//   std::cout << "z = " << z << std::endl;
//   std::cout << "z_pred =" << z_pred << std::endl;  
  VectorXd y = z - z_pred;
//   std::cout << "KalmanFilter UpdateEKF Tack2:" << std::endl;
  
  phi = y(1);
  while(fabs(phi) > 6.2831)
  {
    if(phi > 6.2831)
    {
      phi -= 6.2831;
    }
    else
    {
      phi += 6.2831;
    }
  }
  // normalize ϕ in the y vector so that its angle is between −π and π
  if(phi > 3.14159)
  {
    phi -=6.2831;
  }
  else if(phi < -3.14159)
  {
    phi +=6.2831;
  }
  y(1) = phi;
    // bool in_range = false;
    // while (in_range == false) {
    //   if (y(1) > 3.14159) {
    //     y(1) = y(1) - 6.2831;
    //   }
    //   else if (y(1) < -3.14159) {
    //     y(1) = y(1) + 6.2831;
    //   } 
    // else {
    //     in_range = true;
    //   }
    // }
      
  
  
//   std::cout << "KalmanFilter UpdateEKF Tack3:" << std::endl; 
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

//   std::cout << "KalmanFilter UpdateEKF Tack4:" << std::endl; 
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
//   std::cout << "KalmanFilter UpdateEKF Exit:" << std::endl;   
}
