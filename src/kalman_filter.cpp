#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;


KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Predict() {
 
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_= F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {


  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();

  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}




void KalmanFilter::UpdateEKF(const VectorXd &z) {

  Tools tools;
  MatrixXd Hj(3,4);

  Hj = tools.CalculateJacobian(x_); 
   
  float px = x_(0);
  if (fabs(px) < 0.001)
    if (px < 0)
       px= - 0.001;
    else
       px = 0.001;

  float py = x_(1);
  if (fabs(py) < 0.001)
    if (py < 0)
       py= - 0.001;
    else
       py = 0.001;


  float vx = x_(2);
  float vy = x_(3);
  
  float c0 = px*px + py*py;
  if (c0 < 0.01)
    c0 = 0.01;
  float c1 = sqrt(c0); 
  
  float c2 = 0;
  //if (fabs(px) > 0.00001) {
  c2 = atan2(py,px);
  while (c2 > 3.14)
    c2 =c2 - 2*3.14;
  while (c2 < -3.14)
    c2 = c2 + 3.14*2;
  
  
  
  float c3 = 0;
  c3 = (px * vx + py * vy) / c1;  

  
   

  VectorXd z_pred = VectorXd(3);
  z_pred << c1, c2, c3; 
  

  VectorXd y = z - z_pred;
  
  while (y(1) > 3.14)
    y(1) =y(1) - 2*3.14;
  while (y(1) < -3.14)
    y(1) = y(1) + 3.14*2;
  


  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
 }
