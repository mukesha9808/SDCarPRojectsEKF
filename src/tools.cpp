#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 
  VectorXd rmse(4);
  rmse <<  0,0,0,0;
  
  
  /*Estimation szie is not zero and size of estimtion and ground truth are  same
    for then only calculate RMSE 
    rmse =  (((estimation - ground truth)^2)/n)^0.5*/
  if ((estimations.size() != 0) && (estimations.size()==ground_truth.size())) {
		for(unsigned int i=0; i< estimations.size() ;++i) {
		  VectorXd residue;
		  residue=estimations[i]-ground_truth[i];
		  residue= residue.array() * residue.array() ;
          rmse += residue;
		}
    rmse=rmse/estimations.size() ;
    rmse= rmse.array().sqrt();
	} else {
		std::cout << "Inalid input data " << std::endl;
    }
	return rmse;
}


/* Jacobian Calculation*/
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd h_j(3,4);
  
  double px=x_state(0);
  double py=x_state(1);
  double vx=x_state(2);
  double vy=x_state(3);
  
  /* term Px^2 +py^2 */
  double p_ressqr= px*px + py*py;
  /* term  (Px^2 +py^2)^(1/2)*/
  double p_res=sqrt(p_ressqr);
  
  double h11,h12,h21,h22, h31,h32,h3_temp;
  if  (p_res < 1e-5) {
  /*term px/((Px^2 +py^2)^(1/2))*/
  h11= 1e-5;
  /*term py/((Px^2 +py^2)^(1/2))*/
  h12= 1e-5;
  
  /*term -py/((Px^2 +py^2))*/
  h21= -1e-5;
  /*term px/((Px^2 +py^2))*/
  h22= 1e-5;
  } else{
   /*term px/((Px^2 +py^2)^(1/2))*/
  h11= px/p_res;
  /*term py/((Px^2 +py^2)^(1/2))*/
  h12= py/p_res;
  
  /*term -py/((Px^2 +py^2))*/
  h21= -py/p_ressqr;
  /*term px/((Px^2 +py^2))*/
  h22= px/p_ressqr; 
    
  }
  
 /*term (vypx - vxpy) */
  h3_temp =(vy*h11 - vx*h12);
  /*term py(vxpy - vypx)
        --------------------------
         (Px^2 +py^2)^(3/2) */
  h31= h3_temp*h21;
  
  /*term px(vypx - vxpy)
        --------------------------
         (Px^2 +py^2)^(3/2) */
  h32= h3_temp*h21;
  
 h_j <<  h11, h12, 0, 0,
  		 h21, h22, 0, 0,
  		 h31, h32, h11, h12 ;
  
 return h_j;
  
  
  
}
