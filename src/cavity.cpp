#include <iostream>
#include <fstream>
#include <cmath>

#include "cavity.hpp"

#define NGHOST 2
#define PI 100000 // 1 bar, initial pressure
#define TI 300 // K, initial temp

using Eigen::MatrixXd;

int main(int argc, char* argv[]) {  

  struct flowParams flow;
  struct integParams integ;
  struct flowQuant U;

  double tFinal = 1; // just a guess
  double L = 1; // TODO: MAKE THIS VARIABLE FOR FUTURE VARIATIONS IN MACH NUMBER

  // set up flow parameters
  flow.re = 100;    // Reynolds
  flow.ma = 0.25;   // Mach
  flow.pr = 0.7;    // Prandtl
  flow.gamma = 1.4; // cp/cv
  flow.nu = 0.1;    // kinematic viscosity
  flow.R = 287;     // gas const
  flow.uw = 1;


  // set up integrator parameters
  integ.nx = 3;
  integ.ny = 3;
  integ.ngx = integ.nx + NGHOST;
  integ.ngy = integ.ny + NGHOST;

  // caluculate initial density and total energy
  double rho_i = PI/(flow.R*TI);
  double et_i = flow.R*TI/(flow.gamma-1);

  // initialize flow quantity matrices
  U.rho = MatrixXd::Constant(integ.ngx,integ.ngy,rho_i);
  U.x_mom = MatrixXd::Zero(integ.ngx,integ.ngy);    // lid is stationary at start
  U.y_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  U.et = MatrixXd::Constant(integ.ngx,integ.ngy,et_i);

  // flow quant rhs vectors
  MatrixXd f_rho = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_x_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_y_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_et = MatrixXd::Zero(integ.ngx,integ.ngy);

  rho_rhs(integ,U,f_rho);

  return 0;
}

int rho_rhs(struct integParams integ, struct flowQuant U, Eigen::MatrixXd f) {
  // given the integrator parameters and required quantities, take
  // one timestep of the momentum equation

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      f(k,l) =
        (U.x_mom(k+1,l)-U.x_mom(k-1,l))/(2*integ.dx) +
        (U.y_mom(k,l+1)-U.y_mom(k,l-1))/(2*integ.dy);
    }
  }

  return 0;
}

int stress(struct flowParams flow, struct integParams integ, struct flowQuant q, int i, int j) {

  return 0;
}

double pressure(struct flowParams flow, const int rho, const int et, const int u, const int v) {
  /* calculate pressure given a density and total energy via EOS
  p=rhoRT
  */

  double mag = 0.5*(u*u+v*v);

  return rho*(et-mag)/(flow.gamma-1);
}