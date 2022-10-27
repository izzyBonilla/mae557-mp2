#include <iostream>
#include <fstream>
#include <cmath>

#include "cavity.hpp"

#define NGHOST 2
#define PI 100000 // 1 bar, initial pressure
#define TI 300 // K, initial temp
#define LI 1

using Eigen::MatrixXd;

int main(int argc, char* argv[]) {  

  struct flowParams flow;
  struct integParams integ;
  struct flowQuant U;
  struct Stress S;

  if(argc == 5) {
    // request commandline input for Ma, square gridspacing, time steps,
    // and non-dimensional final time
    flow.ma = atof(argv[1]);
    integ.nx = atoi(argv[2]);
    integ.ny = atoi(argv[2]);
    integ.nt = atoi(argv[3]); 
    flow.to = atof(argv[4]);
  } else {
    flow.ma = 0.25;
    integ.nx = 3;
    integ.ny = 3;
    integ.nt = 10;
    flow.to = 10;
  }

  // set up flow parameters
  flow.re = 100;    // Reynolds
  flow.pr = 0.7;    // Prandtl
  flow.gamma = 1.4; // cp/cv
  flow.R = 287;     // gas const
  flow.uw = flow.ma*sqrt(flow.R*flow.gamma*TI); // wall velocity
  flow.nu = flow.uw*LI/flow.re;                // kinematic viscosity
  flow.L = flow.nu*flow.re/flow.uw;           //! REVISIT LENGTH
  flow.omega = 2*flow.nu/(pow(flow.L,2));    // frequency

  // set up integrator parameters
  integ.ngx = integ.nx + NGHOST;
  integ.ngy = integ.ny + NGHOST;
  integ.tf = flow.to/flow.omega;
  integ.dx = flow.L/integ.nx;
  integ.dy = flow.L/integ.ny;
  integ.dt = integ.tf/integ.nt;
  integ.cd = MatrixXd::Zero(integ.ngx,integ.ngy);

  // Central Difference Operator
  for(int k = 1; k < integ.ngx-1; ++k) {
    integ.cd(k,k-1) = -1;
    integ.cd(k,k+1) = 1;
  }
  integ.cd(1,0) = 0;
  integ.cd(integ.ngx-2,integ.ngx-1) = 0;

  // caluculate initial density and total energy
  double rho_i = PI/(flow.R*TI);
  double et_i = flow.R*TI/(flow.gamma-1);

  // initialize flow quantity matrices
  U.rho = MatrixXd::Constant(integ.ngx,integ.ngy,rho_i);
  U.x_mom = MatrixXd::Zero(integ.ngx,integ.ngy);    // lid is stationary at start
  U.y_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  U.et = MatrixXd::Constant(integ.ngx,integ.ngy,et_i);

  // stresses
  S.sig1 = MatrixXd::Zero(integ.ngx,integ.ngy);
  S.sig2 = MatrixXd::Zero(integ.ngx,integ.ngy);
  S.sig_off = MatrixXd::Zero(integ.ngx,integ.ngx);

  // flow quant rhs vectors
  MatrixXd f_rho = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_x_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_y_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_et = MatrixXd::Zero(integ.ngx,integ.ngy);

  // first step of the lid

  for(int l = 1; l < integ.ngy-1; ++l) {
    U.x_mom(integ.ngx-2,l) = flow.uw*sin(flow.omega*integ.dt);
  }

  f_rho = rho_rhs(integ,U);

  S.sig1 = sig_diag1(flow,integ, U);
  S.sig2 = sig_diag2(flow,integ, U);
  S.sig_off = sig_off(integ, U);

  f_x_mom = x_rhs(integ, U);
  f_y_mom = y_rhs(integ, U);
  f_et = et_rhs(integ, U);

  std::cout << U.x_mom << std::endl << std::endl << f_rho << std::endl;

  return 0;
}

Eigen::MatrixXd rho_rhs(struct integParams integ, struct flowQuant U) {
  // given the integrator parameters and required quantities, take
  // one timestep of the momentum equation

  MatrixXd f = MatrixXd::Zero(integ.ngx,integ.ngy);

  f = -(integ.cd*(U.x_mom/(2*integ.dx)+U.y_mom/(2*integ.dx)));

  return f;
}

Eigen::MatrixXd sig_diag1(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute 1-direction principal stresses on k,l half grid
  // grid 0 is

  MatrixXd sigma = MatrixXd::Zero(integ.ngx,integ.ngy);
  double up;
  double ud;
  double vp;
  double vd;

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      up = (U.x_mom(k+1,l)/U.rho(k+1,l)+U.x_mom(k,l)/U.rho(k,l)); // u_(k+1/2,l)
      ud = (U.x_mom(k,l)/U.rho(k,l)+U.x_mom(k-1,l)/U.rho(k-1,l));
      sigma(k,l) = 0;
    }
  }

  return sigma;
}

Eigen::MatrixXd sig_diag2(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute 2-direction principal stresses on k,l grid

  return 0;
}

Eigen::MatrixXd sig_off(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute off-diagonall stresses on k,l grid 
  // note sigma_12 = sigma_21
  return 0;
}

double pressure(struct flowParams flow, const int rho, const int et, const int u, const int v) {
  /* calculate pressure given a density and total energy via EOS
  p=rhoRT
  */

  double mag = 0.5*(u*u+v*v);

  return rho*(et-mag)/(flow.gamma-1);
}