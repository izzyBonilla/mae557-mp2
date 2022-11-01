#include <iostream>
#include <fstream>
#include <cmath>

#include "cavity.hpp"

#define NGHOST 2
#define PI 100000 // 1 bar, initial pressure
#define TI 300 // K, initial temp
#define LI 1

using Eigen::MatrixXd;
using Eigen::ArrayXXd;

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
  S.sig11 = MatrixXd::Zero(integ.ngx,integ.ngy);
  S.sig22 = MatrixXd::Zero(integ.ngx,integ.ngy);
  S.south = MatrixXd::Zero(integ.ngx,integ.ngy);
  S.west = MatrixXd::Zero(integ.ngx,integ.ngy);

  // flow quant rhs vectors
  MatrixXd f_rho = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_x_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_y_mom = MatrixXd::Zero(integ.ngx,integ.ngy);
  MatrixXd f_et = MatrixXd::Zero(integ.ngx,integ.ngy);

  // first step of the lid
  // note convention: x direction is increasing with column index
  // note convention: y direction is increasing with row index

  for(int k = 1; k < integ.ngx-1; ++k) {
    U.x_mom(k,integ.ngy-2) = flow.uw*sin(flow.omega*integ.dt);
  }

  f_rho = rho_rhs(integ,U);

  S.sig11 = sig11(flow,integ, U);
  S.sig22 = sig22(flow,integ, U);
  S.south = sig_south(flow,integ,U);
  S.west = sig_west(flow,integ,U);

  std::cout << U.x_mom << std::endl << S.sig12 << std::endl;
  // S.sig2 = sig_diag2(flow,integ, U);
  // S.sig_off = sig_off(flow, integ, U);

  return 0;
}

Eigen::MatrixXd rho_rhs(struct integParams integ, struct flowQuant U) {
  // given the integrator parameters and required quantities, take
  // one timestep of the momentum equation

  MatrixXd f = MatrixXd::Zero(integ.ngx,integ.ngy);

  f = -(integ.cd*(U.x_mom/(2*integ.dx)+U.y_mom/(2*integ.dx)));

  return f;
}

Eigen::MatrixXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);

Eigen::MatrixXd sig11(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute 1-direction principal stresses on k,l grid
  // using the following scheme: every k,l index pair corresponds
  // to the stress at the western boundary; i-e sigma(1,1) corresponds to
  // sigma(1/2,1)

  MatrixXd sigma = MatrixXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd u = U.x_mom.array()/U.rho.array(); // get u velocity field
  ArrayXXd v = U.y_mom.array()/U.rho.array(); // get v velocity field

  double press; // placeholder for pressure term
  double mu;    // placeholder for dynamic viscosity term
  double et_w;
  double rho_w;
  double u_w;
  double v_nw;
  double v_sw;

  for(int k = 1; k < integ.ngx; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      rho_w = interp2(U.rho(k,l),U.rho(k-1,l));
      et_w = interp2(U.et(k,l),U.et(k-1,l));
      mu = rho_w*flow.nu; // mu at half gridpoint
      u_w = interp2(u(k,l),u(k-1,l)); // western u velocity
      press = pressure(flow,et_w,u_w,v(k,l));
      v_nw = interp4(v(k-1,l+1),v(k-1,l),v(k,l+1),v(k,l));
      v_sw = interp4(v(k-1,l),v(k-1,l-1),v(k,l),v(k,l-1));
      sigma(k,l) = mu*((4/(3*integ.dx))*(u(k,l)-u(k-1,l))-(2/(3*integ.dx))*(v_nw-v_sw)) - press;
    }
  }

  return sigma;
}

Eigen::MatrixXd sig22(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute 2-direction principal stresses on k,l grid
  // using the following scheme: every k,l index pair corresponds
  // to the stress at the southern boundary; i-e sigma(1,1) corresponds to
  // sigma(1,1/2)

  MatrixXd sigma = MatrixXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd u = U.x_mom.array()/U.rho.array(); // get u velocity field
  ArrayXXd v = U.y_mom.array()/U.rho.array(); // get v velocity field

  double press; // placeholder for pressure term
  double mu;    // placeholder for dynamic viscosity term
  double et_s;  // interpolated total energy at southern border
  double rho_s; // interpolated density at southern border
  double v_s;   // interpolated v velocity at southern border
  double u_sw;  // interpolated u velocity at southwestern corner
  double u_se;  // interpolated u velocity at southeastern corner

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy; ++l) {
      rho_s = interp2(U.rho(k,l),U.rho(k,l-1));
      et_s = interp2(U.et(k,l),U.et(k,l-1));
      mu = rho_s*flow.nu; // mu at half gridpoint
      v_s = interp2(v(k,l),v(k,l-1)); // southern v velocity
      press = pressure(flow,et_s,v_s,u(k,l));
      u_sw = interp4(u(k,l),u(k,l-1),u(k-1,l-1),u(k-1,l));
      u_se = interp4(u(k,l),u(k+1,l),u(k+1,l-1),u(k,l-1));
      sigma(k,l) = mu*((4/(3*integ.dx))*(u(k,l)-u(k-1,l))-(2/(3*integ.dx))*(u_sw-u_se)) - press;
    }
  }

  return sigma;
}

Eigen::MatrixXd sig_south(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute off-diagonall stresses on k,l grid 
  // note sigma_12 = sigma_21

  MatrixXd sigma = MatrixXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd u = U.x_mom.array()/U.rho.array();
  ArrayXXd v = U.y_mom.array()/U.rho.array();
  
  double rho_s;
  double mu;
  double v_sw;
  double v_se;

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy; ++l) {
      rho_s = interp2(U.rho(k,l),U.rho(k,l-1)); // interpolate density
      mu = rho_s*flow.nu; // calculate mu at this gridpoint
      v_sw = interp4(v(k,l),v(k-1,l),v(k-1,l-1),v(k,l-1)); // interpolate
      v_se = interp4(v(k,l),v(k+1,l),v(k+1,l-1),v(k,l-1)); // corner velocities

      // now integrate
      sigma(k,l) = mu*((u(k,l)-u(k,l-1))/integ.dy + (v_se-v_sw)/integ.dx);
    }
  }

  return sigma;
}

Eigen::MatrixXd sig_west(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  MatrixXd sigma = MatrixXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd u = U.x_mom.array()/U.rho.array();
  ArrayXXd v = U.y_mom.array()/U.rho.array();

  double rho_w;
  double mu;
  double u_nw;
  double u_sw;

  for(int k = 1; k < integ.ngx; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      rho_w = interp2(U.rho(k,l),U.rho(k-1,l));
      mu = rho_w*flow.nu;
      u_nw = interp4(u(k,l),u(k-1,l),u(k-1,l+1),u(k,l+1));
      u_sw = interp4(u(k,l),u(k,l-1),u(k-1,l-1),u(k-1,l));

      sigma(k,l) = mu*((u_nw-u_sw)/integ.dy+(v(k,l)-v(k-1,l))/integ.dx);
    }
  }

  return sigma;
}


double pressure(struct flowParams flow, double et, double u, double v) {
  /* calculate pressure given a density and total energy via EOS
  p=rhoRT
  */

  double p;

  return p;

}

double interp2(const double q1, const double q2) {
  return (q1+q2)/2;
}

double interp4(const double q1, const double q2, const double q3, const double q4) {
  return (q1+q2+q3+q4)/4;
}
