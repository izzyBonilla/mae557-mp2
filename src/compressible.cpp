#include <fstream>
#include <iostream>
#include <cmath>

#include "cavity.hpp"

#define TI 300 // 300 Kelvin initial temperature
#define PI 100000 // 100000 Pa
#define NGHOST 2

using Eigen::ArrayXXd;

int main(int argc, char* argv[]) {

  // define structs to carry important quantitites
  struct flowParams n; // parameters defining flow like Re and Ma
  struct integParams integ; // integrator parameters
  struct flowQuant U; // flow quantities: density, u velocity, v velocity, total energy
  struct Stress S; // stress grids

  // input flow quantities
  n.L = 1;      // TODO: Length should be variable in the future to allow for different mach
  n.ma = 0.025; // TODO: Mach will eventually be commandline input
  n.re = 100;
  n.pr = 0.7;
  n.gamma = 1.4;
  n.R = 287;
  n.uw = n.ma*sqrt(n.gamma*n.R*TI);
  n.nu = n.uw*n.L/n.re;
  n.cp = n.gamma/(n.gamma-1)*n.R;
  n.omega = pow(1/n.L,2)*2*n.nu; // * default value of this is 0.173594

  // input integrator quantities
  integ.nt = 10;
  integ.tf = 1;
  integ.dt = integ.tf/integ.nt;

  integ.nx = 10;
  integ.ngx = integ.nx + NGHOST;
  integ.ny = 10;
  integ.ngy = integ.ny + NGHOST;

  integ.dx = n.L/integ.nx;
  integ.dy = n.L/integ.ny;

  // get initial density and enery of flow field
  double rho_i = PI/(n.R*TI);
  double et_i = n.R*TI/(n.gamma-1);

  // initialize flow quantity arrays
  U.rho = ArrayXXd::Constant(integ.ngx,integ.ngy,rho_i);
  U.et = ArrayXXd::Constant(integ.ngx,integ.ngy,et_i);

  U.u = ArrayXXd::Zero(integ.ngx,integ.ngy); // * flow is quiescent to start
  U.v = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // ! Array indices indicate direction such that increasing k represents increasing
  // ! x while increasing l represents increasing y


  // ! Stress scheme is to store the -1/2th stress at k,l
  // ! i.e. sig_11|k-1/2,l is stored at S.sig11(k,l)
  // initialize stress arrays
  S.sig11 = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.sig22 = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // non-principal stresses on southern and western cell faces
  S.south = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.west = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double t;
  double wall_velo;

  // time integration
  for(int s = 0; s < integ.nt; ++s) {

    t = s*integ.dt;
    wall_velo = n.uw*sin(n.omega*t);

    // handle boundary conditions first
    // leave the corner values out of this
    // * top and bottom walls
    for(int k = 1; k < integ.ngx-1; ++k) {
      // density
      U.rho(k,integ.ngy-1) = U.rho(k,integ.ngy-2);
      U.rho(k,0) = U.rho(k,1);

      // top velocity boundary conditions
      U.u(k,integ.ngy-1) = 2*wall_velo-U.u(k,integ.ngy-2);
      U.v(k,integ.ngy-1) = -U.v(k,integ.ngy-2);

      // bottom velocity boundary conditions
      U.u(k,0) = -U.u(k,1);
      U.v(k,0) = -U.v(k,1);

      // energy boundary conditions
      U.et(k,integ.ngy-1) = 2*(et_i+pow(wall_velo,2))-U.et(k,integ.ngy-2);
      U.et(k,0) = 2*et_i-U.et(k,1); // bottom
    }
    
    // * left and right walls
    for(int l = 1; l < integ.ngy-1; ++l) {

      // density
      U.rho(integ.ngx-1,l) = U.rho(integ.ngx-2,l);
      U.rho(0,l) = U.rho(1,l);

      // right velocity boundary conditions
      U.u(integ.ngx-1,l) = U.u(integ.ngx-2,l);
      U.v(integ.ngx-1,l) = U.v(integ.ngx-2,l);

      // left velocity boundary conditions
      U.u(0,l) = U.u(1,l);
      U.v(0,l) = U.v(1,l);

      // energy
      U.et(integ.ngx-1,l) = 2*et_i-U.et(integ.ngx-2,l);
      U.et(0,l) = 2*et_i-U.et(1,l);
    }

    S.sig11 = sig11(n,integ,U);

    U.rho = U.rho + integ.dt*rho_rhs(integ,U);
    U.u = U.u + integ.dt*x_rhs(n,integ,U,S);
    U.v = U.v + integ.dt*y_rhs(n,integ,U,S);
    U.et = U.et + integ.dt*et_rhs(n,integ,U,S);

  }

  std::cout << S.sig11 << std::endl;

  return 0;
}


Eigen::ArrayXXd rho_rhs(struct integParams integ, struct flowQuant U) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      f(k,l) =
        -(U.rho(k+1,l)*U.u(k+1,l)-U.rho(k-1,l)*U.u(k-1,l))/(2*integ.dx)
        -(U.rho(k,l+1)*U.v(k,l+1)-U.rho(k,l-1)*U.v(k,l-1))/(2*integ.dy);
    }
  }

  return f;
}

Eigen::ArrayXXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      f(k,l) =
        - (U.rho(k+1,l)*pow(U.u(k+1,l),2)-U.rho(k-1,l)*pow(U.u(k-1,l),2))/(2*integ.dx)
        - (U.rho(k,l+1)*U.u(k,l+1)*U.v(k,l+1)-U.rho(k,l-1)*U.u(k,l-1)*U.v(k,l-1))/(2*integ.dy)
        + (S.sig11(k+1,l)-S.sig11(k,l))/integ.dx
        + (S.south(k,l+1)-S.south(k,l-1))/integ.dy;
    }
  }

  return f;
}

Eigen::ArrayXXd y_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  return f;
}

Eigen::ArrayXXd et_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  return f;
}


Eigen::ArrayXXd sig11(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double interp_rho;
  double pressure;
  double rt;
  double mag_v;
  double mu;

  double interp_v_plus;
  double interp_v_minus;

  for(int k = 1; k < integ.ngx; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      interp_rho = (U.rho(k-1,l)+U.rho(k,l))/2; // interpolate rho between k-1 and k to get k-1/2
      mag_v = 0.5*(pow((U.u(k-1,l)+U.u(k,l))/2,2)+pow((U.v(k-1,l)+U.v(k,l))/2,2)); // interp u and v
      rt = ((U.et(k-1,l)+U.et(k,l))/2-mag_v)*(flow.gamma-1); // RT = et(gamma-1) - 0.5*|V|^2 
      pressure = interp_rho*rt; // eos
      mu = interp_rho*flow.nu;
      interp_v_plus = (U.v(k,l)+U.v(k-1,l)+U.v(k-1,l+1)+U.v(k,l+1))/4; // top left corner v velo
      interp_v_minus = (U.v(k,l)+U.v(k,l-1)+U.v(k-1,l-1)+U.v(k-1,l))/4; // bottom left corner v velo

      sigma(k,l) = -pressure + mu * (
        4/3*(U.u(k,l)-U.u(k-1,l))/integ.dx - 2/3*(interp_v_plus-interp_v_minus)/integ.dy
      );
    }
  }

  // ! handle top corners separately
  // *  top left 
  sigma(1,integ.ngy-2) = -U.rho(1,integ.ngy-2)*flow.R*TI+U.rho(1,integ.ngy-2)*flow.nu*(
    4/3*(U.u(1,integ.ngy-2)-U.u(0,integ.ngy-2))/integ.ngx
  );

  // * top right
  sigma(integ.ngx-2,integ.ngy-2) = -U.rho(integ.ngx-2,integ.ngy-2)*flow.R*TI+U.rho(integ.ngx-2,integ.ngy-2)*flow.nu*(
    4/3*(U.u(integ.ngx-1,integ.ngy-2)-U.u(integ.ngx-2,integ.ngy-2))/integ.ngx
  );

  return sigma;
}

Eigen::ArrayXXd sig22(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  return sigma;
}

Eigen::ArrayXXd sig_south(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  return sigma;
}

Eigen::ArrayXXd sig_west(struct flowParams flow, struct integParams integ, struct flowQuant U) {

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // for(int k = 1; k < integ.ngx-1; ++k) {
  //   for(int l = 1; l < integ.ngy; ++l) {

  //   }
  // }

  return sigma;
}