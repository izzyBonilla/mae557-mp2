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

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, "\t", "\n");

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
    integ.tf = atof(argv[4]);
  } else {
    flow.ma = 0.025;
    integ.nx = 25;
    integ.ny = 25;
    integ.nt = 100;
    integ.tf = 0.00001;
  }

  // set up flow parameters
  flow.re = 100;    // Reynolds
  flow.pr = 0.7;    // Prandtl
  flow.gamma = 1.4; // cp/cv
  flow.R = 287;     // gas const
  flow.uw = flow.ma*sqrt(flow.R*flow.gamma*TI); // wall velocity
  flow.nu = flow.uw*LI/flow.re;                // kinematic viscosity
  flow.L = flow.nu*flow.re/flow.uw;           //! REVISIT LENGTH
  flow.omega = 2*flow.nu/(pow(flow.L,2));    // frequency, default 0.1735954
  flow.cp = flow.gamma/(1-flow.gamma)*flow.R;

  // set up integrator parameters
  integ.ngx = integ.nx + NGHOST;
  integ.ngy = integ.ny + NGHOST;
  integ.dx = flow.L/integ.nx;
  integ.dy = flow.L/integ.ny;
  integ.dt = integ.tf/integ.nt;

  // caluculate initial density and total energy
  double rho_i = PI/(flow.R*TI);
  double et_i = flow.R*TI/(flow.gamma-1);

  // initialize flow quantity matrices
  U.rho = ArrayXXd::Constant(integ.ngx,integ.ngy,rho_i);
  U.u = ArrayXXd::Zero(integ.ngx,integ.ngy);    // lid is stationary at start
  U.v = ArrayXXd::Zero(integ.ngx,integ.ngy);
  U.et = ArrayXXd::Constant(integ.ngx,integ.ngy,et_i);

  // stresses
  // Note on convention: the stresses are calculated such that 
  // sig11(k,l) corresponds to sigma_11|k-1/2,l
  // sig22(k,l) corresponds to sigma_11|k,l-1/2
  // south(k,l) corresponds to sigma_12|k-1/2,l
  // west(k,l)  corresponds to sigma_21|k,l-1/2 
  S.sig11 = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.sig22 = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.south = ArrayXXd::Zero(integ.ngx,integ.ngy);
  S.west = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // flow quant rhs vectors
  ArrayXXd f_rho = ArrayXXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd f_x_mom = ArrayXXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd f_y_mom = ArrayXXd::Zero(integ.ngx,integ.ngy);
  ArrayXXd f_et = ArrayXXd::Zero(integ.ngx,integ.ngy);

  // first step of the lid
  // note convention: x direction is increasing with column index
  // note convention: y direction is increasing with row index

  double t;
  double e_ghost;

  double wall_velo;

  for(int s = 0; s < integ.nt; ++s) {

    t = integ.dt*s;
    wall_velo = flow.uw*sin(flow.omega*t);

    // unknowns at cell centers
    // horizontal wall boundary conditions

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
      U.et(k,integ.ngy-1) = (flow.R*TI-RT(flow,U.et(k,integ.ngy-2),U.u(k,integ.ngy-2),U.v(k,integ.ngy-2)))/(flow.gamma-1)
                            + 0.5*(pow(U.u(k,integ.ngy-1),2)+pow(U.v(k,integ.ngy-1),2));

      U.et(k,0) = (flow.R*TI-RT(flow,U.et(k,1),U.u(k,1),U.v(k,1)))/(flow.gamma-1)
                            + 0.5*(pow(U.u(k,0),2)+pow(U.v(k,0),2)); // bottom
    }
    
    // * left and right walls
    for(int l = 1; l < integ.ngy-1; ++l) {

      // density
      U.rho(integ.ngx-1,l) = U.rho(integ.ngx-2,l);
      U.rho(0,l) = U.rho(1,l);

      // right velocity boundary conditions
      U.u(integ.ngx-1,l) = -U.u(integ.ngx-2,l);
      U.v(integ.ngx-1,l) = -U.v(integ.ngx-2,l);

      // left velocity boundary conditions
      U.u(0,l) = -U.u(1,l);
      U.v(0,l) = -U.v(1,l);

      // energy
      U.et(integ.ngx-1,l) = (flow.R*TI-RT(flow,U.et(integ.ngx-2,l),U.u(integ.ngx-2,l),U.v(integ.ngx-2,l)))/(flow.gamma-1)
                            + 0.5*(pow(U.u(integ.ngx-1,l),2)+pow(U.v(integ.ngx-1,l),2));
      
      U.et(0,l) = (flow.R*TI-RT(flow,U.et(1,l),U.u(1,l),U.v(1,l)))/(flow.gamma-1)
                  + 0.5*(pow(U.u(0,l),2)+pow(U.v(0,l),2));
    }

    S.sig11 = sig11(flow,integ, U);
    // S.sig22 = sig22(flow,integ, U);
    // S.south = sig_south(flow,integ,U);
    // S.west = sig_west(flow,integ,U);

    f_rho = rho_rhs(integ,U);
    f_x_mom = x_rhs(flow,integ,U,S);
    f_y_mom = y_rhs(flow,integ,U,S);
    f_et = et_rhs(flow,integ,U,S);

    U.rho = U.rho + integ.dt*f_rho;
    U.u = U.u + integ.dt*f_x_mom/U.rho;
    U.v = U.v + integ.dt*f_y_mom/U.rho;

    // U.et = U.et + integ.dt*f_et;

  }

  // write density to file
  std::ofstream rho_file("rho.csv");
  if(rho_file.is_open()) {
      rho_file << S.sig11.format(CSVFormat);
  } else {
      std::cout << "Can't write to file rho.csv \n";
      exit(1);
  }

  // write u to file
  std::ofstream u_file("u.csv");
  if(u_file.is_open()) {
      u_file << U.u; // .format(CSVFormat);
  } else {
      std::cout << "Can't write to file u.csv \n";
      exit(1);
  }

  // write v to file
  std::ofstream v_file("v.csv");
  if(v_file.is_open()) {
      v_file << U.v.format(CSVFormat);
  } else {
      std::cout << "Can't write to file v.csv \n";
      exit(1);
  }

  // write et to file
  std::ofstream et_file("et.csv");
  if(et_file.is_open()) {
      et_file << U.et.format(CSVFormat);
  } else {
      std::cout << "Can't write to file et.csv \n";
      exit(1);
  }

  return 0;
}

Eigen::ArrayXXd rho_rhs(struct integParams integ, struct flowQuant U) {
  // given the integrator parameters and required quantities, take
  // one timestep of the momentum equation

  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      f(k,l) = -(U.rho(k+1,l)*U.u(k+1,l)-U.rho(k-1,l)*U.u(k-1,l))/(2*integ.dx)
               -(U.rho(k,l+1)*U.v(k,l+1)-U.rho(k,l-1)*U.v(k,l-1))/(2*integ.dy);
    }
  }

  return f;
}

Eigen::ArrayXXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      f(k,l) = (S.sig11(k+1,l)-S.sig11(k,l))/integ.dx + (S.south(k,l+1)-S.south(k,l))/integ.dy -
               (U.rho(k+1,l)*pow(U.u(k+1,l),2)-U.rho(k-1,l)*pow(U.u(k-1,l),2))/(2*integ.dx) -
               (U.rho(k,l+1)*U.u(k,l+1)*U.v(k,l+1)-U.rho(k,l-1)*U.u(k,l-1)*U.v(k,l-1))/(2*integ.dy);              
    }
  }

  return f;
}

Eigen::ArrayXXd y_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {
  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      f(k,l) = (S.west(k+1,l)-S.west(k,l))/integ.dx + (S.sig22(k,l+1)-S.sig22(k,l))/integ.dy -
               (U.rho(k+1,l)*U.u(k+1,l)*U.v(k+1,l)-U.rho(k-1,l)*U.u(k-1,l)*U.v(k-1,l))/(2*integ.dx) -
               (U.rho(k,l+1)*pow(U.v(k,l+1),2)-U.rho(k,l-1)*pow(U.v(k,l-1),2))/(2*integ.dy);
    }
  }

  return f;
}

Eigen::ArrayXXd et_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S) {

  ArrayXXd f = ArrayXXd::Zero(integ.ngx,integ.ngy);

  struct Interps u_i;
  struct Interps v_i;
  struct Interps lambda;

  double Tkp;
  double Tkm;
  double Tc;
  double Tlp;
  double Tlm;

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      

      // collect interpolated velocities
      u_i.east = interp2(U.u(k,l),U.u(k+1,l));
      u_i.north = interp2(U.u(k,l),U.u(k,l+1));
      u_i.south = interp2(U.u(k,l),U.u(k,l-1));
      u_i.west = interp2(U.u(k,l),U.u(k-1,l));

      v_i.east = interp2(U.v(k,l),U.v(k+1,l));
      v_i.north = interp2(U.v(k,l),U.v(k,l+1));
      v_i.south = interp2(U.v(k,l),U.v(k,l-1));
      v_i.west = interp2(U.v(k,l),U.v(k-1,l));

      lambda.east = interp2(U.rho(k,l),U.rho(k+1,l))*flow.cp*flow.nu/flow.pr;
      lambda.north = interp2(U.rho(k,l),U.rho(k,l+1))*flow.cp*flow.nu/flow.pr;
      lambda.south = interp2(U.rho(k,l),U.rho(k,l-1))*flow.cp*flow.nu/flow.pr;
      lambda.west = interp2(U.rho(k,l),U.rho(k-1,l))*flow.cp*flow.nu/flow.pr;

      Tkp = RT(flow,U.et(k+1,l),U.u(k+1,l),U.v(k+1,l));
      Tkm = RT(flow,U.et(k-1,l),U.u(k-1,l),U.v(k-1,l));
      Tc = RT(flow,U.et(k,l),U.u(k,l),U.v(k,l));
      Tlp = RT(flow,U.et(k,l+1),U.u(k,l+1),U.v(k,l+1));
      Tlm = RT(flow,U.et(k,l-1),U.u(k,l-1),U.v(k,l-1));


      f(k,l) = 
        (S.sig11(k+1,l)*u_i.east+S.west(k+1,l)*v_i.east-S.sig11(k,l)*u_i.west+S.west(k,l)*v_i.west)/integ.dx +
        (S.south(k,l+1)*u_i.north+S.sig22(k,l+1)*v_i.north-S.south(k,l)*u_i.south+S.sig22(k,l)*v_i.south)/integ.dy +
        (Tkp*lambda.east-Tc*lambda.west-Tc*lambda.east+Tkm*lambda.south)/(integ.dx*integ.dx) +
        (Tlp*lambda.north - Tc*lambda.south - Tc*lambda.north + Tlm*lambda.south)/(integ.dy*integ.dy) -
        (U.rho(k+1,l)*U.u(k+1,l)*U.et(k+1,l)-U.rho(k-1,l)*U.u(k-1,l)*U.et(k-1,l))/(2*integ.dx) -
        (U.rho(k,l+1)*U.u(k,l+1)*U.et(k,l+1)-U.rho(k,l-1)*U.u(k,l-1)*U.et(k,l-1))/(2*integ.dy);
    }
  }

  return f;
}

Eigen::ArrayXXd sig11(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute 1-direction principal stresses on k,l grid
  // using the following scheme: every k,l index pair corresponds
  // to the stress at the western boundary; i-e sigma(1,1) corresponds to
  // sigma(1/2,1)

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double press; // placeholder for pressure term
  double mu;    // placeholder for dynamic viscosity term
  double et_w;
  double rho_w;
  double u_w;
  double v_w;
  double v_nw;
  double v_sw;

  for(int k = 1; k < integ.ngx; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      rho_w = interp2(U.rho(k,l),U.rho(k-1,l));
      et_w = interp2(U.et(k,l),U.et(k-1,l));
      mu = rho_w*flow.nu; // mu at half gridpoint
      u_w = interp2(U.u(k,l),U.u(k-1,l)); // western u velocity
      v_w = interp2(U.v(k,l),U.v(k-1,l)); // western v velocity
      press = rho_w*RT(flow,et_w,u_w,v_w);
      v_nw = interp4(U.v(k-1,l+1),U.v(k-1,l),U.v(k,l+1),U.v(k,l));
      v_sw = interp4(U.v(k-1,l),U.v(k-1,l-1),U.v(k,l),U.v(k,l-1));
      sigma(k,l) = mu*((4/(3*integ.dx))*(U.u(k,l)-U.u(k-1,l))-(2/(3*integ.dx))*(v_nw-v_sw)) - press;
    }
  }

  // top left corner stress
  rho_w = interp2(U.rho(1,integ.ngy-2),U.rho(0,integ.ngy-2));
  et_w = interp2(U.et(1,integ.ngy-2),U.et(0,integ.ngy-2));
  mu = rho_w*flow.nu;
  press = rho_w*RT(flow,et_w,0,0);
  sigma(1,integ.ngy-2) = mu*4/(3*integ.dx)*(U.u(1,integ.ngy-2)-U.u(0,integ.ngy-2))-press;

  // top right corner stress
  rho_w = interp2(U.rho(integ.ngx-1,integ.ngy-2),U.rho(integ.ngx-2,integ.ngy-2));
  et_w = interp2(U.et(integ.ngx-1,integ.ngy-2),U.et(integ.ngx-2,integ.ngy-2));
  mu = rho_w*flow.nu;
  press = rho_w*RT(flow,et_w,0,0);
  sigma(integ.ngx-1,integ.ngy-2) = mu*4/(3*integ.dx)*(U.u(integ.ngx-1,integ.ngy-2)-U.u(integ.ngx-2,integ.ngy-2))-press;

  return sigma;
}

Eigen::ArrayXXd sig22(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute 2-direction principal stresses on k,l grid
  // using the following scheme: every k,l index pair corresponds
  // to the stress at the southern boundary; i-e sigma(1,1) corresponds to
  // sigma(1,1/2)

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double press; // placeholder for pressure term
  double mu;    // placeholder for dynamic viscosity term
  double et_s;  // interpolated total energy at southern border
  double rho_s; // interpolated density at southern border
  double u_s;   // interpolated u velocity at southern border
  double v_s;   // interpolated v velocity at southern border
  double u_sw;  // interpolated u velocity at southwestern corner
  double u_se;  // interpolated u velocity at southeastern corner

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy; ++l) {
      rho_s = interp2(U.rho(k,l),U.rho(k,l-1));
      et_s = interp2(U.et(k,l),U.et(k,l-1));
      mu = rho_s*flow.nu; // mu at half gridpoint
      u_s = interp2(U.u(k,l),U.u(k,l-1)); // southern u velocity
      v_s = interp2(U.v(k,l),U.v(k,l-1)); // southern v velocity
      press = rho_s*RT(flow,et_s,v_s,u_s);
      u_sw = interp4(U.u(k,l),U.u(k,l-1),U.u(k-1,l-1),U.u(k-1,l));
      u_se = interp4(U.u(k,l),U.u(k+1,l),U.u(k+1,l-1),U.u(k,l-1));
      sigma(k,l) = mu*((4/(3*integ.dy))*(U.v(k,l)-U.v(k-1,l))-(2/(3*integ.dx))*(u_sw-u_se)) - press;
    }
  }

  // top left corner stress
  rho_s = interp2(U.rho(1,integ.ngy-1),U.rho(1,integ.ngy-2));
  et_s = interp2(U.et(1,integ.ngy-1),U.et(1,integ.ngy-2));
  mu = rho_s*flow.nu;
  u_s = interp2(U.u(1,integ.ngy-1),U.u(1,integ.ngy-2));
  v_s = interp2(U.v(1,integ.ngy-1),U.v(1,integ.ngy-2));
  press = rho_s*RT(flow,et_s,v_s,u_s);
  sigma(1,integ.ngy-1) = 4*mu/(3*integ.dy)*(U.v(1,integ.ngy-1)-U.v(1,integ.ngy-2))- press;

  // top right corner stress
  rho_s = interp2(U.rho(1,integ.ngy-1),U.rho(1,integ.ngy-2));
  et_s = interp2(U.et(1,integ.ngy-1),U.et(1,integ.ngy-2));
  mu = rho_s*flow.nu;
  u_s = interp2(U.u(1,integ.ngy-1),U.u(1,integ.ngy-2));
  v_s = interp2(U.v(1,integ.ngy-1),U.v(1,integ.ngy-2));
  press = rho_s*RT(flow,et_s,v_s,u_s);
  sigma(integ.ngx-2,integ.ngy-1) = 4*mu/(3*integ.dy)*(U.v(1,integ.ngy-1)-U.v(1,integ.ngy-2))- press;

  return sigma;
}

Eigen::ArrayXXd sig_south(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  // compute off-diagonall stresses on k,l grid 
  // note sigma_12 = sigma_21

  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);
  
  double rho_s;
  double mu;
  double v_sw;
  double v_se;

  for(int k = 1; k < integ.ngx-1; ++k) {
    for(int l = 1; l < integ.ngy; ++l) {
      rho_s = interp2(U.rho(k,l),U.rho(k,l-1)); // interpolate density
      mu = rho_s*flow.nu; // calculate mu at this gridpoint
      v_sw = interp4(U.v(k,l),U.v(k-1,l),U.v(k-1,l-1),U.v(k,l-1)); // interpolate
      v_se = interp4(U.v(k,l),U.v(k+1,l),U.v(k+1,l-1),U.v(k,l-1)); // corner velocities

      // now integrate
      sigma(k,l) = mu*((U.u(k,l)-U.u(k,l-1))/integ.dy + (v_se-v_sw)/integ.dx);
    }
  }

  // top left corner stress
  rho_s = interp2(U.rho(1,integ.ngy-2),U.rho(0,integ.ngy-2));
  mu = rho_s*flow.nu;
  sigma(1,integ.ngy-2) = mu/integ.dx*(U.v(1,integ.ngy-2),U.v(0,integ.ngy-2));

  // top right corner stress
  rho_s = interp2(U.rho(integ.ngx-1,integ.ngy-2),U.rho(integ.ngx-2,integ.ngy-2));
  mu = rho_s*flow.nu;
  sigma(integ.ngx-1,integ.ngy-2) = mu/integ.dx*(U.v(integ.ngx-1,integ.ngy-2),U.v(integ.ngx-2,integ.ngy-2));

  return sigma;
}

Eigen::ArrayXXd sig_west(struct flowParams flow, struct integParams integ, struct flowQuant U) {
  ArrayXXd sigma = ArrayXXd::Zero(integ.ngx,integ.ngy);

  double rho_w;
  double mu;
  double u_nw;
  double u_sw;

  for(int k = 1; k < integ.ngx; ++k) {
    for(int l = 1; l < integ.ngy-1; ++l) {
      rho_w = interp2(U.rho(k,l),U.rho(k-1,l));
      mu = rho_w*flow.nu;
      u_nw = interp4(U.u(k,l),U.u(k-1,l),U.u(k-1,l+1),U.u(k,l+1));
      u_sw = interp4(U.u(k,l),U.u(k,l-1),U.u(k-1,l-1),U.u(k-1,l));
      sigma(k,l) = mu*((u_nw-u_sw)/integ.dy+(U.v(k,l)-U.v(k-1,l))/integ.dx);
    }
  }

  // top left corner stress
  rho_w = interp2(U.rho(1,integ.ngy-1),U.rho(1,integ.ngy-2));
  mu = rho_w*flow.nu;
  sigma(1,integ.ngy-1) = mu/integ.dx*(U.v(1,integ.ngy-1),U.v(1,integ.ngy-2));

  // top left corner stress
  rho_w = interp2(U.rho(integ.ngx-2,integ.ngy-1),U.rho(integ.ngx-2,integ.ngy-2));
  mu = rho_w*flow.nu;
  sigma(integ.ngx-2,integ.ngy-2) = mu/integ.dx*(U.v(integ.ngx-2,integ.ngy-1),U.v(integ.ngx-2,integ.ngy-2));

  return sigma;
}

double RT(struct flowParams flow, double et, double u, double v) {
  // return the value of RT from total energy and velocity
  return (et-0.5*(u*u+v*v))*(flow.gamma-1);
}

double interp2(const double q1, const double q2) {
  return (q1+q2)/2;
}

double interp4(const double q1, const double q2, const double q3, const double q4) {
  return (q1+q2+q3+q4)/4;
}
