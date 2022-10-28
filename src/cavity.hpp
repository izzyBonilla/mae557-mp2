#ifndef CAVITY_HPP_
#define CAVITY_HPP_

#include <eigen3/Eigen/Dense>

struct flowParams {
  int re;
  double ma;
  double pr;
  double gamma;
  double R;
  double uw;
  double nu;
  double L;
  double omega;
  double to;
};

struct integParams {
  int nx;
  int ny;
  int ngx;
  int ngy;
  int nt;
  double tf;
  double dx;
  double dy;
  double dt;
  Eigen::MatrixXd cd;
};

struct flowQuant {
  Eigen::MatrixXd rho;
  Eigen::MatrixXd x_mom;
  Eigen::MatrixXd y_mom;
  Eigen::MatrixXd et;
};

struct Stress {
  Eigen::MatrixXd sig1;
  Eigen::MatrixXd sig2;
  Eigen::MatrixXd sig_off;
};

Eigen::MatrixXd rho_rhs(struct integParams integ, struct flowQuant U);
Eigen::MatrixXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);
Eigen::MatrixXd y_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);
Eigen::MatrixXd et_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);

Eigen::MatrixXd sig_diag1(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::MatrixXd sig_diag2(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::MatrixXd sig_off(struct flowParams flow, struct integParams integ, struct flowQuant U);

double pressure(struct flowParams flow, struct flowQuant U, int k, int l, double u, double v);

#endif // CAVITY_HPP_