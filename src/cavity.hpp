#ifndef CAVITY_HPP_
#define CAVITY_HPP_

#include <eigen3/Eigen/Dense>

struct flowParams {
  int re;
  double ma;
  double pr;
  double gamma;
  double nu;
  double R;
  double uw;
};

struct integParams {
  int nx;
  int ny;
  int ngx;
  int ngy;
  int nt;
  double dx;
  double dy;
  double dt;
};

struct flowQuant {
  Eigen::MatrixXd rho;
  Eigen::MatrixXd x_mom;
  Eigen::MatrixXd y_mom;
  Eigen::MatrixXd et;
};

int rho_rhs(struct integParams integ, struct flowQuant U, Eigen::MatrixXd f);
int x_rhs(struct integParams integ, struct flowQuant U, Eigen::MatrixXd f);
int y_rhs(struct integParams integ, struct flowQuant U, Eigen::MatrixXd f);
int et_rhs(struct integParams integ, struct flowQuant U, Eigen::MatrixXd f);

int sig_diag1(struct flowParams flow, struct integParams integ, struct flowQuant U, Eigen::MatrixXd sig);
int sig_diag2(struct flowParams flow, struct integParams integ, struct flowQuant U, Eigen::MatrixXd sig);
int sig_off(struct flowParams flow, struct integParams integ, struct flowQuant U, Eigen::MatrixXd sig);

double pressure(struct flowParams flow, const int rho, const int et, const int u, const int v);

#endif // CAVITY_HPP_