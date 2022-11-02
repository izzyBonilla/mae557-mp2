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
  double cp;
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
  Eigen::MatrixXd u;
  Eigen::MatrixXd v;
  Eigen::MatrixXd et;
};

struct Stress {
  Eigen::MatrixXd sig11;
  Eigen::MatrixXd sig22;
  Eigen::MatrixXd south;
  Eigen::MatrixXd west;
};

struct Interps {
  double east;
  double west;
  double south;
  double north;
};

Eigen::MatrixXd rho_rhs(struct integParams integ, struct flowQuant U);
Eigen::MatrixXd x_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);
Eigen::MatrixXd y_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);
Eigen::MatrixXd et_rhs(struct flowParams flow, struct integParams integ, struct flowQuant U, struct Stress S);

Eigen::MatrixXd sig11(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::MatrixXd sig22(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::MatrixXd sig_south(struct flowParams flow, struct integParams integ, struct flowQuant U);
Eigen::MatrixXd sig_west(struct flowParams flow, struct integParams integ, struct flowQuant U);

double RT(struct flowParams flow, double et, double u, double v);

double interp2(const double q1, const double q2);
double interp4(const double q1, const double q2, const double q3, const double q4);

#endif // CAVITY_HPP_