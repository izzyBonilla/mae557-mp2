#ifndef CAVITY_HPP_
#define CAVITY_HPP_

struct Params {
  int nx;
  int ny;
  int ng;
  int nsteps;
  double re;
  double nu;
};

int main(int argc, char* argv[]);

#endif // CAVITY_HPP_