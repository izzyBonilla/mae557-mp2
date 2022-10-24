#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

#include "cavity.hpp"

int main(int argc, char* argv[]) {  

  return 0;
}


/*

Gameplan:

solve:

using r as stand-in for density

dr/dt + d(ru_i)/dx_i = 0    continuity

d(ru_i)/dt + d(ru_iu_j)/dx_j = d(sigma_ij)/dx_j   momentum

d(re_t)/dt + d(ru_je_t)/dx_j = d(sigma_ij u_i)/dx_j + d/dx_j(l dT/dx_j)   total energy

p = rRT
e = 1/(1-g) 

*/