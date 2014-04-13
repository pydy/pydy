#include <iostream>
#include <fstream>
#include "rattleback.h"
#include <Eigen/Dense>

int main(int argc, char *argv[]) {
  using namespace Eigen;
  rattleback_params p;
  p.a = 0.2;
  p.b = 0.03;
  p.c = 0.02;
  p.d = p.e = 0.0;
  p.f = 0.01;
  p.m = 1.0;
  p.g = 9.81;
  p.Ixx =  0.0002;
  p.Iyy =  0.0016;
  p.Izz =  0.0017;
  p.Ixy = -0.00002;
  p.Iyz = p.Ixz = 0.0;
  p.s = 1e-4;

  // Initial state
  double state[] = {0.0,   // Yaw
                    0.0,   // Roll
                    0.0,   // Pitch
                    0.0,   // x
                    0.0,   // y
                    0.0,   // u0
                    0.0,   // u1
                    0.0};  // u2
  double dfdx[64];
  double dfdt[8];
  Matrix<double, 4, 4> J;

  int N = 11;
  double spin = 4.0;
  for (int k = 0; k < N; ++k) {
    state[7] = -spin + 2*spin/(N - 1) * k;
    rattleback_jacobian(0.0, state, dfdx, dfdt, &p);
    J(0, 0) = dfdx[8*1 + 1];
    J(0, 1) = dfdx[8*1 + 2];
    J(0, 2) = dfdx[8*1 + 5];
    J(0, 3) = dfdx[8*1 + 6];
    J(1, 0) = dfdx[8*2 + 1];
    J(1, 1) = dfdx[8*2 + 2];
    J(1, 2) = dfdx[8*2 + 5];
    J(1, 3) = dfdx[8*2 + 6];
    J(2, 0) = dfdx[8*5 + 1];
    J(2, 1) = dfdx[8*5 + 2];
    J(2, 2) = dfdx[8*5 + 5];
    J(2, 3) = dfdx[8*5 + 6];
    J(3, 0) = dfdx[8*6 + 1];
    J(3, 1) = dfdx[8*6 + 2];
    J(3, 2) = dfdx[8*6 + 5];
    J(3, 3) = dfdx[8*6 + 6];
    std::cout << state[7] << std::endl;
    std::cout << J << std::endl;
    std::cout << J.eigenvalues() << std::endl << std::endl;
  }
  return 0;
}
