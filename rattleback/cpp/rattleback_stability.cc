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
  p.Ixy = -0.00000;
  p.Iyz = p.Ixz = 0.0;

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
  Matrix<double, 8, 8> J;

  int N = 11;
  double spin = 10.0;
  for (int k = 0; k < N; ++k) {
    state[7] = -spin + 2*spin/(N - 1) * k;
    rattleback_jacobian(0.0, state, dfdx, dfdt, &p);
    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {
        J(i, j) = dfdx[8*i + j];
      }
    }
    std::cout << state[7] << std::endl;
    std::cout << J << std::endl;
    std::cout << J.eigenvalues() << std::endl << std::endl;
  }
  return 0;
}
