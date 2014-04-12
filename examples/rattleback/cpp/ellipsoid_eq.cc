#include <complex>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>

#include "ellipsoid_no_slip.h"

using namespace Eigen;
using namespace std;

int main()
{
  using namespace Mechanics::Rattleback::Ellipsoid::NoSlip;

  parameters p;
  p.a = 0.2;                // x semi-diameter
  p.b = 0.03;               // y semi-diameter
  p.c = 0.02;               // z semi-diameter 
  p.d = p.e = 0.0;          // COM offset in x, y directions
  p.f = 0.01;               // COM offset in z direction
  p.m = 1.0;                // mass
  p.g = 9.81;               // graviational constant
  p.Ixx =  0.0002;          // xx moment of inertia
  p.Iyy =  0.0016;          // yy moment of inertia
  p.Izz =  0.0017;          // zz moment of inertia
  p.Ixy = -0.00002;         // xy product of inertia
  p.Iyz = p.Ixz = 0.0;      // yz, xz products of inertia
  p.s = 0.001;             // sigma, viscous air damping

  double x[8] = {0.0,       // Yaw (ignorable)
                 0.0,       // Roll
                 0.0,       // Pitch
                 0.0, 0.0,  // x, y of contact (ignorable)
                 0.0,       // u0
                 0.0,       // u1
                 0.0};      // u2  (spin)
  double dfdx[64], dfdt[8];

  Matrix<double, 4, 4> A;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      A(i, j) = 0.0;
  EigenSolver<Matrix<double, 4, 4> > evalSolver;
  Matrix<complex<double>, 4, 1> eigs;
  
  // Open a file for writing
  std::ofstream f("eigenvalues.dat", std::ios::binary | std::ios::out);

  int N = 5001;
  double minSpin = -5.0;
  double maxSpin = 50.0;
  double delta = (maxSpin - minSpin)/(N - 1);

  for (int i = 0; i < N; ++i) {
     x[7] = minSpin + i*delta;
     jacobian(0.0, x, dfdx, dfdt, &p);

     A(0, 0) = dfdx[1*8 + 1];
     A(0, 1) = dfdx[1*8 + 2];
     A(0, 2) = dfdx[1*8 + 5];
     A(0, 3) = dfdx[1*8 + 6];
     A(1, 0) = dfdx[2*8 + 1];
     A(1, 1) = dfdx[2*8 + 2];
     A(1, 2) = dfdx[2*8 + 5];
     A(1, 3) = dfdx[2*8 + 6];
     A(2, 0) = dfdx[5*8 + 1];
     A(2, 1) = dfdx[5*8 + 2];
     A(2, 2) = dfdx[5*8 + 5];
     A(2, 3) = dfdx[5*8 + 6];
     A(3, 0) = dfdx[6*8 + 1];
     A(3, 1) = dfdx[6*8 + 2];
     A(3, 2) = dfdx[6*8 + 5];
     A(3, 3) = dfdx[6*8 + 6];

     evalSolver.compute(A, false);
     eigs = evalSolver.eigenvalues();

     // std::cout << eigs << std::endl;
     f.write((char *) &(x[7]), sizeof(double));
     for (int j = 0; j < 4; ++j) {
       double re_im[2] = {eigs(j).real(), eigs(j).imag()};
       f.write((char *) re_im, 2*sizeof(double));
     }
  }

  f.close();
  return 0;
}
