#include <cmath>
#include <gsl/gsl_errno.h>
#include "rattleback.h"

int rattleback_ode(double t,
                   const double x[],
                   double dxdt[],
                   void *params)
{
  rattleback_params *p = static_cast<rattleback_params *>(params);
  const double a = p->a,
               b = p->b,
               c = p->c,
               d = p->d,
               e = p->e,
               f = p->f,
               m = p->m,
               Ixx = p->Ixx,
               Iyy = p->Iyy,
               Izz = p->Izz,
               Ixy = p->Ixy,
               Iyz = p->Iyz,
               Ixz = p->Ixz,
               g = p->g;

  double mu0 = -sin(x[2])*cos(x[1]),
         mu1 = sin(x[1]),
         mu2 = cos(x[1])*cos(x[2]);

  // Equations 23, checked DLP
  double mu0_dot = x[7]*mu1 - x[6]*mu2,
         mu1_dot = x[5]*mu2 - x[7]*mu0,
         mu2_dot = x[6]*mu0 - x[5]*mu1;
  double epsilon = sqrt((a*mu0)*(a*mu0) + (b*mu1)*(b*mu1) + (c*mu2)*(c*mu2)),
         epsilon_dot = (a*a*mu0*mu0_dot + b*b*mu1*mu1_dot + c*c*mu2*mu2_dot)/epsilon;
  // Equations 17, 18, checked DLP
  double r0 = a*a*mu0/epsilon,
         r1 = b*b*mu1/epsilon,
         r2 = c*c*mu2/epsilon;
  // Equations 20, checked DLP
  double rd0 = a*a*(epsilon*mu0_dot - epsilon_dot*mu0)/(epsilon*epsilon),
         rd1 = b*b*(epsilon*mu1_dot - epsilon_dot*mu1)/(epsilon*epsilon),
         rd2 = c*c*(epsilon*mu2_dot - epsilon_dot*mu2)/(epsilon*epsilon);

  // Begin Copy Paste
  dxdt[0] = (-x[5]*sin(x[2]) + x[7]*cos(x[2]))/cos(x[1]);
  dxdt[1] = x[5]*cos(x[2]) + x[7]*sin(x[2]);
  dxdt[2] = (x[5]*sin(x[2]) - x[7]*cos(x[2]))*tan(x[1]) + x[6];
  dxdt[3] = rd0*(-sin(x[0])*sin(x[1])*sin(x[2]) + cos(x[0])*cos(x[2])) - rd1*sin(x[0])*cos(x[1]) + rd2*(sin(x[0])*sin(x[1])*cos(x[2]) + sin(x[2])*cos(x[0]));
  dxdt[4] = rd0*(sin(x[0])*cos(x[2]) + sin(x[1])*sin(x[2])*cos(x[0])) + rd1*cos(x[0])*cos(x[1]) + rd2*(sin(x[0])*sin(x[2]) - sin(x[1])*cos(x[0])*cos(x[2]));


  dxdt[5] = (-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((-d + r0)*cos(x[1])*cos(x[2]) - (f - r2)*sin(x[2])*cos(x[1])) + m*(-d + r0)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(f - r2)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - (-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))*(-g*m*((d - r0)*sin(x[1]) - (-e + r1)*sin(x[2])*cos(x[1])) + m*(d - r0)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) + m*(-e + r1)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixz - m*(d - r0)*(-f + r2))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - (-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))*(-g*m*((-d + r0)*cos(x[1])*cos(x[2]) - (f - r2)*sin(x[2])*cos(x[1])) + m*(-d + r0)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(f - r2)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) + (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[7] - (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[5])/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2))) - (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[6] + (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[5])/(-Izz - m*pow(d - r0, 2) - m*pow(-e + r1, 2) - pow(-Ixz - m*(d - r0)*(-f + r2), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - pow(-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)), 2)/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))) + (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[7] - (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[5])/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2))) - (-Ixz - m*(d - r0)*(-f + r2))*(-g*m*((d - r0)*sin(x[1]) - (-e + r1)*sin(x[2])*cos(x[1])) + m*(d - r0)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) + m*(-e + r1)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixz - m*(d - r0)*(-f + r2))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - (-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))*(-g*m*((-d + r0)*cos(x[1])*cos(x[2]) - (f - r2)*sin(x[2])*cos(x[1])) + m*(-d + r0)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(f - r2)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) + (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[7] - (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[5])/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2))) - (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[6] + (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[5])/(-Izz - m*pow(d - r0, 2) - m*pow(-e + r1, 2) - pow(-Ixz - m*(d - r0)*(-f + r2), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - pow(-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)), 2)/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2));
  dxdt[6] = (-g*m*((-d + r0)*cos(x[1])*cos(x[2]) - (f - r2)*sin(x[2])*cos(x[1])) + m*(-d + r0)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(f - r2)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - (-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))*(-g*m*((d - r0)*sin(x[1]) - (-e + r1)*sin(x[2])*cos(x[1])) + m*(d - r0)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) + m*(-e + r1)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixz - m*(d - r0)*(-f + r2))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - (-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))*(-g*m*((-d + r0)*cos(x[1])*cos(x[2]) - (f - r2)*sin(x[2])*cos(x[1])) + m*(-d + r0)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(f - r2)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) + (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[7] - (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[5])/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2))) - (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[6] + (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[5])/(-Izz - m*pow(d - r0, 2) - m*pow(-e + r1, 2) - pow(-Ixz - m*(d - r0)*(-f + r2), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - pow(-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)), 2)/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))) + (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[7] - (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[5])/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)));
  dxdt[7] = (-g*m*((d - r0)*sin(x[1]) - (-e + r1)*sin(x[2])*cos(x[1])) + m*(d - r0)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) + m*(-e + r1)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixz - m*(d - r0)*(-f + r2))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - (-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)))*(-g*m*((-d + r0)*cos(x[1])*cos(x[2]) - (f - r2)*sin(x[2])*cos(x[1])) + m*(-d + r0)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(f - r2)*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7]) - (-Ixy - m*(-d + r0)*(e - r1))*(-g*m*((e - r1)*cos(x[1])*cos(x[2]) + (-f + r2)*sin(x[1])) + m*(e - r1)*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6]) + m*(-f + r2)*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7]) - (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[7] + (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[6])/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) + (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[7] - (Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[5])/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2))) - (Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[6] + (Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[5])/(-Izz - m*pow(d - r0, 2) - m*pow(-e + r1, 2) - pow(-Ixz - m*(d - r0)*(-f + r2), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)) - pow(-Iyz - m*(-e + r1)*(f - r2) - (-Ixy - m*(-d + r0)*(e - r1))*(-Ixz - m*(d - r0)*(-f + r2))/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2)), 2)/(-Iyy - m*pow(-d + r0, 2) - m*pow(f - r2, 2) - pow(-Ixy - m*(-d + r0)*(e - r1), 2)/(-Ixx - m*pow(e - r1, 2) - m*pow(-f + r2, 2))));
  // End Copy Paste

  return GSL_SUCCESS;
}

void rattleback_outputs(simdata *s, rattleback_params *p)
{
  double dxdt[8];
  double *x = s->x;
  const double a = p->a,
               b = p->b,
               c = p->c,
               d = p->d,
               e = p->e,
               f = p->f,
               m = p->m,
               Ixx = p->Ixx,
               Iyy = p->Iyy,
               Izz = p->Izz,
               Ixy = p->Ixy,
               Iyz = p->Iyz,
               Ixz = p->Ixz,
               g = p->g;

  double mu0 = -sin(x[2])*cos(x[1]),
         mu1 = sin(x[1]),
         mu2 = cos(x[1])*cos(x[2]);
  // Equations 23, checked DLP
  double mu0_dot = x[7]*mu1 - x[6]*mu2,
         mu1_dot = x[5]*mu2 - x[7]*mu0,
         mu2_dot = x[6]*mu0 - x[5]*mu1;
  double epsilon = sqrt((a*mu0)*(a*mu0) + (b*mu1)*(b*mu1) + (c*mu2)*(c*mu2)),
         epsilon_dot = (a*a*mu0*mu0_dot + b*b*mu1*mu1_dot + c*c*mu2*mu2_dot)/epsilon;
  // Equations 17, 18, checked DLP
  double r0 = a*a*mu0/epsilon,
         r1 = b*b*mu1/epsilon,
         r2 = c*c*mu2/epsilon;
  // Equations 20, checked DLP
  double rd0 = a*a*(epsilon*mu0_dot - epsilon_dot*mu0)/(epsilon*epsilon),
         rd1 = b*b*(epsilon*mu1_dot - epsilon_dot*mu1)/(epsilon*epsilon),
         rd2 = c*c*(epsilon*mu2_dot - epsilon_dot*mu2)/(epsilon*epsilon);
  // Compute dxdt
  rattleback_ode(s->t, s->x, dxdt, static_cast<void *>(p));

  // Begin copy paste
  s->CF[0] = m*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6])*sin(x[2]) + m*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7])*cos(x[2]) - (-m*dxdt[5]*(e - r1) - m*dxdt[6]*(-d + r0))*sin(x[2]) - (-m*dxdt[6]*(f - r2) - m*dxdt[7]*(-e + r1))*cos(x[2]);
  s->CF[1] = -m*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6])*sin(x[1])*cos(x[2]) + m*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7])*cos(x[1]) + m*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7])*sin(x[1])*sin(x[2]) + (-m*dxdt[5]*(e - r1) - m*dxdt[6]*(-d + r0))*sin(x[1])*cos(x[2]) - (-m*dxdt[5]*(-f + r2) - m*dxdt[7]*(d - r0))*cos(x[1]) - (-m*dxdt[6]*(f - r2) - m*dxdt[7]*(-e + r1))*sin(x[1])*sin(x[2]);
  s->CF[2] = -g*m + m*(rd0*x[6] - rd1*x[5] + ((d - r0)*x[7] - (f - r2)*x[5])*x[5] - (-(e - r1)*x[7] + (f - r2)*x[6])*x[6])*cos(x[1])*cos(x[2]) + m*(-rd0*x[7] + rd2*x[5] - (-(d - r0)*x[6] + (e - r1)*x[5])*x[5] + (-(e - r1)*x[7] + (f - r2)*x[6])*x[7])*sin(x[1]) - m*(rd1*x[7] - rd2*x[6] + (-(d - r0)*x[6] + (e - r1)*x[5])*x[6] - ((d - r0)*x[7] - (f - r2)*x[5])*x[7])*sin(x[2])*cos(x[1]) - (-m*dxdt[5]*(e - r1) - m*dxdt[6]*(-d + r0))*cos(x[1])*cos(x[2]) - (-m*dxdt[5]*(-f + r2) - m*dxdt[7]*(d - r0))*sin(x[1]) + (-m*dxdt[6]*(f - r2) - m*dxdt[7]*(-e + r1))*sin(x[2])*cos(x[1]);
  s->ke = 0.5*m*(pow(-(d - r0)*x[6] + (e - r1)*x[5], 2) + pow((d - r0)*x[7] - (f - r2)*x[5], 2) + pow(-(e - r1)*x[7] + (f - r2)*x[6], 2)) + 0.5*(Ixx*x[5] + Ixy*x[6] + Ixz*x[7])*x[5] + 0.5*(Ixy*x[5] + Iyy*x[6] + Iyz*x[7])*x[6] + 0.5*(Ixz*x[5] + Iyz*x[6] + Izz*x[7])*x[7];
  s->pe = -g*m*(-(d - r0)*sin(x[2])*cos(x[1]) + (e - r1)*sin(x[1]) + (f - r2)*cos(x[1])*cos(x[2]));
  s->te = s->ke + s->pe;
  s->delta = acos(cos(x[1])*cos(x[2]));
  // End copy paste

  s->alpha[0] = dxdt[5];
  s->alpha[1] = dxdt[6];
  s->alpha[2] = dxdt[7];
}
