#include <cmath>
#include <gsl/gsl_errno.h>
#include "reproduceKane1982.h"

int rattleback_ode(double t,
                   const double x[],
                   double dxdt[],
                   void *params)
{
  // x = {alpha, beta, gamma, omega_1, omega_2, omega_3}
  double alpha = x[0],
         beta =  x[1],
         gamma = x[2],
         w1 = x[3],
         w2 = x[4],
         w3 = x[5];
  // parameters
  rattleback_params *p = static_cast<rattleback_params *>(params); 
  double a = p->a,
         b = p->b,
         c = p->c,
         h = p->h,
         M = p->M,
         A = p->A,
         B = p->B,
         C = p->C,
         D = p->D,
         g = p->g;
  
  // Kane & Levinson 1982, equations 40, 41, 42, checked DLP
  double alpha_dot = w3*sin(beta) + w1*cos(beta),
         beta_dot  = (-w3*cos(beta) + w1*sin(beta))*tan(alpha) + w2,
         gamma_dot = (w3*cos(beta) - w1*sin(beta))/cos(alpha);

  // Air resistance
  double sigma = 0.0;

  // Equations 37, checked DLP
  double mu1 = -cos(alpha)*sin(beta),
         mu2 = sin(alpha),
         mu3 = cos(alpha)*cos(beta);
  
  // Equations 23, checked DLP
  double mu1_dot = w3*mu2 - w2*mu3,
         mu2_dot = w1*mu3 - w3*mu1,
         mu3_dot = w2*mu1 - w1*mu2;

  // Equations 16, 19, checked DLP
  double epsilon = sqrt((a*mu1)*(a*mu1) + (b*mu2)*(b*mu2) + (c*mu3)*(c*mu3)),
         epsilon_dot = (a*a*mu1*mu1_dot + b*b*mu2*mu2_dot + c*c*mu3*mu3_dot)/epsilon;

  // Equations 17, 18, checked DLP
  double x1 = a*a*mu1/epsilon,
         x2 = b*b*mu2/epsilon,
         x3 = c*c*mu3/epsilon;

  // Equations 26, checked DLP
  double v1 = w2*(h - x3) + w3*x2,
         v2 = -w3*x1 - w1*(h - x3),
         v3 = -w1*x2 + w2*x1;

  // Equations 20, checked DLP
  double x1_dot = a*a*(epsilon*mu1_dot - epsilon_dot*mu1)/(epsilon*epsilon),
         x2_dot = b*b*(epsilon*mu2_dot - epsilon_dot*mu2)/(epsilon*epsilon),
         x3_dot = c*c*(epsilon*mu3_dot - epsilon_dot*mu3)/(epsilon*epsilon);

  // Equations 32, 33, 34, checked DLP
  double zeta1 = w2*(v2 - x3_dot) - w3*(v2 - x2_dot),
         zeta2 = w3*(v1 - x1_dot) - w1*(v3 - x3_dot),
         zeta3 = w1*(v2 - x2_dot) - w2*(v1 - x1_dot);

  // Equations 47, 48, 49, checked DLP
  double F1 = -sigma*w1 + M*g*((x3 - h)*mu2 - x2*mu3),
         F2 = -sigma*w2 + M*g*((h - x3)*mu1 + x1*mu3),
         F3 = -sigma*w3 + M*g*(x2*mu1 - x1*mu2);

  // Equations 53, 54, 55, checked DLP
  double R1 = (D*w1 + (B - C)*w2)*w3,
         R2 = ((C - A)*w1 - D*w2)*w3,
         R3 = D*(w2*w2 - w1*w1) + (A - B)*w1*w2;
  
  // Equations 58 - 63, checked DLP
  double I11 = A + M*(x2*x2 + (h - x3)*(h - x3)),
         I22 = B + M*(x1*x1 + (h - x3)*(h - x3)),
         I33 = C + M*(x1*x1 + x2*x2),
         I12 = D - M*x1*x2,
         I23 = M*(h - x3)*x2,
         I31 = M*(h - x3)*x1;
  double I21 = I12, I32 = I23, I13 = I31;
  
  // Equations 64, 65, 66, checked DLP
  double S1 = M*((h - x3)*zeta2 + x2*zeta3),
         S2 = M*((x3 - h)*zeta1 - x1*zeta3),
         S3 = M*(x1*zeta2 - x2*zeta1);

  // Equation 68, checked DLP
  double Q1 = F1 + R1 + S1,
         Q2 = F2 + R2 + S2,
         Q3 = F3 + R3 + S3;

  // Equation 69
  double G = I11*I22*I33 + I12*I23*I31 + I13*I21*I32 - I13*I22*I31 - I12*I21*I33 - I11*I23*I32;
  
  // Equations 70
  double E1 = Q1*I22*I33 + I12*I23*Q3 + I13*Q2*I32 - I13*I22*Q3 - I12*Q2*I33 - Q1*I23*I32;
  double E2 = I11*Q2*I33 + Q1*I23*I31 + I13*I21*Q3 - I13*Q2*I31 - Q1*I21*I33 - I11*I23*Q3;
  double E3 = I11*I22*Q3 + I12*Q2*I31 + Q1*I21*I32 - Q1*I22*I31 - I12*I21*Q3 - I11*Q2*I32;

  double w1_dot = E1/G,
         w2_dot = E2/G,
         w3_dot = E3/G;

  dxdt[0] = alpha_dot;
  dxdt[1] = beta_dot;
  dxdt[2] = gamma_dot;
  dxdt[3] = w1_dot;
  dxdt[4] = w2_dot;
  dxdt[5] = w3_dot;
  return GSL_SUCCESS;
}

void rattleback_outputs(simdata *s, rattleback_params *p)
{
  double dxdt[6];
  double *x = s->x;
  // x = {alpha, beta, gamma, omega_1, omega_2, omega_3}
  double alpha = x[0],
         beta =  x[1],
         gamma = x[2],
         w1 = x[3],
         w2 = x[4],
         w3 = x[5];
  // parameters
  double a = p->a,
         b = p->b,
         c = p->c,
         h = p->h,
         M = p->M,
         A = p->A,
         B = p->B,
         C = p->C,
         D = p->D,
         g = p->g;
  
  // Kane & Levinson 1982, equations 40, 41, 42, checked DLP
  double alpha_dot = w3*sin(beta) + w1*cos(beta),
         beta_dot  = (-w3*cos(beta) + w1*sin(beta))*tan(alpha) + w2,
         gamma_dot = (w3*cos(beta) - w1*sin(beta))/cos(alpha);

  // Air resistance
  double sigma = 0.0;

  // Equations 37, checked DLP
  double mu1 = -cos(alpha)*sin(beta),
         mu2 = sin(alpha),
         mu3 = cos(alpha)*cos(beta);
  
  // Equations 23, checked DLP
  double mu1_dot = w3*mu2 - w2*mu3,
         mu2_dot = w1*mu3 - w3*mu1,
         mu3_dot = w2*mu1 - w1*mu2;

  // Equations 16, 19, checked DLP
  double epsilon = sqrt((a*mu1)*(a*mu1) + (b*mu2)*(b*mu2) + (c*mu3)*(c*mu3)),
         epsilon_dot = (a*a*mu1*mu1_dot + b*b*mu2*mu2_dot + c*c*mu3*mu3_dot)/epsilon;

  // Equations 17, 18, checked DLP
  double x1 = a*a*mu1/epsilon,
         x2 = b*b*mu2/epsilon,
         x3 = c*c*mu3/epsilon;

  // Equations 26, checked DLP
  double v1 = w2*(h - x3) + w3*x2,
         v2 = -w3*x1 - w1*(h - x3),
         v3 = -w1*x2 + w2*x1;

  // Compute dxdt
  rattleback_ode(s->t, s->x, dxdt, static_cast<void *>(p));

  // Begin copy paste
  s->delta = acos(cos(x[0])*cos(x[1]));

  // End copy paste
  s->alpha[0] = dxdt[3];
  s->alpha[1] = dxdt[4];
  s->alpha[2] = dxdt[5];

  // energy
  s->ke = 0.5*M*(v1*v1 + v2*v2 + v3*v3) + 0.5*(A*w1*w1 + B*w2*w2 + C*w3*w3 + 2*D*w1*w2);
  // s->pe = -g*M*((-pow(b, 2)*sin(x[1])/sqrt(pow(a, 2)*pow(sin(x[2]), 2)*pow(cos(x[1]), 2) + pow(b, 2)*pow(sin(x[1]), 2) + pow(c, 2)*pow(cos(x[1]), 2)*pow(cos(x[2]), 2)) + e)*sin(x[1]) - (pow(a, 2)*sin(x[2])*cos(x[1])/sqrt(pow(a, 2)*pow(sin(x[2]), 2)*pow(cos(x[1]), 2) + pow(b, 2)*pow(sin(x[1]), 2) + pow(c, 2)*pow(cos(x[1]), 2)*pow(cos(x[2]), 2)) + d)*sin(x[2])*cos(x[1]) + (-pow(c, 2)*cos(x[1])*cos(x[2])/sqrt(pow(a, 2)*pow(sin(x[2]), 2)*pow(cos(x[1]), 2) + pow(b, 2)*pow(sin(x[1]), 2) + pow(c, 2)*pow(cos(x[1]), 2)*pow(cos(x[2]), 2)) + f)*cos(x[1])*cos(x[2]));
}
