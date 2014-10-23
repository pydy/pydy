function dx=plane_equationsofmotion(t,x,par)
%function dx=plane_equationsofmotion(t,x,par)
%This function contains the equations of motion of a model airplane
%connected to the ceiling by a linear spring, as used in homework
%assignments of the course "Advanced Dynamics" at TUD.
%
%Input: t (time) and 
%x: Cartesian positions, Euler angles, linear speeds, omega
%Output: dx: derivatives thereof
%
%H. Vallery, October 2014

%----------------------------
%extract parameters:
%----------------------------

g=par.g;%[m/s^2], acceleration of gravity (points in positive z direction)

%spring parameters:
l0=par.l0;%m, resting length of spring
k=par.k;%N/m, stiffness of spring
rp_B=par.rp_B;%position vector of point P (spring attachment point) in B frame

%mass parameters of the plane:
m=par.m;% mass of plane
Ixx=par.Ixx;% moment of inertia about x axis
Iyy=par.Iyy;% moment of inertia about y axis
Izz=par.Izz;% moment of inertia about z axis

%----------------------------
%read current states:
%----------------------------

sX=x(1);
sY=x(2);
sZ=x(3);
psi=x(4);
theta=x(5);
phi=x(6);

dsX=x(7);
dsY=x(8);
dsZ=x(9);
omegax=x(10);
omegay=x(11);
omegaz=x(12);

rs_N=[sX;sY;sZ];%position of point S (center of mass) in N coordinates

%----------------------------
%calculate forces and moments in body frame:
%----------------------------

%construction of rotation matrix that maps vectors that are
%expressed in the N frame to their representation in the B frame:
R_psi = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];
R_theta = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
R_phi = [1, 0, 0; 0, cos(phi), sin(phi); 0, -sin(phi), cos(phi)];
R_total = R_phi * R_theta * R_psi;

%position vector of point P in the N frame:
rp_N = rs_N + R_total'*rp_B;

%external forces acting on the plane:
l = sqrt(rp_N'*rp_N);
Fspring_N = -k*(l - l0) * rp_N/l; %spring force in the N frame
Fg_N = [0; 0; m*g]; %gravitational force in the N frame
Ftot_N = (Fspring_N + Fg_N); %total force in the N frame

%map the spring force to the body frame:
Fspring_B = R_total*Fspring_N; %spring force in the body frame

%calculate the external moments acting about the plane's center of mass:
M = cross(rp_B, Fspring_B); %moment of the spring in the body frame
%components about the body-fixed axes:
Mx = M(1); My = M(2); Mz = M(3);

%----------------------------
%apply Newton-Euler:
%----------------------------

%Newton equations (in space-fixed frame):
ddsX = Ftot_N(1)/m;
ddsY = Ftot_N(2)/m;
ddsZ = Ftot_N(3)/m;

%Euler equations (in body-fixed frame):
domegax = (Mx - (Izz - Iyy)*omegay*omegaz)/Ixx;
domegay = (My - (Ixx - Izz)*omegaz*omegax)/Iyy;
domegaz = (Mz - (Iyy - Ixx)*omegax*omegay)/Izz;

%----------------------------
%calculate state derivatives:
%----------------------------

%Euler angle derivatives, Greenwood (3.16)
dpsi=sec(theta)* (omegay*sin(phi)+omegaz*cos(phi));
dtheta = omegay*cos(phi)-omegaz*sin(phi);
dphi = omegax + dpsi *sin(theta);

dx=[dsX;dsY;dsZ;dpsi;dtheta;dphi;ddsX;ddsY;ddsZ;domegax;domegay;domegaz];
