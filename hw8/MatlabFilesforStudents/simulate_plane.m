%This script simulates a model plane connected to the ceiling by a linear
%spring, as treated in homework assignments of the course "Advanced
%Dynamics" at TUD. 
%Author: H. Vallery, October 2014
%Modified: D. Lemus, Oct 2014

%----------------------------
%define constant parameters:
%----------------------------
Ts=.01;%[s], sampling time
endtime=0.5;%[s] %end time of integration
par.Ts_anim=.01;%[s] pause between frames during animation (can be used for slow motion if larger than Ts)

%gravity:
par.g=9.81;%[m/s^2] acceleration of gravity (points in positive z direction)

%mass properties:
Ip=1*.1^2;%[kgm^2] moment of inertia of the plane about the z axis
par.Ixx=.5*Ip;
par.Iyy=.5*Ip;
par.Izz=Ip;

par.m=1;%[kg], mass of the plane

%geometry:
px=.1;%[m] local x coordinate of point P (spring attachment point) in B frame
pz=-.1;%[m] local z coordinate of P
par.rp_B=[px;0;pz];%coordinates of point P in B frame

par.l0=1;%[m], resting length of spring
par.k=20;%[N/m], stiffness of spring

%additional geometric data for aninmation of the plane:
par.length_plane=2;%[m], length of the plane's body
par.wingspan=4;%[m] wingspan
par.height_fin=.5;%[m] height of the fin
par.r_wingcenter=[.1;0;-.01];%[m] vector from CoM to wing center in body frame


%----------------------------
%set initial conditions:
%----------------------------

%Cartesian positions of the plane's center of mass:
sX=-px;
sY=.2;%.2;
sZ=par.l0+par.m*par.g/par.k;
%corresponding velocities:
dsX=2;
dsY=0;
dsZ=0;
%Euler angles (type I) to describe the orientation of the plane with
%respect to the inertial N frame (XYZ):
psi=0*pi/3;%rotation about Z axis
theta=0*pi/20;%rotation about intermediate Y' axis
phi=0*pi/10;%rotation about new x axis
%angular velocities (expressed in the body-fixed B frame (xyz) of the plane):
omegax=0;
omegay=0;
omegaz=pi/4/10;

%----------------------------
%integrate:
%----------------------------
x0=[sX,sY,sZ,psi,theta,phi,dsX,dsY,dsZ,omegax,omegay,omegaz];%vector of initial conditions:
options = odeset('AbsTol',1e-10,'RelTol',1e-8);
%[t,y]=ode45(@plane_equationsofmotion,[0:Ts:endtime],x0,options,par); 
[t,y0] = Integrate_Euler(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
%[t,y] = Integrate_ModifiedEuler(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
%[t,y] = Integrate_RungeKutta2(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
[t,y1] = Integrate_RungeKutta4(@plane_equationsofmotion,[0,endtime],Ts, x0,par);

%extract the generalized coordinates as time series:
qmatrix0=y0(:,1:6);
qmatrix1=y1(:,1:6);

%----------------------------
%animate the result:
%----------------------------
planeanimation(qmatrix0,par);
planeanimation(qmatrix1,par);

% Define Rotation Matrix for Euler angles type I
R_psi = @(psi) [cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0; 0 0 1];
R_theta= @(theta) [cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
R_phi= @(phi) [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];


x0 = y0(end,1:3);
psi = y0(end,4);
theta = y0(end,5);
phi = y0(end,6);

% Rotation Matrix final orientation euler's integration method
R_N2E= R_phi(phi)*R_theta(theta)*R_psi(psi);

x1 = y1(end,1:3);
psi = y1(end,4);
theta = y1(end,5);
phi = y1(end,6);

% Rotation Matrix final orientation RK4 integration method
R_N2RK= R_phi(phi)*R_theta(theta)*R_psi(psi);

% Rotation Matrix from Euler to RK4 final orientations
R = R_N2RK*R_N2E';

% principal angle of rotation 
phi = acos((trace(R)-1)*0.5);

% Distance between final center of mass postions Euler and RK4 integration
% methods
dist = norm(x1-x0);

disp(sprintf('instantaneous angle of rotation = %g',phi))
disp(sprintf('final disance = %g',dist))