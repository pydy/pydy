%This script simulates a model plane connected to the ceiling by a linear
%spring, as treated in homework assignments of the course "Advanced
%Dynamics" at TUD. 
%Author: H. Vallery, October 2014


%----------------------------
%define constant parameters:
%----------------------------
Ts=.01;%[s], sampling time
endtime=2;%[s] %end time of integration
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
%if strcmp(solver, 'ode45')
%    [t,y]=ode45(@plane_equationsofmotion,[0:Ts:endtime],x0,options,par); 
%elseif strcmp(solver, 'euler')
%    [t,y] = Integrate_EulersMethod(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
%elseif strcmp(solver, 'modifiedeuler')
%    [t,y] = Integrate_ModifiedEulerMethod(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
%elseif strcmp(solver, 'rk2')
%    [t,y] = Integrate_RungeKutta2(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
%elseif strcmp(solver, 'rk4')
%    [t,y] = Integrate_RungeKutta4(@plane_equationsofmotion,[0,endtime],Ts, x0,par);
%else
%    error('invalid solver')
%end

[t1, y1] = Integrate_EulersMethod(@plane_equationsofmotion,...
    [0, endtime], Ts, x0, par);
[t2, y2] = Integrate_RungeKutta4(@plane_equationsofmotion,...
    [0, endtime], Ts, x0, par);

R_z = @(x) [cos(x), sin(x), 0; -sin(x), cos(x), 0; 0, 0, 1];
R_y = @(x) [cos(x), 0, -sin(x); 0, 1, 0; sin(x), 0, cos(x)];
R_x = @(x) [1, 0, 0; 0, cos(x), sin(x); 0, -sin(x), cos(x)];
R_zyx = @(x) R_x(x(3)) * R_y(x(2)) * R_z(x(1));

dx = y1(end, 1:3) - y2(end, 1:3);
dist = sqrt(dx*dx');
disp(sprintf('distance between final positions = %g', dist))

R1 = R_zyx(y1(end, 4:6));
R2 = R_zyx(y2(end, 4:6));

R_21 = R2*R1';
angle_rot = acos((trace(R_21) - 1)/2);
disp(sprintf('angle of rotation between final orientations = %g', angle_rot))

%extract the generalized coordinates as time series:
%qmatrix=y1(:,1:6);

%----------------------------
%animate the result:
%----------------------------
%planeanimation(qmatrix,par);
