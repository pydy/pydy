%This script simulates a three-wheeled cart, 
%as treated in week 7 of the course "Advanced
%Dynamics" at TUD. 
%Author: H. Vallery, October 2014


%----------------------------
%define constant parameters:
%----------------------------
Ts=.04;%[s], sampling time
endtime=10;%[s] %end time of integration
par.Ts_anim=.04;%[s] pause between frames during animation 
%(can be used for slow motion if larger than Ts)


%geometry:

par.length_cart=2;%[m], length of the cart
par.width_cart=1;%[m] width of the cart

%mass properties:
par.m=20;%[kg], mass of the cart
par.Is=par.m*1/12*(par.length_cart^2+par.width_cart^2);%[kgm^2] moment of inertia 
%of the cart about the z axis

%----------------------------
%set initial conditions:
%----------------------------

%Cartesian positions of the cart's center of mass:
sX=0;%[m]
sY=.2;%[m]
%corresponding velocities:
dsX=2;%[m/s]
dsY=2;%[m/s]

%orientation of the cart with
%respect to the inertial N frame (XYZ):
theta=0;%[rad], rotation about z axis
%angular velocity: 
omega=2;%[rad/s]
%Remark: These initial conditions fulfill the given constraint

%----------------------------
%integrate:
%----------------------------
x0=[sX,sY,theta,dsX,dsY,omega];%vector of initial conditions:
options = odeset('AbsTol',1e-10,'RelTol',1e-8);
%[t,y]=ode45(@cart_equationsofmotion,[0:Ts:endtime],x0,options,par); 
%[t,y] = Integrate_Euler(@cart_equationsofmotion,[0,endtime],Ts, x0,par);
[t,y] = Integrate_ModifiedEuler(@cart_equationsofmotion,[0,endtime],Ts, x0,par);
%[t,y] = Integrate_RungeKutta2(@cart_equationsofmotion,[0,endtime],Ts, x0,par);
%[t,y] = Integrate_RungeKutta4(@cart_equationsofmotion,[0,endtime],Ts, x0,par);


%extract the generalized coordinates as time series:
qmatrix=y(:,1:3);

%----------------------------
%animate the result:
%----------------------------

cartanimation(qmatrix,par);