%This script simulates a model cart connected to the ceiling by a linear
%spring, as treated in homework assignments of the course "Advanced
%Dynamics" at TUD. 
%Author: H. Vallery, October 2014


%----------------------------
%define constant parameters:
%----------------------------
Ts=.01;%[s], sampling time
endtime=5;%[s] %end time of integration
par.Ts_anim=.04;%[s] pause between frames during animation (can be used for slow motion if larger than Ts)


%geometry:

par.length_cart=2;%[m], length of the cart
par.width_cart=1;%[m] width of the cart

%mass properties:
par.m=20;%[kg], mass of the cart
par.Is=par.m*1/12*(par.length_cart^2+par.width_cart^2);%[kgm^2] moment of inertia of the cart about the z axis

% calculate forces
dt = 0.01;
N = 500;
t = (0:(N - 1))*dt;

f = 2;
v0 = 20;
l = par.length_cart;

% constraint:
% dsX*sin(theta) - dsY*cos(theta) + l/2*omega = 0

theta_t = pi/8 * sin(2*pi*f*t);
omega_t = pi/8 * 2*pi*f * cos(2*pi*f*t);
domega_t = -(2*pi*f)^2 * theta_t;
dsX_t = -v0*ones(size(t));
ddsX_t = zeros(size(t));

dsY_t = l/2*omega_t + dsX_t.*sin(theta_t);
ddsY_t = -1./cos(theta_t) .* (-dsX_t.*omega_t.*cos(theta_t) -...
    dsY_t.*omega_t.*sin(theta_t) - l/2.*domega_t - sin(theta_t).*ddsX_t);

[Fx, Fy] = inversedynamics_cart(theta_t, ddsX_t, ddsY_t, domega_t, par);
disp(sprintf('Fx(end) = %g', Fx(end)))
disp(sprintf('Fy(end) = %g', Fy(end)))

F = [Fx, Fy];
par.F = F;


%----------------------------
%set initial conditions:
%----------------------------

%Cartesian positions of the cart's center of mass:
sX=0;
sY=.2;
%corresponding velocities:
dsX = -v0;
dsY = pi^2*2/8;
%orientation of the cart with
%respect to the inertial N frame (XYZ):
theta = 0;
%angular velocity 
omega = pi/8 * 2*pi*f;

%----------------------------
%integrate:
%----------------------------
x0=[sX,sY,theta,dsX,dsY,omega];%vector of initial conditions:
options = odeset('AbsTol',1e-10,'RelTol',1e-8);
%[t,y]=ode45(@cart_equationsofmotion,[0:Ts:endtime],x0,options,par); 
%[t,y] = Integrate_EulersMethod(@cart_equationsofmotion,[0,endtime],Ts, x0,par);
[t,y] = Integrate_ModifiedEulerMethod(@cart_equationsofmotion,[0,endtime-Ts],Ts, x0,par);
%[t,y] = Integrate_RungeKutta2(@cart_equationsofmotion,[0,endtime],Ts, x0,par);
%[t,y] = Integrate_RungeKutta4(@cart_equationsofmotion,[0,endtime],Ts, x0,par);


%extract the generalized coordinates as time series:
qmatrix=y(:,1:3);

x_sim = y(end, 1);
x_mdl = -v0*(N-1)*dt;
disp(sprintf('error in sX: %g - %g = %g', x_sim, x_mdl, x_sim - x_mdl))

% calculate constraint force:
% given the Lagrange equations:
%  m*ddsX = Fx*cos(theta) - Fy*sin(theta) + sin(theta)*lambda
%  m*ddsY = Fx*sin(theta) + Fy*cos(theta) - cos(theta)*lambda
%  Is*ddtheta = Fy*l/2*lambda

dy = zeros(size(y));
for i = 1:length(y) - 1
    dy(i, :) = cart_equationsofmotion(t(i), y(i, :), par);
end
Cx = par.m*dy(:, 4) - par.F(:, 1).*cos(y(:, 3)) + par.F(:, 2).*sin(y(:, 3));
Cy = -1*(par.m*dy(:, 5) - par.F(:, 1).*sin(y(:, 3)) -...
         par.F(:, 2).*cos(y(:, 3)));

C_force = sqrt(Cx.^2 + Cy.^2);
F_fr = 0.7*par.m/2*9.81;
disp(sprintf('static friction force = μ*m/2*g = %g', F_fr));
figure(); plot(t, C_force, 'b', [t(1), t(end)], F_fr*[1, 1], 'r-'); pause

%----------------------------
%animate the result:
%----------------------------

cartanimation(qmatrix,par);
