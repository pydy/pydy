% Program to simulate a ladder falling due to gravity. The ladder leans
% against a frictionless wall and is supported by a frictionless floor
function [] = Root_SlidingLadder()
clc; clear all; close all;

% Define Parameters
m = 1; L = 1; g = 10; I = m*L^2/12; % I for a line about CoM

% Pack parameters
params.m = m; params.I = I; params.L = L; params.g = g;

% Initial conditions and time list
StateVar0 = [3*pi/4; 0];

tmax  = 10;
tlist = linspace(0,tmax,10000);

% Call to ODE45
options = odeset('reltol',1e-9,'abstol',1e-9,'Events',@HitFloor);

TargetODE = @(t,statevars) SlidingLadder_ODE(t,statevars,params);
[tlist,statevars] = ode45(TargetODE,tlist,StateVar0,options);

%% Unpack output
thetalist = statevars(:,1);
dthetalist = statevars(:,2);

% Find coordinates of the end-points P1, P2 and the Center of Mass G
x_P1 = zeros(size(thetalist));
y_P1 = L*sin(thetalist);

x_G = -L*cos(thetalist)/2;
y_G = L*sin(thetalist)/2;

x_P2 = -L*cos(thetalist);
y_P2 = zeros(size(thetalist));

%% Calculate Energy 
% Contact forces do no work and should lead to constant energy
vx_G = L*sin(thetalist).*dthetalist/2;
vy_G = L*cos(thetalist).*dthetalist/2;

KE_linear = m*(vx_G.^2 + vy_G.^2)/2;
KE_rotational = I*(dthetalist.^2)/2;

PE_gravity = m*y_G*g;

TE = KE_linear + KE_rotational + PE_gravity;

figure;
plot(tlist,TE);
title('Total Energy vs Time');
xlabel('Time (s)');
ylabel('Total Energy (Joules)');

%% Plot theta with time
figure;
plot(tlist,thetalist*180/pi);
title('Theta vs Time');
xlabel('Time (s)');
ylabel('Theta (Degrees)');

%% Animate Ladder
figure;

% subplot(3,1,1)
point_P1 = plot(x_P1(1),y_P1(1),'ko','markerfacecolor','k');
hold on;
% subplot(3,1,2)
point_G = plot(x_G(1),y_G(1),'ro','markerfacecolor','r');
% subplot(3,1,3)
point_P2 = plot(x_P2(1),y_P2(1),'ko','markerfacecolor','k');

% Draw the ladder, wall and floor
Ladder = line([x_P1(1), x_P2(1)], [y_P1(1), y_P2(1)]);
Wall = line([0 0],[-1.2*L 1.2*L]);
Floor = line([-1.2*L 1.2*L], [0 0]);
   
axis([-1.2*L 1.2*L -1.2*L 1.2*L]);
axis square;

for i = 1:length(tlist)
    set(point_P1,'xdata',x_P1(i),'ydata',y_P1(i));
    set(point_G,'xdata',x_G(i),'ydata',y_G(i));
    set(point_P2,'xdata',x_P2(i),'ydata',y_P2(i));
    set(Ladder,'xdata',[x_P1(i),x_P2(i)],'ydata',[y_P1(i),y_P2(i)]);
    pause(0.01)
end

end

% ODE file for Sliding Ladder
function [dstatevar] = SlidingLadder_ODE(~,statevar,params)

% Unpack params and statevars
m = params.m; g = params.g; I = params.I; L = params.L;

theta = statevar(1);
dtheta = statevar(2);

% Determine derivatives
thetadot = dtheta;

thetadotdot = -(g*L*cos(theta)/2)/(I + (m*(L^2)/4)); 

% Pack derivatives
dstatevar = [thetadot; thetadotdot];

end

% Event function
function [value,isterminal,direction] = HitFloor(~,statevar,~)

    theta = statevar(1); 
    
    value = pi - theta; % angle reaches 180 degrees
    isterminal = 1; % Terminates program
    direction = -1; % pi - theta should be decreasing

end