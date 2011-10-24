function chaos_pendulum_integrate()
% function chaos_pendulum_integrate()
% Integrates the equations of motion of the chaos pendulum and plots the
% results. Edit the values in this function to change the time histories.

% Define the system's parameters in a structure for easy passing.
p.g = 9.81; % 9.81 m/s^2
p.lA = 0.075; % m
p.mA = 0.01; % kg
p.IAxx = 5E-6; % kg * m^2
p.lB = 0.2; % m
p.mB = 0.1; % kg
p.IBxx = 50E-6; % kg * m^2
p.IByy = 250E-6; % kg * m^2
p.IBzz = 200e-6; % kg * m^2

% Create a function handle to an anonymous function which which can pass the
% parameters.
f = @(t, x) state_derivatives(t, x, p);

% Integrate the equations from 0 to 10 seconds at a 0.02 second time step.
timeSpan = linspace(0.0, 10.0, 500);

% Stable %
theta0 = deg2rad(45);
% Set initial angles in radians and the initial speeds to zero.
initialConditions = [theta0, deg2rad(0.5), 0.0, 0.0];
% Integrate the equations of motion with default integration settings.
[t, x] = ode45(f, timeSpan, initialConditions);
% Plot the results for the first integration.
fig1 = figure(1);
plot(t, rad2deg(x(:, 2)))
% Now change the initial condition of phi by half a degree
initialConditions = [theta0, deg2rad(1.0), 0.0, 0.0];
% Integrate the equations of motion with the specified parameters.
[t, x] = ode45(f, timeSpan, initialConditions);
hold all
plot(t, rad2deg(x(:, 2)))
title('Plate angle as a function of time, \theta_o = 45 degrees')
xlabel('Time [s]')
ylabel('\phi [deg]')
legend('\phi_o = 0.5 degrees', '\phi_o = 1.0 degrees')
% Save the figure as a png.
saveas(fig1, 'chaos-pendulum-stable-matlab-plot.png', 'png')

% Chaotic %
theta0 = deg2rad(90.0);
% Set initial angles in radians and the initial speeds to zero.
initialConditions = [theta0, deg2rad(0.5), 0.0, 0.0];
% Integrate the equations of motion with default integration settings.
[t, x] = ode45(f, timeSpan, initialConditions);
% Plot the results for the first integration.
fig2 = figure(2);
plot(t, rad2deg(x(:, 2)))
% Now change the initial condition of phi by half a degree
initialConditions = [theta0, deg2rad(1.0), 0.0, 0.0];
% Integrate the equations of motion with the specified parameters.
[t, x] = ode45(f, timeSpan, initialConditions);
hold all
plot(t, rad2deg(x(:, 2)))
title('Plate angle as a function of time, \theta_o = 90 degrees')
xlabel('Time [s]')
ylabel('\phi [deg]')
legend('\phi_0 = 0.5 degrees', '\phi_0 = 1.0 degrees')
% Save the figure as a png.
saveas(fig2, 'chaos-pendulum-chaotic-matlab-plot.png', 'png')

function xd = state_derivatives(t, x, p)
%function xd = state_derivatives(t, x, p)
% Returns the derivatives of the states as a function of the current state
% and time.
%
% Parameters
% ----------
% t : double
%   Current time.
% x : vector, (4, 1)
%   Current state [theta, phi, omega, alpha].
% p : structure
%   Contains all of the model parameters.
%
% Returns
% -------
% xd : matrix, 4 x 1
%   The derivative of the current state.

% Unpack the variables so that you can use the sympy equations as is.
theta = x(1);
phi = x(2);
omega = x(3);
alpha = x(4);

% Initialize a vector for the derivatives.
xd = zeros(4, 1);
% Calculate the derivatives of the states. These equations can be copied
% directly from the sympy output but be sure to print with `mprint` in
% sympy.physics.mechanics to remove the `(t)` and use Matlab's find and
% replace to change the python `**` to the matlab `^`. Also note that the
% structure `p` was used to pass in the parameters and each parameter must
% be prepended with `p.`.

% theta'
xd(1) = omega;
% phi'
xd(2) = alpha;
% omega'
xd(3) = (-2 * p.IBxx * alpha * omega * sin(phi) * cos(phi) + 2 * ...
    p.IByy * alpha * omega * sin(phi) * cos(phi) - p.g * p.lA * p.mA * ...
    sin(theta) - p.g * p.lB * p.mB * sin(theta)) / (p.IAxx + p.IBxx * ...
    sin(phi)^2 + p.IByy * cos(phi)^2 + p.lA^2 * p.mA + p.lB^2 * p.mB);
% alpha'
xd(4) = (p.IBxx - p.IByy) * omega^2 * sin(phi) * cos(phi) / p.IBzz;
