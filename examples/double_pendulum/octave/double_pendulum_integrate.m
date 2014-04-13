function double_pendulum_integrate()
% function double_pendulum_integrate()
% Simulates the double pendulum.

% Integrate the equations over from 0 to 5 seconds with 50 steps.
timeSpan = linspace(0.0, 5.0, 50)';

% Set initial angles in radians and the initial speeds to zero.
initialConditions = [0.1; 0.2; 0.0; 0.0];

% Define particles' mass, pendulums' length, and the acceleration due to
% gravity.
m = 1.0;
l = 1.0;
g = 9.8;

% Create a function handle to an anonymous function of which we can pass the
% parameters.
f = @(x, t) double_pendulum_state_derivatives(x, t, m, l, g);

% Integrate the equations of motion with default integration settings.
x = lsode(f, initialConditions, timeSpan);

% Plot the results.
fig = figure(1);
plot(timeSpan, x)
title('Double pendulum states as a function of time')
xlabel('Time [s]')
legend('q1 [rad]', 'q2 [rad]', 'u1 [rad/s]', 'u2 [rad/s]')

% Save the figure as a png.
print -dpng double-pendulum-octave-plot.png
