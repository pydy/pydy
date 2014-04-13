function xd = double_pendulum_state_derivatives(x, t, m, l, g)
%function xd = double_pendulum_state_derivatives(x, t, m, l, g)
% Returns the derivatives of the states as a function of the current state
% and time.
%
% Parameters
% ----------
% x : vector, 4 x 1
%   Current state.
% t : double
%   Current time.
% m : double
%   The mass of each particle.
% l : double
%   Length of each link.
% g : double
%   Acceleratoin due to gravity.
%
% Returns
% -------
% xd : vector, 4 x 1
%   The derivatives of the states in this order [q1, q2, u1, u2].

% Unpack the variables so that you can use the sympy equations as is.
q1 = x(1);
q2 = x(2);
u1 = x(3);
u2 = x(4);

% Initialize a vector for the derivatives.
xd = zeros(4, 1);
% Calculate the derivatives of the states. These equations can be copied
% directly from the sympy output but be sure to print with `mprint` in
% sympy.physics.mechanics to remove the `(t)` and use Matlab's find and
% replace and regular expressions to change the python `**` to the matlab
% `^`.
xd(1) = u1;
xd(2) = u2;
xd(3) = (-g * sin(q1) * sin(q2)^2 + 2 * g * sin(q1) - g * sin(q2) * ...
    cos(q1) * cos(q2) + 2 * l * u1^2 * sin(q1) * cos(q1) * ...
    cos(q2)^2 - l * u1^2 * sin(q1) * cos(q1) - 2 * l * u1^2 * ...
    sin(q2) * cos(q1)^2 * cos(q2) + l * u1^2 * sin(q2) * cos(q2) + ...
    l * u2^2 * sin(q1) * cos(q2) - l * u2^2 * sin(q2) * cos(q1)) / ...
    (l * (sin(q1)^2 * sin(q2)^2 + 2 * sin(q1) * sin(q2) * cos(q1) * ...
    cos(q2) + cos(q1)^2 * cos(q2)^2 - 2));
xd(4) = (-sin(q1) * sin(q2) / 2 - cos(q1) * cos(q2) / 2) * (2 * g * ...
    l * m * sin(q1) - l^2 * m * (-sin(q1) * cos(q2) + sin(q2) * ...
    cos(q1)) * u2^2) / (l^2 * m * (sin(q1) * sin(q2) / 2 + cos(q1) * ...
    cos(q2) / 2) * (sin(q1) * sin(q2) + cos(q1) * cos(q2)) - l^2 * m) + ...
    (g * l * m * sin(q2) - l^2 * m * (sin(q1) * cos(q2) - sin(q2) * ...
    cos(q1)) * u1^2) / (l^2 * m * (sin(q1) * sin(q2) / 2 + cos(q1) * ...
    cos(q2) / 2) * (sin(q1) * sin(q2) + cos(q1) * cos(q2)) - l^2 * m);
