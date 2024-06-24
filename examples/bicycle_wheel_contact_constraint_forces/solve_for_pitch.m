function [q5, lam] = solve_for_pitch(q4, q7, d1, d2, d3, rf, rr)
% SOLVE_FOR_PITCH - Returns the pitch angle, q5, for any given roll and
% steer angles, q4 and q7.
%
% Syntax: q5 = solve_for_pitch(q4, q7, d1, d2, d3, rf, rr)
%
% Inputs:
%   q4 - Roll steer angle [rad]
%   q7 - Steer angle [rad]
%   d1 - Distance from steer axis to rear wheel center [m]
%   d2 - Distance between wheel centers along steer axis [m]
%   d3 - Distance from steer axis to front wheel center [m]
% Outputs:
%   q5 - Pitch angle [rad]
%   lam - Steer axis tilt [rad]

f = @(q5) eval_holonomic([q4, q5, q7], [d1, d2, d3, rf, rr]);
g = @(lam) sin(lam) - (rf - rr + d2*cos(lam))/(d1 + d3);
guess = atan(d2 / (d1 + d3));  % guess based on equal wheel radii
lam = fsolve(g, guess); % steer axis tilt
q5 = fsolve(f, lam);

end
