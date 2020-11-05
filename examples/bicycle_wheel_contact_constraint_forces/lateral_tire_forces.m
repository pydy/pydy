function [Ff, Fr] = lateral_tire_forces(q, u, up, p)
% LATERAL_TIRE_FORCES - Returns the lateral tire constraint forces of the
% Carvallo-Whipple model.
%
% Syntax: [Ff, Fr] = lateral_tire_forces(q, u, up, p)
%
% Inputs:
%   q - Coordinates at time t, size 3x1 [q4 (roll angle); q7 (steer angle)]
%   u - Independent generalized speeds at time t, size 3x1 [u4 (roll rate);
%       u6 (rear wheel rate); u7 (steer rate)]
%   up - Derivatives of the independent generalized speeds at time t, size
%        3x1 [u4p (roll ang acc); u6p (rear ang acc); u7p (steer ang acc)]
%   p - Constant parameter structure with p items where p is the number of
%       parameters. Includes: [d1, d2, d3, g, ic11, ic22, ic31, ic33, id11,
%       id22, ie11, ie22, ie31, ie33, if11, if22, l1, l2, l3, l4, mc, md,
%       me, mf, rf, rr]
% Outputs:
%   Ff - Lateral constraint force at the front wheel contact.
%   Fr - Lateral constraint force at the rear wheel contact.

% unpack the function inputs
q4 = q(1);
q7 = q(2);

u4 = u(1);
u6 = u(2);
u7 = u(3);

u4p = up(1);
u6p = up(2);
u7p = up(3);

% all constants in alphabetic order
const = [p.d1, p.d2, p.d3, p.g, p.ic11, p.ic22, p.ic31, p.ic33, p.id11, ...
    p.id22, p.ie11, p.ie22, p.ie31, p.ie33, p.if11, p.if22, p.l1, p.l2, ...
    p.l3, p.l4, p.mc, p.md, p.me, p.mf, p.rf, p.rr];

% solve for the dependent pitch angle
f = @(q5) eval_holonomic([q4, q5, q7], [p.d1, p.d2, p.d3, p.rf, p.rr]);
% TODO : Add the correct lam calculate here.
lam = pi/10.0; % steer axis tilt
q5 = fsolve(f, lam);

% calculate the dependent generalized speeds
[A_nh, b_nh] = eval_dep_speeds([q4, q5, q7], ...
                               [u4, u6, u7], ...
                               [p.d1, p.d2, p.d3, p.rf, p.rr]);
res_nh = A_nh\b_nh;
u3 = res_nh(1);
u5 = res_nh(2);
u8 = res_nh(3);

% calculate the time derivatives of the dependent generalized speeds
[Ap_nh, bp_nh] = eval_dep_speeds_derivs([q4, q5, q7], ...
                                        [u3, u4, u5, u6, u7, u8], ...
                                        [u4p, u6p, u7p], ...
                                        const);
res = Ap_nh\bp_nh;
u3p = res(1);
u5p = res(2);
u8p = res(3);

% finally calculate the constraint forces
[A, b] = eval_lat_forces([q4, q5, q7], ...
                         [u3, u4, u5, u6, u7, u8], ...
                         [u3p, u4p, u5p u6p, u7p, u8p], ...
                         const);
x = A\b;
Ff = x(1);
Fr = x(2);

end
