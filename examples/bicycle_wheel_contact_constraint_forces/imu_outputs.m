function [C_angvel, E_angvel, P_acc, Q_acc] = imu_outputs(q, u, up, p)
% IMU_OUTPUTS - Returns the IMU angular velocity and acceleration measures
% given the generalized variables of the Carvallo-Whipple model. The IMU on
% the rear frame is assumed to have its x axis parallel to the ground and
% pointing forward, y axis parallel to the ground and pointing rightward,
% and the z axis normal to the ground pointing downward when in the no roll,
% no steer configuration. The IMU on the front frame is assumed to have it's
% z axes parallel to the steer axis pointing downward and the x and y axes
% normal to the steer axis with the x axis pointing mostly in the forward
% direction and the y rightward.
%
% Syntax: [C_angvel, E_angvel, P_acc, Q_acc] = imu_outputs(q, u, up, p)
%
% Inputs:
%   q - Coordinates at time t, size 2x1, [q4 (roll angle), q7 (steer angle)]
%   u - Generalized speeds at time t, size 4x1 [u3 (yaw rate), u4 (roll
%       rate), u6 (rear wheel rate), u7 (steer rate)]
%   up - Derivatives of the generalized speeds at time t, size 4x1 [u3p (yaw
%        ang acceleration), u4p (roll ang acceleration), u6p (rear ang
%        acceleration), u7p (steer ang acceleration)]
%   p - Constant parameter structure which includes: [bx, by, bz, d1, d2,
%       d3, ex, ey, ez, g, rf, rr]. (bx, by, bz) and (ex, ey, ez) are the
%       distances along the IMU axes from the rear wheel or front wheel to
%       the respective IMU.
% Outputs:
%   C_angvel - Body fixed angular velocity of the rear frame as reported by the IMU.
%   E_angvel - Body fixed angular velocity of the front frame as reported by the IMU.
%   P_acc - Body fixed linear acceleration of a point P on the rear frame,
%           distance bx, by, and bz from the rear wheel center, as reported
%           by the IMU. (includes gravity!)
%   Q_acc - Body fixed linear acceleration of a point Q on the front frame,
%           distance ex, ey, and ez from the front wheel center, as reported
%           by the IMU. (includes gravity!)

    q4 = q(1);
    q7 = q(2);
    [q5, lam] = solve_for_pitch(q4, q7, p.d1, p.d2, p.d3, p.rf, p.rr);

    u3 = u(1);
    u4 = u(2);
    u5 = 0.0; % assume neglible pitch motion
    u6 = u(3);
    u7 = u(4);

    u3p = up(1);
    u4p = up(2);
    u5p = 0.0; % assume neglible pitch motion
    u6p = up(3);
    u7p = up(4);

    res = eval_imu([q4, q5, q7], ...
                   [u3, u4, u5, u6, u7], ...
                   [u3p, u4p, u5p, u6p, u7p], ...
                   [p.bx, p.by, p.bz, p.d1, p.d2, p.d3, p.ex, p.ey, p.ez, p.g, lam, p.rr]);

    C_angvel = res(1:3);
    P_acc = res(4:6);
    E_angvel = res(7:9);
    Q_acc = res(10:12);

end
