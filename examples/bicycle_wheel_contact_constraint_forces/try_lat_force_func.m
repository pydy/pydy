% This script shows how to use the lateral_tire_forces function.

% Load some geometry, mass, inertia values. These are the benchmark bicycle
% values converted to the parameter set in [Moore2012].
p.d1 = 0.9534570696121849;
p.d2 = 0.2676445084476887;
p.d3 = 0.03207142672761929;
p.g = 9.81;
p.ic11 = 7.178169776497895;
p.ic22 = 11.0;
p.ic31 = 3.8225535938357873;
p.ic33 = 4.821830223502103;
p.id11 = 0.0603;
p.id22 = 0.12;
p.ie11 = 0.05841337700152972;
p.ie22 = 0.06;
p.ie31 = 0.009119225261946298;
p.ie33 = 0.007586622998470264;
p.if11 = 0.1405;
p.if22 = 0.28;
p.l1 = 0.4707271515135145;
p.l2 = -0.47792881146460797;
p.l3 = -0.00597083392418685;
p.l4 = -0.3699518200282974;
p.mc = 85.0;
p.md = 2.0;
p.me = 4.0;
p.mf = 3.0;
p.rf = 0.35;
p.rr = 0.3;

% place IMUs at the mass centers of front and rear frames
p.bx = p.l1;
p.bz = p.l2;
p.ex = p.l3;
p.ez = p.l4;

% steady turn taken from Table 2 in Basu-Mandal 2007 row 2
%roll_angle = 1.9178291654;
%steer_angle = 0.4049333918;
%rear_wheel_spin_rate = 10.3899258905;
%rear_wheel_traversal_radius = 2.2588798195;

% try a steady turn
roll_angle = 5.0*pi/180;
steer_angle = 5.0*pi/180;
v = 10.0;  % forward speed of rear wheel center
rear_wheel_spin_rate = -v/p.rr;

q = [roll_angle, steer_angle];
u = [0.0, rear_wheel_spin_rate, 0.0];
up = [0.0, 0.0, 0.0];

[Ff, Fr] = lateral_tire_forces(q, u, up, p)

% a constant roll rate should give a centriptial acceleration in the z
% direction of the rear frame IMU (plus the gravitational component).
q = [0.0, 0.0];
u = [0.0, 1.0, 0.0, 0.0];
up = [0.0, 0.0, 0.0, 0.0];
%p.g = 0.0;  % don't include the nominal gravity measure
[C_angvel, E_angvel, P_acc, Q_acc] = imu_outputs(q, u, up, p)
