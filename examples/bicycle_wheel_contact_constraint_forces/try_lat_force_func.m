const = [
 0.9534570696121849, % d1
 0.2676445084476887, % d2
 0.03207142672761929, % d3
 9.81, % g
 7.178169776497895, % ic11
 11.0, % ic22
 3.8225535938357873, % ic31
 4.821830223502103, % ic33
 0.0603, % id11
 0.12, % id22
 0.05841337700152972, % ie11
 0.06, % ie22
 0.009119225261946298, % ie31
 0.007586622998470264, % ie33
 0.1405, % if11
 0.28, % if22
 0.4707271515135145, % l1
 -0.47792881146460797, % l2
 -0.00597083392418685, % l3
 -0.3699518200282974, % l4
 85.0, % mc
 2.0, % md
 4.0, % me
 3.0, % mf
 0.35, % rf
 0.3, % rr
];

rear_wheel_radius = 0.3;
front_wheel_radius = 0.35;
% steady turn taken from Table 2 in Basu-Mandal 2007 row 2
roll_angle = 1.9178291654;
steer_angle = 0.4049333918;
rear_wheel_spin_rate = 10.3899258905;
rear_wheel_traversal_radius = 2.2588798195;

[q4, q5, q7, u11, u12, u4, u6, u7]

[A, b] = eval_dep_speeds(
res = A\b;
u3 = res(1);
u5 = res(2);
u8 = res(3);

rear_wheel_traversal_speed = rear_wheel_radius*rear_wheel_spin_rate;
yaw_rate = rear_wheel_traversal_speed/rear_wheel_traversal_radius;
% this is an estimate (not exact because front wheel doesn't have same
% traversal radius)
front_wheel_traversal_speed = rear_wheel_traversal_speed;
front_wheel_spin_rate = front_wheel_traversal_speed/front_wheel_radius;

var = [
roll_angle, %q4
0.0, %q5
steer_angle, %q7
yaw_rate, %u3
0.0, %u3p
0.0, %u4
0.0, %u4p
0.0, %u5
0.0, %u5p
rear_wheel_spin_rate, %u6
0.0, %u6p
0.0, %u7
0.0, %u7p
front_wheel_spin_rate, %u8
0.0 %u8p
];

[A, b] = eval_lat_forces(var, const);

x = A\b;

Fr = x(1)
Ff = x(2)
