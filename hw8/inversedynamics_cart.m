function [Fx, Fy] = inversedynamics_cart(theta, ddsX, ddsY, domega, par)

if length(theta) > 1
    theta = reshape(theta, [1, 1, length(theta)]);
    ddsX = reshape(ddsX, [1, 1, length(ddsX)]);
    ddsY = reshape(ddsY, [1, 1, length(ddsY)]);
    domega = reshape(domega, [1, 1, length(domega)]);
end

m = par.m*ones(size(theta));
a = 2*par.Is/par.length_cart;

A = [  m, 0*m, -a*sin(theta);
     0*m,   m,  a*cos(theta)];
B = [cos(theta), -2*sin(theta);
     sin(theta),  2*cos(theta)];
q = [ddsX; ddsY; domega];

Fx = zeros(length(theta), 1);
Fy = zeros(length(theta), 1);
for i = 1:length(theta)
    F = B(:, :, i)\A(:, :, i)*q(:, :, i);
    Fx(i) = F(1);
    Fy(i) = F(2);
end
