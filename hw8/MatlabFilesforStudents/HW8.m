%% Exercise 1

m = 2;
l = 9.81;
g = 9.81;

phidd = g*sin(phi)/l;

%% Exercise 2
clc 
clear all
[t,y0] = Integrate_Euler(@vanderpoldemo,[0 15],0.01,[2 0],1);
[t,y1] = Integrate_ModifiedEuler(@vanderpoldemo,[0 15],0.01,[2 0],1);
[t,y2] = Integrate_RungeKutta4(@vanderpoldemo,[0 15],0.01,[2 0],1);

disp('Final value = [Euler , Modified Euler, RK4]')
[y0(end,:) y1(end,:) y2(end,:)]

plot(t,y0(:,1))
hold all
plot(t,y1(:,1))
plot(t,y2(:,1))
legend('Euler','Modified Euler','RK4')

%% Exercise 3

