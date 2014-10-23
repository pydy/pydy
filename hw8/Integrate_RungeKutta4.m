function [tout,yout] = Integrate_RungeKutta4(odefun,tspan,Ts,y0,par)
%[tout,yout] = Integrate_ModifiedEulerMethod(odefun,tspan,Ts,y0,par) 
%with tspan = [T0 TFINAL] integrates the system of differential equations 
%dy/dt = f(t,y) from time T0 to TFINAL, with initial conditions 
%specified by the row vector y0 and step size Ts, 
%using the modified Euler method 
%(D. Greenwood, "Advanced Dynamics", page 341).
%
%odefun must be a function handle to the function that calculates state
%derivatives dy/dt (as a column vector), when provided time t and a row or 
%column vector y. 
%Parameters can be passed on to odefun using the additional input par.
%
%Outputs:
%Each row in the solution array yout corresponds to a time 
%returned in the column vector tout.
%
%usage resembles the Matlab function ode45
%H. Vallery, Oct 2014

%time vector (t_1 to t_numsteps):
tout=(tspan(1):Ts:tspan(end))';

numsteps=length(tout);%number of steps 
numstates=length(y0);%number of states

%initialize outputs to allocate memory:
yout=zeros(numsteps,numstates);

%set initial values:
y_n=y0;%initial condition

for n=1:numsteps
%first generate outputs:
t_n=tout(n);
yout(n,:)=y_n;

%(1)
dy_n=feval(odefun,t_n,y_n',par)';%transpose output to produce row vector

%(2)
ytilde_n_plus_half = y_n + 0.5*Ts*dy_n;

%(3)
dytilde_n_plus_half = feval(odefun, t_n + 0.5*Ts, ytilde_n_plus_half', par)';

%(4)
ytilde2_n_plus_half = y_n + 0.5*Ts*dytilde_n_plus_half;

%(5)
dytilde2_n_plus_half = feval(odefun, t_n + 0.5*Ts, ytilde2_n_plus_half', par)';

%(6)
ytilde_n_plus_1 = y_n + Ts*dytilde2_n_plus_half;

%(7)
dytilde_n_plus_1 = feval(odefun, t_n + Ts, ytilde_n_plus_1', par)';

%(8)
y_n = y_n + Ts/6*(dy_n + 2*dytilde_n_plus_half + 2*dytilde2_n_plus_half +...
    dytilde_n_plus_1);
end


