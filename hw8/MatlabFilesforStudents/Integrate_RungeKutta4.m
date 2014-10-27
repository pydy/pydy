function [tout,yout] = Integrate_RungeKutta4(odefun,tspan,Ts,y0,par)
%[tout,yout] = Integrate_RungeKutta4(odefun,tspan,Ts,y0,par) 
%with tspan = [T0 TFINAL] integrates the system of differential equations 
%dy/dt = f(t,y) from time T0 to TFINAL, with initial conditions 
%specified by the row vector y0 and step size Ts, 
%using the modified Euler method 
%(D. Greenwood, "Advanced Dynamics", page 343).
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
%Modified by Daniel Lemus, Oct 2014

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
    t_nplus1=t_n+Ts;
    yout(n,:)=y_n;

    %(1) calculate function derivative:
    dy_n=feval(odefun,t_n,y_n',par)';%transpose output to produce row vector

    %(2) calculate 1st predictor:
    y_tildeplus1half= y_n+0.5*Ts*dy_n ;
    
    %(3) calculate function derivative evaluated in the 1st predictor:
    dy_tildeplus1half= feval(odefun,t_n+0.5*Ts,y_tildeplus1half',par)';
    
    %(4) calculate 2nd predictor:
    y_tildetildeplus1half= y_n + 0.5*Ts*dy_tildeplus1half ;
    
    %(5) calculate function derivative evaluated in the 2nd predictor:
    dy_tildetildeplus1half= feval(odefun,t_n+0.5*Ts,y_tildetildeplus1half',par)';
    
    %(6) calculate 3rd predictor:
    y_tildeplus1= y_n + Ts*dy_tildetildeplus1half;
    
    %(7) calculate function derivative evaluated in the 3rd predictor:
    dy_tildeplus1 = feval(odefun,t_n+Ts,y_tildeplus1',par)';
    
    %(8) calculate function derivative evaluated in the 3rd predictor:
    y_nplus1 = y_n + (Ts/6)*(dy_n+2*dy_tildeplus1half+ ...
               2*dy_tildetildeplus1half+dy_tildeplus1);

    %re-define n+1 back to n:
    y_n=y_nplus1;
end


