% -------------------------------------------------------------------------
% This is a script that calculates how heat diffuses in a insulated
% circular rod heated at one end. 
% 
% The insulated end and the heated end are Neumann boundaries, i.e., the
% heat flux across them (dependent on the temperature gradient) is constant
%
% Original code: http://people.uncw.edu/hermanr/pde1/NumHeatEqn.pdf
% Author: R. L. Herman
% 
% Modified: 2021-07-11
% -------------------------------------------------------------------------

clear all; clc; close all;

% Parameter definitions ---------------------------------------------------
L= 1; % x in (0,L)
Tf =1; % t in (0,T)
k=1; % conductivity 
N = 100; % number of grid points
M = 20000; % number of time points
dx=L/N; dt=Tf/M; % grid spacing
alpha = k*dt/dx^2; 

% Initialization ----------------------------------------------------------
% solution grid
x = linspace(0,L,N+1);
T = ones(N+1,M);

% Initial Condition
T(:,1) = sin(pi*x);

% Dirichlet Boundary conditions at either end
T(1,:) = 1;
T(end,:) = -0.5;

% PDE Solution ------------------------------------------------------------
% Space: discretized using second order derivatives
% Time: explicit time marching, i.e, k+1 values are determined by k values
% BCs: Dirichlet (constant temperature)
for k=1:M-1 % time
    for i=2:N % space
        T(i,k+1)=T(i,k)+alpha*(T(i+1,k)-2*T(i,k)+T(i-1,k));
    end
end

% Visuaization ------------------------------------------------------------
% line plot
figure(1)
plot(T(:,1)); hold on 
plot(T(:,10))
plot(T(:,100))
plot(T(:,1000))
plot(T(:,10000))
plot(T(:,20000))
xlabel('x'); ylabel('T(x,t)');
legend('t=0','t=0.0005','t=0.005','t=0.05','t=0.5','t=1')

