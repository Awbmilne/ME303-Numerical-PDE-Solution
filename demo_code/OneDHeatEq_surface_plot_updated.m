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
L = 1; % x in (0,L)
Tf = 0.5; % t in (0,T)
k = 1; % conductivity 
N = 20; % number of grid points
M = 5000; % number of time points
dx=L/N; dt=Tf/M; % grid spacing
alpha = k*dt/dx^2; 

% Initialization ----------------------------------------------------------
% solution grid
x = linspace(0,L,N+1);
T = ones(N+1,M);

% Initial Condition
T(:,1) = sin(pi*x);

% Dirichlet Boundary conditions at either end
T(1,:) = 0;
T(end,:) = 0;

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
% surface plot
figure(2)
[X,Y] = meshgrid(0:dx:L,dt:dt:Tf);
mesh(X,Y,T')
colormap('hot')
xlabel('x'); ylabel('t'); zlabel('T(x,t)'); 
colorbar

