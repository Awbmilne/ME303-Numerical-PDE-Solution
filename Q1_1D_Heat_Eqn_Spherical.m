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
R = 0.02625; % Radius of egg [M]
T_start = 20; % Room temperature in celsius
T_min = 80; % Minimum cooking temp throughout
T_max = 100;
k = 0.5; % conductivity of egg [W K^-1 M^-1]
rho = 1035; % Density of egg [Kg M^-3]
c_p = 3.2*10e-3; % Specific heat of Egg [J Kg^-1 K^-1]
N = 100; % number of grid points
dt=0.001; % Size of time step
dr=R/N;  % grid spacing
% alpha = k / (rho*c_p) * dt * dr^2 / (2 * dr); % Increment coefficient
alpha = k / (rho*c_p) * dt * dr / (2 * dr); % Increment coefficient
% alpha = k * dt * dr^2 / (2 * dr); % Increment coefficient
% alpha = k * dt * dr / (2 * dr); % Increment coefficient

% Initialization ----------------------------------------------------------
% solution grid
x = linspace(0,R,N+1);
T = ones(N+1,1);

% Initial Condition
T(:,1) = T_start;

% Dirichlet Boundary conditions at either end
T(end,1) = T_max;


% PDE Solution ------------------------------------------------------------
% Space: discretized using second order derivatives
% Time: explicit time marching, i.e, k+1 values are determined by k values
% BCs: Dirichlet (constant temperature)
k = 1;
while any(T(:,k) < T_min) % While the minimum temp isnt reached
    T(1,k+1) = T(2,k); % Set center to same as 1 dr away from center
    T(end,k+1) = T(end,k); % Keep end value
    for i=2:N % space
        T(i,k+1)=T(i,k) + alpha*(T(i+1,k)-T(i-1,k)); % Increment T
        T(i,k+1) = (T(i,k+1) > T_max) * T_max + not(T(i,k+1) > T_max) * T(i,k+1); % Branchless limit T to T_max
    end
    time = k * dt;
    % fprintf("Time: %5.1f - Temp: %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f \n", time, T(1, k), T(N/100, k), T(N/10, k), T(N*5/10,k), T(N*9/10,k), T(N*99/100,k), T(N-1,k))
    if not(mod(time, 1))
        fprintf("Time: %5.1f - Temp: %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f \n", time, T(1, k), T(N/100, k), T(N/10, k), T(N*5/10,k), T(N*9/10,k), T(N*99/100,k), T(N-1,k))
    end
    k = k + 1;
end

size = length(T(1,:))
time = size * dt

% Visuaization ------------------------------------------------------------
% line plot
figure(1)
plot(T(:,1)); hold on 
plot(T(:,round(size/4)))
plot(T(:,round(size/2)))
plot(T(:,round(size*3/4)))
plot(T(:,round(size)))
% plot(T(:,1000))
% plot(T(:,10000))
% plot(T(:,20000))
xlabel('x'); ylabel('T(x,t)');
legend('t=0','t=0.0005','t=0.005','t=0.05','t=0.5','t=1')

