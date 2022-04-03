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
R = .02625; % Radius of egg [M]
T_start = 5; % Room temperature in celsius
T_min = 80; % Minimum cooking temp throughout
T_max = 100;
k = .500; % conductivity of egg [W K^-1 m^-1]
rho = 1035; % Density of egg [kg m^-3]
c_p = 3200; % Specific heat of Egg [J kg^-1 K^-1]
N = 100; % number of grid points
dt=0.001; % Size of time step
dr=R/N;  % grid spacing

alpha = k / (rho*c_p); % Increment coefficient
coeff = 1-2*alpha*dt/(dr^2);
fprintf("alpha = %f, coeff = %f \n", alpha, coeff);

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
while any(T(:,k) < T_min); % While the minimum temp isnt reached
    T(1,k+1) = max(T(2,k) - (T(3,k) - T(2,k)), T_start); % Set center to 1 'slope' away from i=2, Dont go down tho
    T(end,k+1) = T(end,k); % Keep end value
    for i=2:N % space
        r = (i-1) * dr;
        d2T_dr2 = (T(i+1,k)-2*T(i,k)+T(i-1,k))/(dr^2);
        dT_dr = (T(i+1,k)-T(i-1,k))/(2*dr);
        T(i,k+1)=T(i,k) + alpha*dt*(d2T_dr2 + (2/r)*dT_dr); % Increment T
    end
    time = k * dt;
    % fprintf("Time: %5.6f - Temp: %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f, %3.10f \n", time, T(1, k), T(N/100, k), T(N/10, k), T(N*5/10,k), T(N*9/10,k), T(N*99/100,k), T(end-1,k), T(end,k))
    if not(mod(time, 1))
        fprintf("Time: %5.6f - Temp: %3.14f, %3.14f, %3.14f, %3.14f, %3.14f, %3.14f, %3.14f \n", time, T(1, k), T(N/100, k), T(N/10, k), T(N*5/10,k), T(N*9/10,k), T(N*99/100,k), T(N-1,k))
    end
    k = k + 1;
end

size = length(T(1,:))
time = size * dt
T(:,end)

% Visuaization ------------------------------------------------------------
% line plot
indexes = [1, round(size/4), round(size/2), round(size*3/4), round(size)];
times = indexes * dt;
figure(1)
hold on;
plot(T(:,indexes(1)))
plot(T(:,indexes(2)))
plot(T(:,indexes(3)))
plot(T(:,indexes(4)))
plot(T(:,indexes(5)))
% plot(T(:,1000))
% plot(T(:,10000))
% plot(T(:,20000))
xlabel('x'); ylabel('T(x,t)');
legend(split(sprintf('t = %i sec,', round(times(:))), ","));
