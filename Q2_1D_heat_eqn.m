% -------------------------------------------------------------------------
% This is a script that plots the Numerical and Analytical Solutions
% for the 1D Heat Equation.

% Author: Joshua Selvanayagam
% 
% Modified: 2022-04-04
% -------------------------------------------------------------------------

clear all; clc; close all;

% Parameter definitions ---------------------------------------------------
Length = 1; % x in (0,L)                                (*@\label{code:q2_num_param_start}@*)
TotalTime = 10; % t in (0,t) defined by project outline
c = 2; % conductivity 
N = 100; % number of grid points
M = 400000; % number of time points
dx=Length/N; dt=TotalTime/M; % grid spacing
alpha = c*dt/dx^2;                                      (*@\label{code:q2_num_param_end}@*)

% Initialization ----------------------------------------------------------
% solution grid
x = linspace(0,Length,N+1);
T = ones(N+1,M);

% Initial Condition     (*@\label{code:q2_initial_cond_start}@*)
T(:,1) = cos(pi*x);     (*@\label{code:q2_initial_cond_end}@*)

% Dirichlet Boundary conditions at either end   (*@\label{code:q2_boundary_cond_start}@*)
T(1,:) = 0;
T(end,:) = 2;                                   (*@\label{code:q2_boundary_cond_end}@*)

% Matrix Initialization for Analytical Solution
V = ones(N+1,M);
V(:,1) = cos(pi*x);
V(1,:) = 0;
V(end,:) = 2;

% Numerical Solution ------------------------------------------------------
% Space: discretized using second order derivatives (*@\label{code:q2_loop_start}@*)
% Time: explicit time advancement, 
% i.e, k+1 values are determined by k values
% BCs: Dirichlet (constant temperature)
for k=1:M-1 % time
    for i=2:N % space
        T(i,k+1)=T(i,k)+alpha*(T(i+1,k)-2*T(i,k)+T(i-1,k));
    end
end (*@\label{code:q2_loop_end}@*)

% Analytical Solution -----------------------------------------------------
v = 0;
individualSum = 0;
time = 0;
space = 0;
D_const = 0;

for j=1:M-1 %time points
   time = TotalTime * j/M; %scale time based on 0 to 10 seconds
   for grid_space=2:N %space points
       space = grid_space./100; %scale grid based on 0 to 1 unit length
       v = 0;
        for n_sum=2:100 %summation of Fourier Series
            %Calculate D_n constant at each n
            D_const = (2.*((pi.*(n_sum.^3).*(cos(pi.*n_sum)+1)...
            ./((n_sum.^2)-1))-(2.*sin(pi.*n_sum))...
            +(2.*pi.*n_sum.*cos(pi.*n_sum))))./((pi.^2).*(n_sum.^2));
            
            %Calculate the partial sum at each n
            individualSum = D_const.*(exp((-2.*((n_sum.*pi).^2)).*time)...
            .*sin(n_sum.*pi.*space));
        
            %Add the partial sum to the total after each iteration of n
            v = v + individualSum;
        end
        V(grid_space, j+1) = 2.*space + v;
   end
end

% Visuaization ------------------------------------------------------------
% line plot
figure(1)
plot(T(:,4000)); hold on 
plot(T(:,40000))
plot(T(:,400000))
xlabel('x'); ylabel('T(x,t)');
legend('t=0.1','t=1','t=10')

figure(2)
plot(V(:,4000)); hold on
plot(V(:,40000))
plot(V(:,400000))
xlabel('x'); ylabel('T(x,t)');
legend('t=0.1','t=1','t=10')