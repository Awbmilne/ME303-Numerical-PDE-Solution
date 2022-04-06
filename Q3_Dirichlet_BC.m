% -------------------------------------------------------------------------
% This script calculates heat diffusion on a square domain using Dirichlet Boundary Conditions.
% -------------------------------------------------------------------------

clear all; clc; close all;

% Parameter definitions ---------------------------------------------------
N = 30; % numbers of nodes along both x and y directions
L = 1; % Length of square domain
alpha = 1;  % thermal diffusivity
dt = 0.0002; % time step
dx = L/N; % grid spacing in x direction
dy = dx; % grid spacing in y direction (same as x direction)
finalTime = 0.15; % final time given in the problem
M = int16(finalTime/dt); % number of time steps required

% Initialization ----------------------------------------------------------
% Solution Grid
x=linspace(0,L,N);
y=linspace(0,L,N);
[X,Y]=meshgrid(x,y);
T = ones([size(X) M]); % all time steps are stored for demo purposes

% Initial Conditions (*@\label{code:q3_dirichlet_ic_start}@*)
T(:,:,1) = sin(pi*X).*sin(4*pi*Y); % at (x,y,t=0)

% Dirichlet boundary conditions (these could also be applied in the loop)
T(1:end,1,:) = 0;                             % at (x=0,y,t)
T(1,:,:) = repmat(sin(pi.*x)',1,1,M);         % at (x,y=0,t)
T(end,:,:) = repmat((cos(2*pi.*x)-1)',1,1,M); % at (x,y=1,t)
T(:,end,:) = 0;                               % at (x=1,y,t)  % (*@\label{code:q3_dirichlet_ic_end}@*)

% PDE Solution ------------------------------------------------------------  (*@\label{code:q3_dirichlet_solving_start}@*)
% Space: second order derivatives
% Time: explicit time marching, i.e, k+1 values are determined by k values
for k=1:M-1  % time step up to specified final time
    for i=2:length(y)-1 % x node step
        for j=2:length(x)-1 % y node step
            % Solving for temperature at next time step
            T(i,j,k+1) = ((((T(i,j+1,k)-2*T(i,j,k)+T(i,j-1,k))/dx^2)....
            +((T(i+1,j,k)-2*T(i,j,k)+T(i-1,j,k))/dy^2))*dt)+T(i,j,k);
        end        
    end  
end %(*@\label{code:q3_dirichlet_solving_end}@*)

% Visualization -----------------------------------------------------------
% Temperature distribution at t = 0
figure()
surf(x,y,T(:,:,1))
view(2)
axis('equal')
title('time=0 s')
h = colorbar();
title(h,'Temperature')
xlabel('X')
ylabel('Y')

% Temperature distribution at t = 0.05
figure()
surf(x,y,T(:,:,250))
view(2)
axis('equal')
title('time=0.05 s')
h = colorbar();
title(h,'Temperature')
xlabel('X')
ylabel('Y')


% Temperature distribution at t = 0.1
figure()
surf(x,y,T(:,:,500))
view(2)
axis('equal')
title('time=0.10 s')
h = colorbar();
title(h,'Temperature')
xlabel('X')
ylabel('Y')

% Temperature distribution at t = 0.15
figure()
surf(x,y,T(:,:,M))
view(2)
axis('equal')
title('time=0.15 s')
h = colorbar();
title(h,'Temperature')
xlabel('X')
ylabel('Y')
