% -------------------------------------------------------------------------
% This is a script that calculates how heat diffuses across a square domain. 
% The boundary conditions are all Dirichlet, i.e., the boundaries are at a
% constant temperature
%
% Modified: 2021-07-11
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

% Initial Conditions 
T(:,:,1) = sin(5*pi*X).*cos(4*pi*Y); % at (x,y,t=0)

% Dirichlet boundary conditions (these could also be applied in the loop)
T(1:end,1,:) = repmat(sin(pi.*y)',1,1,M);     % at (x=0,y,t)
T(1,:,:) = 0;                                 % at (x,y=0,t)
T(end,:,:) = 0;                               % at (x,y=1,t)
T(:,end,:) = repmat(sin(pi.*y)',1,1,M);       % at (x=1,y,t)

% PDE Solution ------------------------------------------------------------
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
end

% Visualization -----------------------------------------------------------
% Temperature disitribution at t = 0
figure()
surf(x,y,T(:,:,1))
view(2)
axis('equal')
title('time=0 s')
h = colorbar();
title(h,'Temperature')
xlabel('X')
ylabel('Y')

% Temperature disitribution at t = 0.15
figure()
surf(x,y,T(:,:,M))
view(2)
a=sprintf('time=%f s' ,finalTime);
title(a)
h = colorbar();
title(h,'Temperature')
axis('equal')

[gx,gy]=gradient(T(:,:,M));

div=divergence(X,Y,gx,gy);

figure()
quiver(X,Y,gx,gy,3,'red')

figure()
surf(X,Y,div)
