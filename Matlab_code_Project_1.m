%% Analytical Solution
close 
clear 
clc
% Define the grid parameters
x_min = 0;
x_max = 40;
dx = 1;
x = x_min:dx:x_max;
% Allocating memory for the solution
u = zeros(1, length(x));
% Evaluate the function
for j = 1:length(x)
 u(1, j) = 0.5 * (1 + tanh(250 * ((x(j) - 20) - 0.5 *10)));
end
% Plot the solution
figure;
plot(x, u,'-s','linewidth',2);
hold on;
% Defining the domain
L = 40; % Length of domain - here the length considered 40
T = 10; % Final computational time
Nx = 41; % Number of spatial grid points
dx = L/(Nx-1); % Spatial grid size
c = 0.5; % Advection velocity
dt = 0.5; % Time step
Nt = ceil(T/dt); % Number of time steps
x = linspace(0,L,Nx+1); % Spatial grid
u = zeros(Nx+1,Nt+1); % Solution matrix
% Initial condition
u(:,1) = 0.5*(1+tanh(250*(x-20))); % as given in the problem statement 
% Dirichlet boundary condition
u(1,:) = u(1,1);
u(end,:) = u(end,1);
%% UPWIND SCHEME
% Time marching using upwind scheme
for n = 1:Nt
 for i = 2:Nx
    u(i,n+1) = u(i,n) -c*dt/dx * (u(i,n) - u(i-1,n));
 end
end
% Plot the solution at the final time
plot(x,u(:,end),'-*','linewidth',2)
hold on;
%% LAX SCHEME
% Time marching using Lax scheme
for n = 1:Nt
 for i = 2:Nx
    u(i,n+1) = (0.5*(u(i+1,n)+u(i-1,n))) -((c/2)*dt/dx * (u(i+1,n) - u(i-1,n)));
 end
end
% Plot the solution at the final time
plot(x,u(:,end),'-+','linewidth',2)
hold on;
%% LEAP FROG SCHEME
% Time marching using Leap Frog scheme
for n = 2:Nt
 for i = 2:Nx
    u(i,n+1) = u(i,n-1)-(c*dt/dx)*(u(i+1,n) - u(i-1,n));
 end
end
% Plot the solution at the final time
plot(x,u(:,end),'-x','linewidth',2)
hold on;
%% LAX-WENDROFF SCHEME 
% Time marching using Lax-Wendroff scheme
for n = 1:Nt
 for i = 2:Nx
    u(i,n+1) = u(i,n)-(c*dt/(2*dx))*(u(i+1,n)-u(i-1,n))+(c^2*dt^2/(2*dx^2))*(u...
    (i+1,n)-2*u(i,n)+u(i-1,n));
 end
end
% Plot the solution at the final time
plot(x,u(:,end),'-o','linewidth',2)
hold on;
%% MAcCORMACK SCHEME
% Time marching using MacCromack scheme
for n = 1:Nt
 for i = 2:Nx
 up(i,n+1) = u(i,n)- c*dt/dx * (u(i+1,n) - u(i,n));
 u(i,n+1) = 0.5*(up(i,n+1)+u(i,n) - c*dt/dx * (up(i,n+1) - up(i-1,n+1)));
 end
end
% Plot the solution at the final time
plot(x,u(:,end),'-d','linewidth',2)
xlabel('x')
ylabel('u(x,t)')
title(['1D linear advection equation (comparison between analytical ...' ...
    'solution and all schemes at dt = 0.5)'],'FontSize',14)
legend('Analytical Solution','Upwind Scheme','Lax Scheme','Leap Frog scheme'...
    ,'Lax_Wendroff Scheme', 'MacCormack Scheme');
