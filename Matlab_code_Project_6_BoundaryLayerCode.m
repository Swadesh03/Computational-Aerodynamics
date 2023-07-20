% Copyright 2013 The MathWorks, Inc.

%% Parameters

function [dis_thick,cf_fine] = boundaryLayerCode(u_vel,Ma,le,te)
Temperature_freestream = 293; % Kelvin
pressure_freestream = 100000; % N/m^2
gas_constant = 287.058; % J/kg/K
rho_freestream = pressure_freestream/gas_constant/Temperature_freestream;% 1.1889; % kg/m^3
U = 301;   % m/s
mu = 1.79e-5; % kg/(m*s)
nu = mu/rho_freestream; % m^2/s
speed_of_sound = sqrt(1.4*gas_constant*Temperature_freestream);
Mach = Ma;
L = 1; % m
H = 0.01; % m
nx = 40;
ny = 41;

dx = L/(nx-1);
dy = H/(ny-1);
y = (0:dy:H);
x = (0:L/4:L);
del_star = zeros(1,nx);
[X,Y] = meshgrid(0:dx:L,0:dy:H);
Re_L = U*L/nu;
delta_L = 5/sqrt(U/nu/L);


%% Edge Velocity, U_e (one-dimensional vector from TSD code)
% U_e = ones(ny);
% U_e = Ma*speed_of_sound*U_e; 
U_e = u_vel(le+1:te,1);

v = zeros(ny,nx);
u = zeros(ny,nx);
%% Boundary Conditions
u(:,1) = U_e(1); % Incoming
v(:,1) = 0; % Incoming
u(1,:) = 0; % Bottom
v(1,:) = 0; % Bottom
u(ny,:) = U_e; % Top, free stream
%% Solve with Gaussian Elimination(TDMA specifically)
i = 0;
% March in the x direction
while i < nx-1
    i = i+1;
    % Determine parameters A,B,C,D in TDMA method
    for j = 2:ny-1
        if j == 2
            A(j) = 0;
            B(j) = (2*nu/dy^2) + u(j,i)/dx;
            C(j) = (v(j,i)/2/dy) - nu/dy^2;
            D(j) = u(j,i)^2/dx-(-nu/dy^2-v(j,i)/2/dy)*u(j-1,i);
        elseif j > 2 && j < ny-1
            A(j) = - nu/dy^2 - (v(j,i)/2/dy);
            B(j) = (2*nu/dy^2) + u(j,i)/dx;
            C(j) = (v(j,i)/2/dy) - nu/dy^2;
            D(j) = u(j,i)^2/dx;
        elseif j == ny-1
            A(j) = - nu/dy^2 - v(j,i)/2/dy;
            B(j) = 2*nu/dy^2 + u(j,i)/dx;
            C(j) = 0;
            D(j) = u(j,i)^2/dx-(v(j,i)/(2*dy)- nu/dy^2)*u(j+1,i);
        end
    end
    % solve for u with TDMA method
    usol = tdma(A(2:end),B(2:end),C(2:end),D(2:end));
    u(2:ny-1,i+1) = usol;
    % solve for v(j,i+1) based on known u
    for j = 2:ny
        v(j,i+1) = v(j-1,i+1) - dy/2/dx*(u(j,i+1)-u(j,i)+u(j-1,i+1)-u(j-1,i));
    end
end
% u_check = u;
% for i = 1:nx
%     for j = 1:ny
% 
%     if u_check>=0.99*U
%         break; 
%     end
%     end
% end

dis_thick = zeros(1,nx);
for i = 1:nx
    for j = 1:ny-1
        dis_thick(i) = dis_thick(i)+((1-u(j,i)/U_e(i))+(1-u(j+1,i)/U_e(i)))...
            *dy/2;
        %dis_thick(i) = dis_thick(i)+((1-u(j,i)/U)+(1-u(j+1,i)/U))...
            %*dy/2;
    end
end

cf_fine = zeros(1,nx);
for i = 1:nx
    dudy = (u(2,i)-u(1,i))/dy;
    tw = dudy*mu;
    cf_fine(i) = tw/(0.5*rho_freestream*U^2);
end

% del_point = ny;
% del = dy*del_point;
% del_star(1,i) = displacementThickness(u,U,H,nx,ny); 
% del_star(1,i) = 0;


%% Plotting
%u/U velocity contour
figure(1),
h1 = subplot(221);
set(h1,'XLim',[0 L],'YLim',[0 H],'NextPlot','replacechildren');
contourf(h1,X,Y,u/U_e(1),[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99]);
colorbar('peer',h1,'SouthOutside');
title(h1,'u/U Velocity Contour');
xlabel(h1,'x (m)'),ylabel(h1,'y (m)');
% velocity vector profile
h2 = subplot(222);
set(h2,'XLim',[0 L],'YLim',[0 H],'NextPlot','replacechildren');
[~,hdelta] = contour(h2,X,Y,u/U,0.99);
set(hdelta,'Color','r','LineWidth',3,'ShowText','On');hold(h2,'on');
quiver(h2,X,Y,u,v,'b');
title(h2,'Velocity Profile and Boundary Layer Thickness (0.99U)');
xlabel(h2,'x (m)'),ylabel(h2,'y (m)');hold(h2,'off');
% streamlines
h3 = subplot(223);
[strx,stry] = meshgrid(0,0:2*dy:H);
set(h3,'XLim',[0 L],'YLim',[0 H],'NextPlot','replacechildren');
streamline(h3,X,Y,u,v,strx,stry);
title(h3,'Streamlines');
xlabel(h3,'x (m)'),ylabel(h3,'y (m)');

% u vs y coordinate
% figure(2);
% plot (u(:,10),y,'linewidth',2);
% hold on;
% plot (u(:,20),y,'linewidth',2);
% plot (u(:,30),y,'linewidth',2);
% grid on; grid minor;
% xlabel('Velocity','FontSize',12); ylabel('y-coordinate','FontSize',12);
% title('Velocity vs y-corrodinate','FontSize',14);
% legend('25% chord length','50% chord length','75% chord length');
end