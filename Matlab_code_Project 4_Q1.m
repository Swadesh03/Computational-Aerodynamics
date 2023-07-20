close 
clear
clc

%% Defining the meshing parameters and variables values
L = 1; % length of computational domain
nx = 50; % number of grid points
h = 0.15; % bump height
t4 = 0.8; t2 = 3;  % t1 = locates the maximum location of the bump
                   % t2 = controls the width of the bump
dx = L/nx;
%dz = 1;

% variables values

g = 1.4; % gamma
%r = 1716;  % gas constant in ft.lb/slug.R
%Tt = 531.2*r; % inlet total temperature
r = 286.9594;
Tt = 295.1111;             % third change
%Pt = 2117;  % inlet total pressure in lb/ft^2
Pt = 101362.5;    % fourth change
Pex = 0.72*Pt; % exit pressure

%% Defining the computational domain
% % faces
% x_faces = 0:dx:L+dx;

% face nodes in x direction
C_x = zeros(nx+1,1);
C_x(1) = 0;
for i = 2:nx+1
    C_x(i) =C_x(i-1)+dx;
end
%C_x = C_x*L/C_x(end);
node_x = zeros(nx+2,1);
for i = 2:nx+1
    node_x(i) = (C_x(i)+C_x(i-1))/2;
end    
node_x(end) = L;

% to generate the nozzle profile
% s_xf store the surface values at control volume faces
s_xf = zeros(nx+1,1);
s_xfn = zeros(nx+1,1);
for i = 1:nx+1
    s_xf(i) = 1-h*(sin(pi*(C_x(i))^t4))^t2;
end
s_xfn = -s_xf;
% s_xn store the surface values at control volume nodes
s_xn = zeros(nx+2,1);
for i = 1:nx+2
    s_xn(i) = 1-h*(sin(pi*(node_x(i))^t4))^t2;
end

%% Plotting the nozzle shape
% figure(1);
% plot(C_x,s_xf);
% hold on;
% plot(C_x,s_xfn)
% % axis([0 1 0.5 1]);


%% Initializing the state vectors and pressure
pho = zeros(nx+2,1); % Density
%M = zeros(nx+2,1);  % Momentum
E = zeros(nx+2,1);   % Energy
Ps =zeros(nx+2,1); % static pressure
Ts = zeros(nx+2,1);  % static temperature
Ma = zeros(nx+2,1);  % local mach no
u = zeros(nx+2,1);  % local air speed
c = zeros(nx+2,1);  % local speed of sound
t = zeros(nx+2,1);  % del_t
W = zeros(nx+2,3); 
F = zeros(nx+2,3);
Fp = zeros(nx+2,3);
Fm = zeros(nx+2,3);
Q = zeros(nx+2,3);
sp = zeros(nx+2,1);
sm = zeros(nx+2,1);
v = zeros(nx+1,1);
R = zeros(nx+2,3);
CFL = 0.3;   % CFL value
z = 0.2;    % Epsilon value
Cv = r/(g-1);
R_o = R;
pho_r = pho;
Ps_r = Ps;
Ma_r = Ma;
% Initializing temperature & Mach no
for i = 1:nx+2
    Ts(i) = Tt;    % Temperature Initialization
    Ma(i) = 0.1;   % Mach number initialization
end 

% Initializing Pressure
for i = 1:nx+2
    if i == nx+2
        Ps(i) = Pex;
    else
        Ps(i) = Pt;
    end
end

% Initializing Density, Momentum and Energy
for i = 1:nx+2
    c(i) = sqrt(g*r*Ts(i));  % local speed of sound
    u(i) = Ma(i)*c(i);   % air speed
    pho(i) = Ps(i)/(r*Ts(i));  % density
   % M(i) = pho(i)*u(i);   % momentum
    E(i) = pho(i)*Cv*Ts(i)+(0.5*u(i)^2);  % energy
    %E(i) = pho(i)*Cv*Tt;                   % first change - Ts - Tt
end
%% Initialization W, F, & Q
W (:,1 ) = pho; W (:,2 ) = pho.*u; W (:,3 ) = E;
F (:,1 ) = pho.*u; F (:,2 ) = (pho.*u.*u)+Ps; F (:,3 ) = (E+Ps).*u;
%F (:,1 ) = pho*u; F (:,2 ) = (pho*u*u)+Ps; F (:,3 ) = (E+Ps)*u;
for i = 2:nx+1
    sp(i) = 0.5*(s_xn(i)+s_xn(i+1)); 
    sm(i) = 0.5*(s_xn(i)+s_xn(i-1));
    Q(i,2) = Ps(i)*(sp(i)-sm(i))/(s_xn(i)*dx);
    %Q(i,2) = (Ps(i)/s_xn(i))*h*t4*(sin(pi*(node_x(i))^t4))^(t2-1)*...
        %cos(pi*(node_x(i))^t4)*t4*pi*node_x(i)^(t4-1);
end


%% Iteration Parameters
itr = 0;
tol = 1E-14;
error = 1;

error_store = zeros(itr,1); % storing the errors
error_store_pho = zeros(itr,1); % storing the error in density
error_store_Ps = zeros(itr,1); % storing the error in Pressure
error_store_Ma = zeros(itr,1); % storing the error in Mach
itr_store = zeros(itr,1); 

%% Analysis Loop
while error>tol
    itr = itr+1;

    % Analysis steps
    for j = 1:3
        for i = 2:nx+1
            sp(i) = 0.5*(s_xn(i)+s_xn(i+1)); 
            sm(i) = 0.5*(s_xn(i)+s_xn(i-1));

            t(i) = (CFL*(C_x(i)-C_x(i-1)))/(u(i)+c(i));

            % Finding the Eigen values
            y1 = 0.5*(u(i)+u(i+1));
            y2 = (0.5*(u(i)+u(i+1)) + 0.5*(c(i)+c(i+1)));
            y3 = (0.5*(u(i)+u(i+1)) - 0.5*(c(i)+c(i+1)));
            yp = max([y1 y2 y3]);

            y4 = 0.5*(u(i)+u(i-1));
            y5 = (0.5*(u(i)+u(i-1)) + 0.5*(c(i)+c(i-1)));
            y6 = (0.5*(u(i)+u(i-1)) - 0.5*(c(i)+c(i-1)));
            ym = max([y4 y5 y6]);

            % Computation of the flux
            Fp(i,j) = 0.5*(F(i,j)+F(i+1,j))-(0.5*z*abs(yp)*(W(i+1,j)-W(i,j)));
            Fm(i,j) = 0.5*(F(i,j)+F(i-1,j))-(0.5*z*abs(ym)*(W(i,j)-W(i-1,j)));
            
            v(i) = s_xn(i)*dx;
            R(i,j) = ((Fp(i,j)*sp(i))-(sm(i)*Fm(i,j)))...
                -(s_xn(i)*dx*Q(i,j));  % 
            W(i,j) = W(i,j) - ((t(i)*R(i,j))/v(i)); 
              
        end
    end
    % Finding the updated value of Flow parameters
    for i = 2:nx+1
        pho(i) = W(i,1);   % density
        u(i) = W(i,2)/pho(i);  % velocity
        E(i) = W(i,3);   % Energy
        Ps(i) = (g-1)*pho(i)*((E(i)/pho(i))-(u(i)^2/2));
        Ts(i)= (E(i)-0.5*u(i)^2)/(pho(i)*Cv);    % temperature  
        c(i) = sqrt(g*r*Ts(i));   % local speed of sound
        Ma(i) = u(i)/c(i);  % Mach number
    end
    %% Inlet Boundary Condition
    if Ma(1) < 1
        a = 2*g*((g-1)/(g+1))*Cv*Tt;
 
        % Compute dp/du
        a1 = (-2*(g-1)*u(1))/((g+1)*a);
        b1 = ((g-1)*u(1)^2)/((g+1)*a);
        c1 = (1-b1)^(1/(g-1));
        DPu = Pt*(g/(g-1))*c1*a1;

        % compute del_u
        t1 = (CFL*dx)/(u(1)+c(1));
        d1 = (((u(2)+u(1))/2) -((c(2)+c(1))/2));
        e1 = d1*t1/dx;
        du_i = -e1*(Ps(2)-Ps(1)-(pho(1)*c(1)*(u(2)-u(1))))/(DPu-(pho(1)*c(1)));

        % updating the flow properties at inlet
        b1 = ((g-1)*u(1)^2)/((g+1)*a);
        Ts(1) = Tt*(1-b1);   % temperature
        u(1) = u(1)+du_i;    % velocity
        Ps(1) = Pt*(Ts(1)/Tt)^(g/(g-1));  % pressure
        pho(1) = Ps(1)/(r*Ts(1));   % density
        E(1) = pho(1)*((Cv*Ts(1))+(0.5*u(1)^2));  % Energy
        c(1) = sqrt((g*Ps(1))/pho(1));   % local speed of sound
        Ma(1) = u(1)/c(1);
    end
    %% Outlet Boundary Condition

    % Compute eigen values
    t_end = (CFL*dx)/(u(end)+c(end));
    f1 = ((u(end)+u(end-1))/2)*(t_end/dx);
    f2 = (((u(end)+u(end-1))/2)+((c(end)+c(end-1))/2))*(t_end/dx);
    f3 = (((u(end)+u(end-1))/2)-((c(end)+c(end-1))/2))*(t_end/dx);

    % Compute Characteristic relations
    R1 = -f1*(pho(end)-pho(end-1)-((1/c(end)^2)*(Ps(end)-Ps(end-1))));
    R2 = -f2*(Ps(end)-Ps(end-1)+(pho(end)*c(end)*(u(end)-u(end-1))));
    R3 = -f3*(Ps(end)-Ps(end-1)-(pho(end)*c(end)*(u(end)-u(end-1))));

    % Compute Mach no - 
    Ma(end) = ((u(end)+u(end-1))/2)/((c(end)+c(end-1))/2);

    % Compute del_p
    if Ma(end) > 1
        dp = 0.5*(R2+R3);
    else
        dp = 0;
    end

    % Update dpho & du_outlet
    dpho = R1 + (dp/c(end)^2);
    du_o = (R2 - dp)/(pho(end)*c(end));

    % Update the flow properties 
    pho(end) = pho(end)+dpho;
    u(end) = u(end)+du_o;
    Ps(end) = Ps(end)+dp;
    Ts(end) = Ps(end)/(pho(end)*r);
    E(end) = pho(end)*((Cv*Ts(end))+(0.5*u(end)^2));
    c(end) = sqrt((g*Ps(end))/pho(end));
    Ma(end) = u(end)/c(end);

    % Updating W, E, Q
    W (:,1 ) = pho; W (:,2 ) = pho.*u; W (:,3 ) = E;
    F (:,1 ) = pho.*u; F (:,2 ) = (pho.*u.*u)+Ps; F (:,3 ) = (E+Ps).*u;
    for i = 2:nx+1
        sp(i) = 0.5*(s_xn(i)+s_xn(i+1)); 
        sm(i) = 0.5*(s_xn(i)+s_xn(i-1));
        Q(i,2) = Ps(i)*(sp(i)-sm(i))/(s_xn(i)*dx);
    end

    %% Check for convergence
    error = max(abs(R(:,1)-R_o(:,1)),[],"all");
    error_pho = max(abs(W(:,1)-pho_r),[],"all");
    error_Ps = max(abs(Ps-Ps_r),[],"all");
    error_Ma = max(abs(Ma-Ma_r),[],"all");
    error_store(itr,1) = error; 
    error_store_pho(itr,1) = error_pho;
    error_store_Ps(itr,1) = error_Ps;
    error_store_Ma(itr,1) = error_Ma;
    itr_store(itr,1) = itr;
    fprintf('%d %d\n',itr,error);
    R_o = R;
    pho_r = W(:,1);
    Ps_r = Ps;
    Ma_r = Ma;
end

%% Plotting the results

% Log(error) vs Residual
figure(1);
grid on; 
plot(itr_store,log10(error_store),'LineWidth',2);
xlabel('Iterations','FontSize',12);
ylabel('Log(error)','FontSize',12);
title('Convergence plot of Density residual', 'FontSize',14);

% Pressure ratio plot
figure(2);
grid on; 
plot(node_x,Ps/Pt,'LineWidth',1);
xlabel('X','FontSize',12);
ylabel('Ps/Pt','FontSize',12);
title('Pressure distribution along the channel', 'FontSize',14);

% Mach number plot
figure(3);
grid on; 
plot(node_x,Ma,'LineWidth',1);
xlabel('X','FontSize',12);
ylabel('Mach Number','FontSize',12);
title('Flow Mach number along the channel', 'FontSize',14);
