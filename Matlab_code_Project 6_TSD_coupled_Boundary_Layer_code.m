close
clear
clc

%% Defining the mesh
L = 50;
xf = 4; % factor by which mesh is refined compared to 10 grid points on airfoil
n1 = 80; % number of grid points before the leading edge
n2 = 40; % number of grid points for the airfoil
n3 = 120;  % number of grid points after the the trailing edge
le = n1; % leading edge x position
te = n1+n2; % trailing edge x position
tc = 0.08; % thickness ratio
nx = n1+n2+n3;
n4 = 163;
n = zeros(nx,1);
m = zeros(n4,1);
m(1) = 0;
m(2) = tc/2;  
n(n1) = n1/xf;
dx = 1/n2;
for i = n1+1:n1+n2
    n(i) = n(i-1) + dx; 
end
for i = n1+n2+1:nx
    n(i) = n(i-1) + (n(i-1)-n(i-2))*1.02858^1.05;
end
for i = n1:-1:2
    n(i-1) = n(i) - (n(i+1)-n(i))*1.04436^1.05;
end
% nt1 =zeros(nx,1);
% nt1(n1+n2+1:nx) = n(n1+n2+1:nx) - n1-1; nt2 = -flip(nt1);
% n(1:n1) = n1;
% n(1:n1) = n(1:n1)+(nt2(n1/2+1:n1+n2));
for i = 3:n4
    m(i) = m(i-1) + (m(i-1)-m(i-2))*1.01858^1.1;
end
% n = n/xf;

% n = n*L/n(nx,1);
n(1)=0;
m(end) = 50; n(end) = 50;
n_fine = n;
[x1,y1]= meshgrid(n,m);
x = x1'; y = y1';

%% Defining the boundary conditions and analysis parameters
mach = [0.8,0.84,0.86];
%mach = 0.8;
cp_store = zeros(length(mach),41);
Cdp_store = zeros(1,length(mach));
Clp_store = zeros(1,length(mach));
Clf_store = zeros(1,length(mach));
Cdf_store = zeros(1,length(mach));
Cl_store = zeros(1,length(mach));
Cd_store = zeros(1,length(mach));
pressure_store = zeros(length(mach),41);
%delta_star = zeros(length(mach),te-le+1);
cf_store = zeros(length(mach),te-le+1);
for k = 1:length(mach)
yg = 1.4; % gamma for air
nx = n1+n2+n3;
ny = n4;
Ma = mach(k);   % Mach no of the flow
Ta = 293;   % Temp
pa = 100;   % Pressure
R = 287.058;  % gas constant
tc = 0.08;
p_a = 0;
p = zeros(nx,ny);

ca = sqrt(yg*R*Ta);
ua = ca*Ma; % calculation of u infinity
a = zeros(nx,ny);
b = zeros(nx,ny);
c = zeros(nx,ny);
d = zeros(nx,ny);
e = zeros(nx,ny);
g = zeros(nx,ny);

eps = 1E-3;
error = 1;   % initial error value
er = 1; % error between the TSD and Boundary layer velocity

er_accept = 1E-4; % tolerance for TSD and Boundary layer
itr = 0;   % iteration
itr_b = 0; % iteration for boundary layer code and TSD
error_store = zeros(itr,1);   % storing the errors
er_store = zeros(itr_b,1); % storing the error for boundary layer and TSD
itr_store = zeros(itr,1);    % storing the iteration
itr_b_store = zeros(itr_b,1); % storing the iteration for boundary and TSD
disp_thick_old = zeros(1,te-le);
disp_thick_new = zeros(1,te-le);
delta_star = zeros(itr_b,te-le+1);
y_profile = zeros(1,te-le);
y_profile_new = zeros(1,te-le);
y_profile_old1 = y_profile;
y_profile_old2 = zeros(1,te-le);
dy_profile = zeros(1,te-le);
er_profile = zeros(1,te-le);
ue_storage = zeros(te-le+1,itr_b);
%% Boundary condition
while er>er_accept
    itr_b = itr_b + 1;
for i = 2 : nx-1
    if i <= le
        p(i,2) = p(i,1);
    elseif i > te
        p(i,2) = p(i,1);
    else
        if itr_b == 1
         %x_pos = n1 + (i-le)/(te-le);
         x_pos = n(i);
         p(i,1)=p(i,2)-ua*pos(x_pos)*(m(2)-m(1));
        else
            x_pos = n(i);
            p(i,1) = p(i,2)-ua*dy_profile(1,i-le)*(m(2)-m(1));
        end
    end
end

px = zeros(nx,ny);
u_vel = zeros(nx,ny);
for i = 2 : nx-1
    for j = 1: ny
        px(i,j) = (p(i+1,j)-p(i-1,j))/((n(i+1)-n(i-1)));
        u_vel(i,j) = ua+px(i,j);
    end
end

A = zeros(nx,ny);
for i = 1 : nx
    for j = 1 : ny
         A(i,j) = (1-Ma^2)-((yg+1)*Ma^2/ua*px(i,j));
    end
end

q = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        if A(i,j)>0
            q(i,j) = 0;
        elseif A(i,j)<0
            q(i,j) = 1;
        end
    end
end

p_gauss = p;
error = 1;

while error>eps
    itr = itr+1;
    for i = 3:nx-1
        for j = 2:ny-1
            c(i,j) = 2/((m(j)-m(j-1))*(m(j+1)-m(j-1)));
            g(i,j) = 2*q(i-1,j)*A(i-1,j)/((n(i-1)-n(i-2))*(n(i)-n(i-2)));
            d(i,j) = (2*(1-q(i,j))*A(i,j))/((n(i)-n(i-1))*(n(i+1)-n(i-1)))+...
            (-2*q(i-1,j)*A(i-1,j))/((n(i)-n(i-1))*(n(i)-n(i-2)))+...
            ((-2*q(i-1,j)*A(i-1,j))/((n(i-1)-n(i-2))*(n(i)-n(i-2))));
            a(i,j) = (-2*(1-q(i,j))*A(i,j))/((n(i+1)-n(i))*(n(i+1)-n(i-1)))+...
            (-2*(1-q(i,j))*A(i,j))/((n(i)-n(i-1))*(n(i+1)-n(i-1)))+...
            2*q(i-1,j)*A(i-1,j)/((n(i)-n(i-1))*(n(i)-n(i-2)))+...
            (-2)/((m(j+1)-m(j))*(m(j+1)-m(j-1)))+...
            (-2)/((m(j)-m(j-1))*(m(j+1)-m(j-1)));
            e(i,j) = (2*(1-q(i,j))*A(i,j))/((n(i+1)-n(i))*(n(i+1)-n(i-1)));
            b(i,j) = 2/((m(j+1)-m(j))*(m(j+1)-m(j-1)));
            p_gauss(i,j) = (-c(i,j)*p_gauss(i,j-1)-d(i,j)*p_gauss(i-1,j)...
                -e(i,j)*p(i+1,j)-b(i,j)*p(i,j+1)-g(i,j)*p_gauss(i-2,j))...
                /a(i,j);
        end
    end
    error = max(abs(p_gauss-p),[],"all");
    error_store(itr,1) = error; 
    itr_store(itr,1) = itr;
    %fprintf('%d %d\n',itr,error);
    for i = 2 : nx-1
        if i <= le
            p_gauss(i,1) = p_gauss(i,2);
        elseif i > te
            p_gauss(i,1) = p_gauss(i,2);
        else
            if itr_b == 1
            %x_pos = n1 + (i-le)/(te-le);
            x_pos = n(i);
            p_gauss(i,1)=p_gauss(i,2)-ua*pos(x_pos)*(m(2)-m(1));
            else
                x_pos = n(i);
                p_gauss(i,1) = p_gauss(i,2)-ua*dy_profile(1,i-le)*(m(2)-m(1));
            end
%             %x_pos = n1 + (i-le)/(te-le);
%             x_pos = n(i);
%             p_gauss(i,1)=p_gauss(i,2)-ua*pos(x_pos)*(m(2)-m(1));
        end
    end
    p = p_gauss;
    % Boundary condition update
    for i = 2 : nx-1
        for j = 1: ny
            px(i,j)= (p_gauss(i+1,j)-p_gauss(i-1,j))/((n(i+1)-n(i-1)));
            u_vel(i,j) = ua+px(i,j); 
        end
    end
    ue_storage(:,itr_b)=u_vel(le:te,1);
    for i = 1 : nx
        for j = 1 : ny
            A(i,j)= (1-Ma^2)-((yg+1)*Ma^2/ua*px(i,j));
        end
    end
    for i = 1:nx
        for j = 1:ny
            if A(i,j)>0
                q(i,j)= 0;
            elseif A(i,j)<0
                q(i,j)= 1;
            end
        end
    end

end
%% integrating the boundary layer code in TSD
[disp_thick_new,cf_fine] = boundaryLayerCode(u_vel,Ma,le,te);
delta_star(itr_b,1:end-1) = disp_thick_new;
delta_star(itr_b,end) = disp_thick_new(1,1);
    for i = 1:n2
        x_pos = n(i-1+le);
        if itr_b == 1
        y_profile(1,i) = tc*((-2*x_pos^2)+(82*x_pos)-840);
        y_profile_new(1,i) = y_profile(1,i)+ disp_thick_new(1,i);
        else
            %y_profile_old1(1,i) = y_profile(1,i);
        y_profile_new(1,i) = y_profile(1,i)+ disp_thick_new(1,i);
        end
    end
    for i = 2:n2-1
        
            dy_profile(1,i) = (y_profile_new(1,i+1) - y_profile_new(1,i-1))...
                /(n(i+le)- n(i-2+le));
            %er_profile(1,i) = y_profile_new(1,i)-y_profile_old2(1,i);
    end
    er = max(abs(disp_thick_new-disp_thick_old),[],"all");
    %er = max(abs(y_profile_new-y_profile_old2),[],"all");
    %y_profile_old1 = y_profile_new;
    er_store(itr_b,1) = er; 
    itr_b_store(itr_b,1) = itr_b;
    fprintf('%d %d\n',itr_b,er);
    disp_thick_old = disp_thick_new;
    %y_profile_old2 = y_profile_new;
end
% delta_star(k,1:end-1) = disp_thick_new;
% delta_star(k,end) = disp_thick_new(1,end);
cp_fine = zeros(1,te-le);
cf_store(k,1:end-1) = cf_fine;
cf_store(k,end) = cf_fine(1,end);
pressure_fine = zeros(1,te-le);
pr = zeros(1,te-le);
x_airfoil_fi = zeros(1,te-le);
for i= 1 : te-le+1
    x_airfoil_fi(i) = n1/4 + (i-1)/(te-le);
    cp_fine(i) = -2*px(i+le,1)/ua;
    pressure_fine(i) = ((cp_fine(i)*0.5*yg*Ma^2)+1)*pa; 
end
cp_store(k,:)= cp_fine;
pressure_store(k,:) = pressure_fine;

%% calculation of aerodynamic coefficient
Cpx = 0;
Cpy = 0;
Cfx = 0;
Cfy = 0;
alpha = 0;
alpha_r = alpha*pi/180;

for i = 2:n2-1
     x_pos = n(i-1+le);
end
for i = le+1:te-1
    x_pos = n_fine(i);
    theta = atan(dy_profile(1,i-le));
    sin_a = sin(theta);
    cos_a = cos(theta);
    dx = n_fine(i)-n_fine(i-1);
    dy = y_profile_new(1,i-le+1)- y_profile_new(1,i-le);
    dA = sqrt((dx)^2+(dy)^2);
    n5 = -sin_a; n6 = cos_a;
    t1 = cos_a; t2 = sin_a;

    Cpx = Cpx - (cp_fine(1,i-le)*n5*dA);
    Cpy = Cpy - (cp_fine(1,i-le)*n6*dA);     
    Clp = -Cpx*sin(alpha) + Cpy*cos(alpha);
    Cdp = Cpx*cos(alpha) + Cpy*sin(alpha);

    Cfx = Cfx + cf_store(k,i-le)*t1*dA;
    Cfy = Cfy +  cf_store(k,i-le)*t2*dA;
    Clf = -Cfx*sin(alpha_r) + Cfy*cos(alpha_r);
    Cdf = Cfx*cos(alpha_r) + Cfy*sin(alpha_r);

    Cl = Clf + Clp;
    Cd = Cdf + Cdp;
end

Clp_store(k,1)= Clp;
Cdp_store(k,1)= Cdp;
Clf_store(k,1)= Clf;
Cdf_store(k,1)= Cdf;
Cl_store(k,1)= Cl;
Cd_store(k,1)= Cd;

% %% Cf calculation
% cf_fine = zeros(1,te-le);
% pho = 1.1889;
% for i= 1 : te-le+1
%     mu = 1.79e-5;
%     dudy = (u_vel(i+le,2)-u_vel(i+le,1))/(y(i+le,2)-y(i+le,1)); 
%     tw = dudy*mu;
%     cf_fine(i)= tw/(0.5*pho*(ua)^2);
% end
% cf_store(k,:)= cf_fine;

%% Calculation of Cdp - pressure drag
% Cpx = 0;
% Cpy = 0;
% itr = 0;
% 
% alpha = 0;
% for i = le+1:te
%     x_pos = n_fine(i);
%     theta = atan(pos(x_pos));
%     sin_a = sin(theta);
%     cos_a = cos(theta);
%     dx = n_fine(i)-n_fine(i-1);
%     dy = pos_y(n_fine(i),tc)- pos_y(n_fine(i-1),tc);
%     dA = sqrt((dx)^2+(dy)^2);
%     n5 = -sin_a; n6 = cos_a;
% 
%     Cpx = Cpx - (cp_fine(1,i-le)*n5*dA);
%     Cpy = Cpy - (cp_fine(1,i-le)*n6*dA);
%     Clp = -Cpx*sin(alpha) + Cpy*cos(alpha);
%     Cdp = Cpx*cos(alpha) + Cpy*sin(alpha);
% end
% Cdp_store(1,k) = Cdp;
% end

r = 1;
while m(r) <= 1
    r = r+1;
end

y_airfoil = m(1:r);
%pressure_contour = zeros(size(x_airfoil,2),size(y_airfoil,2));
pressure_contour = zeros(length(x_airfoil_fi),length(y_airfoil));
 
for i = 1 : length(x_airfoil_fi)
    for j = 1: length(y_airfoil)
        pressure_contour(i,j) = (1-(yg*Ma^2*px(i+le,j)/ua))*pa;
        %pressure_contour(i,j) = (-2*px(i+le,j)/ua);
    end
end
end
% for i= 1 : te-le+1
%     x_airfoil(i) = 20 + (i-1)/(te-le);
%     pr(i) = (1+ (yg*Ma^2*px(i+le,1)/ua))*pa;
% end
%%  Plotting the results
 


% figure(3)
% plot (x_airfoil_fi,cp_store(1,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(2,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(3,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(4,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(5,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(6,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% xlabel('X-airfoil','FontSize',12);
% ylabel('Cp - Pressure coefficient','FontSize',12);
% title('Coefficient of Pressure at varying mach no','FontSize',14);
% % legend('M = 0.8','M = 0.82','M = 0.84','M = 0.86','M = 0.88','M = 0.9',...
% %    'FontSize',10)
% legend('Ma = 0.88')

% figure(4)
% %plot (x_airfoil_fi,Cdp_store(1,:),'LineWidth',2);
% plot (mach,Cdp_store,'LineWidth',2);
% %set(gca,'YDir','reverse');
% xlabel('Mach','FontSize',12);
% ylabel('Cdp - Coefficient of Pressure drag','FontSize',12);
% title('Coefficient of Pressure drag at varying mach no','FontSize',14);
% grid on; grid minor;

% figure(5)
% [x2,y2]= meshgrid(x_airfoil_fi,y_airfoil);
% contourf(x2,y2,pressure_contour.');
% xlabel('X-airfoil','FontSize',12);
% ylabel('Y','FontSize',12);
% zlabel('Pressure');
% title ('Pressure contour at M = 0.88','FontSize',14)

%surface velocity plot
figure(6)
plot (x_airfoil_fi,ue_storage(:,1),'LineWidth',2);
% hold on; 
% plot (x_airfoil_fi,ue_storage(:,2),'LineWidth',2);
% plot (x_airfoil_fi,ue_storage(:,3),'LineWidth',2);
% plot (x_airfoil_fi,ue_storage(:,4),'LineWidth',2);
grid on; grid minor;
xlabel('X-airfoil','FontSize',12);
ylabel('U_e - Surface velocity of airfoil','FontSize',12);
title('Surface velocity over the airfoil','FontSize',14);


% Boundary layer velocity profile plot
% figure(7)
% plot (u_vel(90,1:10),m(1:10),'LineWidth',2);
% hold on;
% plot (u_vel(100,1:10),m(1:10),'LineWidth',2);
% hold on;
% plot (u_vel(110,1:10),m(1:10),'LineWidth',2);
% grid on; grid minor;
% xlabel('y-coordinate','FontSize',12);
% ylabel('Velocity along the airfoil','FontSize',12);
% title('Velocity profile in the boundary layer','FontSize',14);
% legend('25% chord','50% chord','75% chord')

% Displacement thickness plot
figure(8)
plot (x_airfoil_fi,delta_star(1,:),'LineWidth',2);
hold on;
plot (x_airfoil_fi,delta_star(2,:),'LineWidth',2);
hold on;
plot (x_airfoil_fi,delta_star(3,:),'LineWidth',2);
 grid on; grid minor;
xlabel('X-airfoil','FontSize',12);
ylabel('Displacement thickness','FontSize',12);
title('Displacement thickness at varying mach no','FontSize',14);
legend('Ma = 0.8','Ma = 0.84','Ma = 0.86')
% 
% figure(9)
% plot (x_airfoil_fi,cf_store(1,:),'LineWidth',2);
% hold on;
% plot (x_airfoil_fi,cf_store(2,:),'LineWidth',2);
% hold on;
% plot (x_airfoil_fi,cf_store(3,:),'LineWidth',2);
%  grid on; grid minor;
% xlabel('X-airfoil','FontSize',12);
% ylabel('Skin Friction Coefficient','FontSize',12);
% title('Skin Friction Coefficient at varying mach no','FontSize',14);
% legend('Ma = 0.8','Ma = 0.84','Ma = 0.86')

%% Cp plot
% figure(10)
% plot (x_airfoil_fi,cp_store(1,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(2,:),'LineWidth',2);
% hold on;
% plot (x_airfoil_fi,cp_store(3,:),'LineWidth',2);
% grid on; grid minor;
% xlabel('X-airfoil','FontSize',12);
% ylabel('Coefficient of Pressure','FontSize',12);
% title('Coefficient of Pressure at varying mach no','FontSize',14);
% legend('Ma = 0.8','Ma = 0.84','Ma = 0.86')

%% Cp comparison between the original TSD and 
% figure(10)
% plot (x_airfoil_fi,cp_store(1,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil_fi,cp_store(2,:),'LineWidth',2);
% hold on;
% plot (x_airfoil_fi,cp_store(3,:),'LineWidth',2);
% plot (x_airfoil_fi,cp_TSD_store(1,:),'LineWidth',2);
% plot (x_airfoil_fi,cp_TSD_store(2,:),'LineWidth',2);
% plot (x_airfoil_fi,cp_TSD_store(3,:),'LineWidth',2);
% grid on; grid minor;
% xlabel('X-airfoil','FontSize',12);
% ylabel('Coefficient of Pressure','FontSize',12);
% title('Coefficient of Pressure at varying mach no','FontSize',14);
% legend('BLC coupled TSD - Ma = 0.8','BLC coupled TSD - Ma = 0.84',...
%     'BLC coupled TSD - Ma = 0.86','TSD - Ma = 0.8','TSD - Ma = 0.84',...
%     'TSD - Ma = 0.86');

figure(11)
plot (x_airfoil_fi,delta_star(1,:),'LineWidth',2);
hold on; 
plot (x_airfoil_fi,delta_star(2,:),'LineWidth',2);
plot (x_airfoil_fi,delta_star(3,:),'LineWidth',2);
plot (x_airfoil_fi,delta_star(4,:),'LineWidth',2);
grid on; grid minor;
xlabel('X-airfoil','FontSize',12);
ylabel('Displacement thickness','FontSize',12);
title('Displacement thickness (Iteration) at Ma - 0.8','FontSize',14);
legend('Iteration - 1','Iteration - 2','Iteration - 3','Iteration - 4');

%% Cd plots
figure(12)
plot (mach,Cd_store(:,1),'LineWidth',2);
%set(gca,'YDir','reverse');
hold on;
plot (mach,Cdf_store(:,1),'LineWidth',2);
hold on;
plot (mach,Cdp_store(:,1),'LineWidth',2);
grid on; grid minor;
xlabel('Mach no','FontSize',12);
ylabel('Drag Coefficients','FontSize',12);
title('Drag Coefficients at varying mach no','FontSize',14);
legend('Cd','Cdf','Cdp','FontSize',12)

%% Cl plots
figure(13)
plot (mach,Cl_store(:,1),'LineWidth',2);
%set(gca,'YDir','reverse');
hold on;
plot (mach,Clf_store(:,1),'LineWidth',2);
hold on;
plot (mach,Clp_store(:,1),'LineWidth',2);
grid on; grid minor;
xlabel('Mach no','FontSize',12);
ylabel('Lift Coefficients','FontSize',12);
title('Lift Coefficients at varying mach no','FontSize',14);
legend('Cl','Clf','Clp','FontSize',12)

%% Drag Polar

figure(14)
plot (Cd_store(:,1),Cl_store(:,1),'LineWidth',2);
grid on; grid minor;
xlabel('Cd - Drag Coefficient','FontSize',12);
ylabel('Cl - Lift Coefficient','FontSize',12);
title('Drag Polar- Cl vs Cd','FontSize',14);


function dydx = pos(x)
dydx = 0.08*((-4*x)+82);
end

function dy_a = pos_y(x,tc)
dy_a = tc*((-2*x^2)+(82*x)-840);
end

% plot(x_airfoil_co,cp_coarse,'linewidth',2)
% hold on;
% plot(x_airfoil_me,cp_medium,'linewidth',2)
% Unrecognized function or variable
% 'cp_medium'.
%  
% plot(x_airfoil_me,cp_med,'linewidth',2)
% set(gca,'YDir','reverse')
% plot(x_airfoil_fi,cp_fine,'linewidth',2)
% plot(x_airfoil_me,cp_med)
% set(gca,'YDir','reverse');
% hold on;
% plot(x_airfoil_fi,cp_fine)
% hold on;
% plot(x_airfoil_suf,cp_suf)
% legend('120*82 mesh grid','180*163 mesh grid','360*324 mesh grid')