close all
clear 
clc

load NACA0012_flowfieldv2.dat

imax = 513;
jmax = 257;
TE_start = 65; % Trailing Edge Lower Point
TE_end = 449; % Trailing Edge Upper Point
LE = 257; % Leading Edge Point
gamma = 1.4;
mach = 0.1;
p0 = 1.; % Freestream static pressure
alpha = 8; % Angle of Attack
Re = 3e6; % Reynolds number



k = 1;
for j=1:jmax
    for i=1:imax
        x(i,j) = NACA0012_flowfieldv2(k,1);
        y(i,j) = NACA0012_flowfieldv2(k,2);
        k = k +1;
    end
end



% Plot the Airfoil Surface
figure(1)
plot(x(TE_start:TE_end,1),...
    y(TE_start:TE_end,1),'k*-')
axis equal
title('Airfoil','FontSize',14)



% Plot the Computational Domain.
figure(2)
axis equal
hold on
for j=1:jmax
    plot(x(1:imax,j),y(1:imax,j),'k-')  
end
for i=1:imax
    plot(x(i,1:jmax),y(i,1:jmax),'k-')
end
title('Computational Domain','FontSize',14)

% 
% % Plot the Computational Domain.
% figure(3)
% axis equal
% hold on
% for j=1:jmax
%     plot(x(1:imax,j),y(1:imax,j),'k-')  
% end
% for i=1:imax
%     plot(x(i,1:jmax),y(i,1:jmax),'k-')
% end
% axis([-.1 1.1 -.2 .2])


Nvariables = 7;
% Transfer all other flow properties and evaluate Pressure
k = 1;
for j=1:jmax
    for i=1:imax
        for n=1:Nvariables
            w(i,j,n) = NACA0012_flowfieldv2(k,2+n);
        end
        k = k +1;
    end
end




% ***********************************%
% Plot Surface Pressure
% ***********************************%
for j=1:jmax
    for i=1:imax
        pressure(i,j) = (gamma -1)*(w(i,j,4) ...
            -.5*(w(i,j,2)^2  +w(i,j,3)^2)/w(i,j,1));
        velocity(i,j) = sqrt((w(i,j,2)/w(i,j,1))^2 ...
                            +(w(i,j,3)/w(i,j,1))^2);
    end 
end


% Evaluate Shear Stress.
% for i=TE_start:TE_end
%     dxi = x(i,1) -x(i-1,1); % dx/dxi
%     dyi = y(i,1) -y(i-1,1); % dy/dxi
%     dxj = x(i,2) -x(i,1); % dx/deta
%     dyj = y(i,2) -y(i,1); % dy/deta
%     dsj = 1/(dxi*dyj -dyi*dxj);
% 
%     dui = w(i,1,2)/w(i,1,1) -w(i-1,1,2)/w(i-1,1,1); % du/dxi
%     dvi = w(i,1,3)/w(i,1,1) -w(i-1,1,3)/w(i-1,1,1); % dv/dxi 
%     duj = w(i,2,2)/w(i,2,1) -w(i,1,2)/w(i,1,1); % du/deta
%     dvj = w(i,2,3)/w(i,2,1) -w(i,1,3)/w(i,1,1); % dv/deta
% 
%     dux = (dui*dyj -duj*dyi)*dsj; % du/dx
%     dvx = (dvi*dyj -dvj*dyi)*dsj; % dv/dx
%     duy = (duj*dxi -dui*dxj)*dsj; % du/dy
%     dvy = (dvj*dxi -dvi*dxj)*dsj; % dv/dy
% end

%% Question 1 
Cp = ((pressure((65:449),1)/p0)-1)/(0.5*gamma*(mach)^2);
figure(4)
plot(x(257:449),Cp(193:385),'LineWidth',1)
hold on;
plot(x(65:257),Cp(1:193),'LineWidth',1)
set(gca,'YDir','reverse');
ylabel('Coefficient of Pressure')
xlabel('X-airfoil')
grid on;
grid minor;
title('Coeffcient of Pressure variation over airfoil')
legend('Upper Surface', 'Lower Surface')

%% Question 2 
% mew = NACA0012_flowfieldv2(TE_start:TE_end,7)+...
%     NACA0012_flowfieldv2(TE_start:TE_end,8);
% Txx = (mew*2*dux)-((2/3)*mew*(dux+dvy));
% Txy = mew*(duy+dvx);
% Tyy = (mew*2*dvy)-((2/3)*mew*(dux+dvy));
Cf = zeros(TE_end-TE_start,1);
Tw = zeros(TE_end-TE_start,1);
for i = TE_start:TE_end
    dyp = 0.5*(y(i+1,1) +y(i,1));
    dym = 0.5*(y(i,1) +y(i-1,1));
    dxp = 0.5*(x(i+1,1) +x(i,1));
    dxm = 0.5*(x(i,1) +x(i-1,1));

    dx = dxp - dxm;
    dy = dyp - dym;
    dA = sqrt((dx)^2+(dy)^2);
    cos_a = dx/dA; sin_a = dy/dA;

    dxi = x(i,1) -x(i-1,1); % dx/dxi
    dyi = y(i,1) -y(i-1,1); % dy/dxi
    dxj = x(i,2) -x(i,1); % dx/deta
    dyj = y(i,2) -y(i,1); % dy/deta
    dsj = 1/(dxi*dyj -dyi*dxj);

    dui = w(i,1,2)/w(i,1,1) -w(i-1,1,2)/w(i-1,1,1); % du/dxi
    dvi = w(i,1,3)/w(i,1,1) -w(i-1,1,3)/w(i-1,1,1); % dv/dxi 
    duj = w(i,2,2)/w(i,2,1) -w(i,1,2)/w(i,1,1); % du/deta
    dvj = w(i,2,3)/w(i,2,1) -w(i,1,3)/w(i,1,1); % dv/deta

    dux = (dui*dyj -duj*dyi)*dsj; % du/dx
    dvx = (dvi*dyj -dvj*dyi)*dsj; % dv/dx
    duy = (duj*dxi -dui*dxj)*dsj; % du/dy
    dvy = (dvj*dxi -dvi*dxj)*dsj; % dv/dy

    mew = NACA0012_flowfieldv2(i,7) + NACA0012_flowfieldv2(i,8);

    Txx = (mew*2*dux)-((2/3)*mew*(dux+dvy));
    Txy = mew*(duy+dvx);
    Tyy = (mew*2*dvy)-((2/3)*mew*(dux+dvy));

    Tx = -Txx*sin_a + Txy*cos_a;
    Ty = -Txy*sin_a + Tyy*cos_a;
    Tw(i-TE_start+1) = Tx*cos_a + Ty*sin_a;
    Cf(i-TE_start+1) = Tw(i-TE_start+1)/(0.5*p0*(mach)^2*gamma);
end
figure(5)
plot(x(65:256),abs(Cf(1:192)),'LineWidth',1)
hold on; 
plot(x(257:449),Cf(193:385),'LineWidth',1)
ylabel('Skin friction Coefficient')
xlabel('X-airfoil')
title('Skin friction coefficient variation over airfoil')
grid on;
grid minor;
legend('Lower Skin friction coefficient ',...
    'Upper Skin friction coefficient')


%% Question 3
Cpx = 0;
Cpy = 0;
Cfx = 0;
Cfy = 0;
alpha_r = alpha*pi/180;
for i = TE_start:TE_end
    dyp = 0.5*(y(i+1,1) +y(i,1));
    dym = 0.5*(y(i,1) +y(i-1,1));
    dxp = 0.5*(x(i+1,1) +x(i,1));
    dxm = 0.5*(x(i,1) +x(i-1,1));

    dx = dxp - dxm;
    dy = dyp - dym;
    dA = sqrt((dx)^2+(dy)^2);
    cos_a = dx/dA; sin_a = dy/dA;
    n1 = -sin_a; n2 = cos_a;
    t1 = cos_a; t2 = sin_a;

    Cpi = ((pressure(i,1)/p0)-1)/(0.5*gamma*(mach)^2);
    Cpx = Cpx - (Cpi*n1*dA);
    Cpy = Cpy - (Cpi*n2*dA);
    Clp = -Cpx*sin(alpha_r) + Cpy*cos(alpha_r);
    Cdp = Cpx*cos(alpha_r) + Cpy*sin(alpha_r);
    
    
    Cfx = Cfx + Cf(i-TE_start+1,1)*t1*dA;
    Cfy = Cfy + Cf(i-TE_start+1,1)*t2*dA;
    Clf = -Cfx*sin(alpha_r) + Cfy*cos(alpha_r);
    Cdf = Cfx*cos(alpha_r) + Cfy*sin(alpha_r); 

    Cl = Clf + Clp;
    Cd = Cdf + Cdp;
end
cl_actual = ((1.12074 - 0.870758)/(10.1891-7.96346))*(8-7.96346) +  0.870758;
cd_actual = 0.00915707;
error_cl = ((cl_actual-Cl)/cl_actual)*100;
error_cd = ((cd_actual-Cd)/cd_actual)*100;
%% Question 4 
% ui_plusx1 = zeros(5,jmax);
% % ui_plusx2 = zeros(5,jmax);
% % ui_plusx3 = zeros(5,jmax);
% % ui_plusx4 = zeros(5,jmax);
% % ui_plusx5 = zeros(5,jmax);
% yi_plus = zeros(5,jmax);
% k = 1;
% for i = [296,304,351,447,408]
% for j = 1:jmax
% %     du_vel = cos_a*NACA0012_flowfieldv2(i,4)/NACA0012_flowfieldv2(i,3);
% %     dv_vel = sin_a*NACA0012_flowfieldv2(i,5)/NACA0012_flowfieldv2(i,3);
% %     vj = sqrt((du_vel)^2+(dv_vel)^2);
%     u_friction = sqrt(abs(Tw((i-TE_start+1)))/w(i,1,1));
%     ui_plusx1(k,j) = velocity(i,j)/u_friction;
%     %mew = NACA0012_flowfieldv2(i,7) + NACA0012_flowfieldv2(i,8);
% 
%     mew =w(i,1,5) + w(i,1,6);
%     vw = mew/w(i,1,1);
%     y_normal = sqrt(w(i,j,7));
%     yi_plus(k,j) = y_normal*u_friction/vw;
% %     ui_plus1 = (1/0.32)*log(yi_plus)-0.8;
% end
% k = k+1;
% end
% yi_plusx = ui_plusx1(1:27);
% figure(6)
% semilogx(yi_plus(1,:),ui_plusx1(1,:),'LineWidth',1)
% axis([1 1000000 0 120])
% grid on;
% grid minor;
% hold on; 
% semilogx(yi_plus(2,:),ui_plusx1(2,:),'LineWidth',1)
% semilogx(yi_plus(3,:),ui_plusx1(3,:),'LineWidth',1)
% semilogx(yi_plus(4,:),ui_plusx1(4,:),'LineWidth',1)
% semilogx(yi_plus(5,:),ui_plusx1(5,:),'LineWidth',1)
% 
% ylabel('u+')
% xlabel('y+')
% title('Plot of u+ vs y+ in laminar and turbulent zones')
% legend('x = 0.029', 'x = 0.048','x = 0.394','x = 0.993','x = 1')

%% Q4
ui_plusx1 = zeros(2,jmax);
yi_plus = zeros(2,jmax);
k = 1;
for i = [276,313]
for j = 1:jmax
    u_friction = sqrt(abs(Tw((i-TE_start+1)))/w(i,1,1));
    ui_plusx1(k,j) = velocity(i,j)/u_friction;

    mew =w(i,1,5) + w(i,1,6);
    vw = mew/w(i,1,1);
    y_normal = sqrt(w(i,j,7));
    yi_plus(k,j) = y_normal*u_friction/vw;
%     ui_plus1 = (1/0.32)*log(yi_plus)-0.8;
end
k = k+1;
end
yi_plusx = ui_plusx1(1:27);
figure(6)
semilogx(yi_plus(1,:),ui_plusx1(1,:),'LineWidth',1)
axis([1 1000000 0 40])
grid on;
grid minor;
hold on; 
semilogx(yi_plus(2,:),ui_plusx1(2,:),'LineWidth',1)
ylabel('u+')
xlabel('y+')
title('Plot of u+ vs y+ in laminar and turbulent zones')
legend('x = 0.005 (laminar region)', 'x = 0.08 (turbulent region)');

%% question 5
ui_plus = zeros(jmax,1);
ui_plus1 = zeros(jmax,1);
yi_plus = zeros(jmax,1);
for j = 1:jmax
    i = 363;
    u_friction = sqrt(abs(Tw((i-TE_start+1)))/w(i,1,1));
    ui_plus(j) = velocity(i,j)/u_friction;

    mew =w(i,1,5) + w(i,1,6);
    vw = mew/w(i,1,1);
    y_normal = sqrt(w(i,j,7));
    yi_plus(j) = y_normal*u_friction/vw;
    ui_plus1 = (1/0.32)*log(yi_plus)-0.8;
end
yi_plus1 = ui_plus(1:27);

figure(7)
semilogx(yi_plus,ui_plus,'LineWidth',1)
axis([1 10000 0 30])
grid on;
grid minor;
hold on; 
semilogx(yi_plus1,ui_plus(1:27),'LineStyle','--','LineWidth',1);
hold on; 
semilogx(yi_plus,ui_plus1,'LineStyle','-.','LineWidth',1);

ylabel('u+')
xlabel('y+')
title('Zones in Turbulent boundary layer')
legend('Velocity profile, Re = 2000000', 'u+ = y+',...
    'u+ = (1/0.32) * ln(y+) - 0.8')

%% Question 6 
k = 1;
v = zeros(8,jmax);
y_ver = zeros(4,jmax);
%for i = [336,363,449,485]
for i = [22,29,36,42,472,478,485,492]
    for j = 1:jmax
        v(k,j) = velocity(i,j);
        y_ver(k,j) = y(i,j)-y(i,1);
    end
    k = k+1;
end
figure(8)
plot(v(1,:),y_ver(1,:),'linewidth',1)
hold on; 
plot(v(2,:),y_ver(2,:),'linewidth',1)
plot(v(3,:),y_ver(3,:),'linewidth',1)
plot(v(4,:),y_ver(4,:),'linewidth',1)
plot(v(5,:),y_ver(5,:),'linewidth',1)
plot(v(6,:),y_ver(6,:),'linewidth',1)
plot(v(7,:),y_ver(7,:),'linewidth',1)
plot(v(8,:),y_ver(8,:),'linewidth',1)
axis([0.075 0.125 -0.6 0.6])
grid on;
grid minor;
ylabel('Distance from centreline')
xlabel('Velocity')
title('Velocity profile through the wake regions')
legend('1/4 chord - below centerline', '1/2 chord- below centerline',...
    '1 chord- below centerline','2 chord- below centerline',...
    '1/4 chord - above centerline','1/2 chord - above centerline',...
    '1 chord - above centeline','2 chord - above centerline')