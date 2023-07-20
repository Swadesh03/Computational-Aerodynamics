% close
% clear
% clc

%% Defining the mesh
%Ma = 0.88;   % Mach no of the flow
L = 50;
xf = 2;% Length of domain
n1 = 40; % number of grid points before the leading edge
n2 = 20; % number of grid points for the airfoil
n3 = 60;  % number of grid points after the the trailing edge
le = n1;
te = n1+n2;
tc = 0.08;
nx = n1+n2+n3;
n4 = 82;
n = zeros(nx,1);
m = zeros(n4,1);
m(1) = 0;
m(2) = tc/2;
n(n1) = n1/2;
dx = 1/n2;
for i = n1+1:n1+n2
    n(i) = n(i-1) + dx; 
end
for i = n1+n2+1:nx
    n(i) = n(i-1) + (n(i-1)-n(i-2))*1.05765^1.05;
end
for i = n1:-1:2
    n(i-1) = n(i) - (n(i+1)-n(i))*1.0872^1.1;
end
% nt1 =zeros(nx,1);
% nt1(n1+n2+1:nx) = n(n1+n2+1:nx) - n1-1; nt2 = -flip(nt1);
% n(1:n1) = n1;
% n(1:n1) = n(1:n1)+(nt2(n1/2+1:n1+n2));
for i = 3:n4
    m(i) = m(i-1) + (m(i-1)-m(i-2))*1.04847^1.1;
end
% n = n/xf;

% n = n*L/n(nx,1);
n(1)=0;
m(end) = 50; n(end) = 50;
n_medium = n;
[x1,y1]= meshgrid(n,m);
x = x1'; y = y1';
% plot(x1,y1,'k-');
% hold on;
% plot(x,y,'k-')
% % axis equal;
% xlabel('X','FontSize',12);
% ylabel('Y','FontSize',12);
% title('Polynomial Grid stretching', 'FontSize',14);



% mach = (0.8:0.02:0.9);
% cp_store = zeros(length(mach),11);
% for k = 1:length(mach)
% Defining the parameters
yg = 1.4; % gamma for air
nx = n1+n2+n3;
ny = n4;
Ma = 0.88;   % Mach no of the flow
Ta = 293;
pa = 100;
R = 287.058;
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

eps = 1E-4;
error = 1;   % initial error value
itr = 0;   % iteration
error_store = zeros(itr,1);   % storing the errors
itr_store = zeros(itr,1);    % storing the iteration
%% Boundary condition
for i = 2 : nx-1
    if i <= le
        p(i,2) = p(i,1);
    elseif i > te
        p(i,2) = p(i,1);
    else
         %x_pos = n1 + (i-le)/(te-le);
         x_pos = n(i);
         p(i,1)=p(2,1)-ua*pos(x_pos)*(m(2)-m(1));
    end
end

px = zeros(nx,ny);
for i = 2 : nx-1
    for j = 1: ny
        px(i,j) = (p(i+1,j)-p(i-1,j))/((n(i+1)-n(i-1)));
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
    end
end
p_gauss = p;

while error>eps
    itr = itr+1;
    for i = 3:nx-1
        for j = 2:ny-1
            p_gauss(i,j) = (-c(i,j)*p_gauss(i,j-1)-d(i,j)*p_gauss(i-1,j)...
                -e(i,j)*p(i+1,j)-b(i,j)*p(i,j+1)-g(i,j)*p_gauss(i-2,j))...
                /a(i,j);
        end
    end
    error = max(abs(p_gauss-p),[],"all");
    error_store(itr,1) = error; 
    itr_store(itr,1) = itr;
    fprintf('%d %d\n',itr,error);
    for i = 2 : nx-1
        if i <= le
            p_gauss(i,1) = p_gauss(i,2);
        elseif i > te
            p_gauss(i,1) = p_gauss(i,2);
        else
            %x_pos = n1 + (i-le)/(te-le);
            x_pos = n(i);
            p_gauss(i,1)=p_gauss(i,2)-ua*pos(x_pos)*(m(2)-m(1));
        end
    end
    p = p_gauss;
    % Boundary condition update
    for i = 2 : nx-1
        for j = 1: ny
            px(i,j)= (p_gauss(i+1,j)-p_gauss(i-1,j))/((n(i+1)-n(i-1)));
        end
    end
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
        end
    end
end
cp_med = zeros(1,te-le);
pr = zeros(1,te-le);
x_airfoil_me = zeros(1,te-le);
for i= 1 : te-le+1
    x_airfoil_me(i) = n1/2 + (i-1)/(te-le);
    cp_med(i) = -2*px(i+le,1)/ua;
end
%cp_store(k,:)= cp_med;
r = 1;
while m(r) <= 1
    r = r+1;
end
% 
y_airfoil = m(1:r);

%% Plotting the results
figure(1);
plot(itr_store,log10(error_store),'LineWidth',2);
xlabel('Iterations','FontSize',12);
ylabel('Log(error)','FontSize',12);
title('Error vs Number of Interations at M=0.8', 'FontSize',14);

figure(2);
plot (x_airfoil_me,cp_med,'LineWidth',2);
set(gca,'YDir','reverse');

function dydx = pos(x)
dydx = 0.08*((-4*x)+82);
end
% figure(3)
% [x2,y2]= meshgrid(x_airfoil_me,y_airfoil);
% contourf(x2,y2,pressure_contour.');
% xlabel('X-airfoil','FontSize',12);
% ylabel('Y','FontSize',12);
% % zlabel('Pressure');
% title ('Pressure contour at M = 0.8','FontSize',14)



% %pressure_contour = zeros(size(x_airfoil,2),size(y_airfoil,2));
% pressure_contour = zeros(length(x_airfoil_me),length(y_airfoil));
%  
% for i = 1 : length(x_airfoil_me)
%     for j = 1: length(y_airfoil)
%         pressure_contour(i,j) = (1- (yg*Ma^2*px(i+le,j)/ua))*pa;
%         %pressure_contour(i,j) = (-2*px(i+le,j)/ua);
%     end
% end

% for i= 1 : te-le+1
%     x_airfoil(i) = 20 + (i-1)/(te-le);
%     pr(i) = (1+ (yg*Ma^2*px(i+le,1)/ua))*pa;
% end


% hold on;
% plot (x_airfoil,cp_store(2,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil,cp_store(3,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil,cp_store(4,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil,cp_store(5,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% plot (x_airfoil,cp_store(6,:),'LineWidth',2);
% set(gca,'YDir','reverse');
% hold on;
% xlabel('X-airfoil','FontSize',12);
% ylabel('Cp - Pressure coefficient','FontSize',12);
% title('Coefficient of Pressure at M = 0.8','FontSize',14);
% %legend('M = 0.8','M = 0.82','M = 0.84','M = 0.86','M = 0.88','M = 0.9',...
%    % 'FontSize',10)
% 

% for i = n1:-1:1
%     n(i) = n(i+1) - (n(n1+n2+j)-n(n1+n2+j-1));
% end
% n(n1+1:n1+n2) = 
% x1max = n1;
% 
% % x1r = 1.20568;
% % x3r = 1.12944;
% % yr = 1.13747;
% xmin = n1;
% xmax = n1+1;
% x3min = n1+1;
% x = rand(nx,n4);
% y = rand(nx,n4);
% x1 = zeros(n1,1);
% x3 = zeros(n3,1);
% y1 = zeros(n4,1);
% x2 = linspace(xmin,xmax,n2)';
% %dx = (xmax-xmin)/n2;
% x1min = dx;
% for i = n1:-1:1
%     if i == n1 
%     x1(i) = x2(n1-(n1-1),1)-dx;
%     elseif i == n1-1
%         x1(i) = x1(i+1) - 2*dx;
%     else
%         x1(i)= x(i+1)+(x(i+1)-x(i+2))*1.1^1.25;
%     end
% %     x1(i) = x1max-0.1*(i-1)^1.2;
%     x1max = x1(i);
% end
% % x1 = x1*n1/x1(1,end);
% %x2 = linspace(xmin,xmax,n2);
% for i = 1:n3
%     x3(i) = x3min+0.1*x3r^(i-1);
% %     x3(i) = x3min+0.1*(i-1)^1.2;
%     x3min = x3(i);
% end
% % x3 = x3*n3/x3(1,end);
% ymin = 0;
% for i = 1:n4
%     if i == 1
%         y1(i)=ymin;
%     else
% %     y(i) = ymin + 0.04*(i-1)^2;
%     y1(i) = ymin + 0.04*yr^(i-1);
%     ymin = y1(i);
%     end
% end
% y2 = y1*L/y1(1,n4);
% % make y vector

% x1 = [flip(x1),x2,x3];
% [x2,y4]= meshgrid (x1,y2);
% a = x2'; b = y4';
%plot(x2,y4,'k-')