close
clear
clc

%% Defining the computational domain
Lx = 1; Ly = 1;  % Overall grid size
Nx = 100; Ny = 100;  % total number of points
nx = Nx+1; ny = Ny+1;
dx = Lx/Nx; dy = Ly/Ny;
x = (0:Nx)*dx; y = (0:Ny)*dy;

eps = 1E-8;  % tolerance to achieve
u = zeros(nx,ny);

u(end,:) = 1; % Defining the boundary condition

%% Defining the variables and arrays to store error and iteration values 
u_jacobi = u; 
u_Gauss = u;
u_SOR = u;
u_old1 = u;
u_old2 = u;
u_old3 = u;

itr1 = 0;
itr2 = 0;
itr3 = 0;
error1 = 1;
error2 = 1;
error3 = 1;
error1_store = zeros(itr1,1); % to store CPU time for Jacobi
itr1_store = zeros(itr1,1); % to store Number of Iterations for Jacobi
cpu1 = zeros(itr1,1);  % to store CPU time for Jacobi
error2_store = zeros(itr2,1); % to store CPU time for Gauss Seidel
itr2_store = zeros(itr2,1); % to store Number of Iterations for Jacobi
cpu2 = zeros(itr2,1); % to store CPU time for Jacobi
error3_store = zeros(itr3,1); % to store CPU time for Jacobi
itr3_store = zeros(itr3,1); % to store Number of Iterations for Gauss Seidel
cpu3 = zeros(itr3,1); % to store CPU time for SOR
%% Jacobi 
tic
while (error1>eps)
    
    itr1 = itr1+1;
    for i = 2:Nx
        for j = 2:Ny
            u_jacobi(i,j) = 0.25*(u_old1(i+1,j)+u_old1(i-1,j)+u_old1(i,j+1)+u_old1(i,j-1));
        end
        
    end

    error1 = sqrt(sum(sum(abs(u_jacobi-u_old1).^2)));
    error1_store(itr1,1) = error1; 
    itr1_store(itr1,1) = itr1;
    cpu1(itr1,1) = toc;
    u_old1 = u_jacobi;
    fprintf('%d %d\n',itr1,error1);
end

%% Gauss-Seidel
tic
while (error2>eps)
    itr2 = itr2+1;
    for i = 2:Nx
        for j = 2:Ny
            u_Gauss(i,j) = 0.25*(u_old2(i+1,j)+u_Gauss(i-1,j)+u_old2(i,j+1)+u_Gauss(i,j-1));

        end
        
    end
    error2 = sqrt(sum(sum(abs(u_Gauss-u_old2).^2)));
    error2_store(itr2,1) = error2; 
    itr2_store(itr2,1) = itr2;
    cpu2(itr2,1) = toc;
    u_old2 = u_Gauss;
    fprintf('%d %d\n',itr2,error2);
end


%% SOR
tic
h = 1.5;
while (error3>eps)
    itr3 = itr3+1;
    for i = 2:Nx
        for j = 2:Ny
            u_SOR(i,j) = (1-h)*u_old3(i,j)+ (h/4)*(u_old3(i+1,j)+u_SOR(i-1,j)+u_old3(i,j+1)+u_SOR(i,j-1));
        end
        
    end
    error3 = sqrt(sum(sum(abs(u_SOR-u_old3).^2)));
    error3_store(itr3,1) = error3; 
    itr3_store(itr3,1) = itr3;
    cpu3(itr3,1) = toc;
    u_old3 = u_SOR;
    fprintf('%d %d\n',itr3,error3);
end


% Plotting the contour
[X,Y] = meshgrid (x,y);
v = [0.8 0.6 0.4 0.2 0.1 0.05 0.01];
figure(1);
contourf(X,Y,u_jacobi); 
colorbar;
axis equal;
xlabel('$X$','Interpreter','latex','FontSize',12);
ylabel('$Y$','Interpreter','latex','FontSize',12);
title('Solution of the Laplace equation','Interpreter','latex','FontSize',14);

% Plotting the error vs iterations
figure(2);
plot(itr1_store,log10(error1_store));
hold on;
plot(itr2_store,log10(error2_store))
hold on; 
plot(itr3_store,log10(error3_store));
xlabel('$Iterations$','Interpreter','latex','FontSize',12);
ylabel('$Log(error)$','Interpreter','latex','FontSize',12);
title('Error vs Number of Interations','Interpreter','latex','FontSize',14);
legend('Jacobi','Gauss-seidel','SOR','Interpreter','latex','FontSize',12);

% Plotting the error vs total CPU time
figure(3);
plot(cpu1,log10(error1_store));
hold on; 
plot(cpu2,log10(error2_store));
hold on; 
plot(cpu3,log10(error3_store));
xlabel('$CPU time$','Interpreter','latex','FontSize',12);
ylabel('$Log(error)$','Interpreter','latex','FontSize',12);
title('Error vs Total CPU time','Interpreter','latex','FontSize',14);
legend('Jacobi','Gauss-seidel','SOR','Interpreter','latex','FontSize',12);
