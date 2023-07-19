# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:39:42 2023

@author: NOTEBOOK
"""

import numpy as np
import matplotlib.pyplot as plt
import time

# Defining the computational domain
Lx = 1
Ly = 1
Nx = 100
Ny = 100
dx = Lx / Nx
dy = Ly / Ny
x = np.arange(0, Lx + dx, dx)
y = np.arange(0, Ly + dy, dy)

eps = 1E-4  # tolerance to achieve
u = np.zeros((Nx + 1, Ny + 1))

u[-1, :] = 1  # Defining the boundary condition

# Defining the variables and arrays to store error and iteration values
u_jacobi = np.copy(u)
u_Gauss = np.copy(u)
u_SOR = np.copy(u)
u_old1 = np.copy(u)
u_old2 = np.copy(u)
u_old3 = np.copy(u)

itr1 = 0
itr2 = 0
itr3 = 0
error1 = 1
error2 = 1
error3 = 1
error1_store = []
itr1_store = []
cpu1 = []
error2_store = []
itr2_store = []
cpu2 = []
error3_store = []
itr3_store = []
cpu3 = []

# Jacobi
start_time = time.time()
while error1 > eps:
   
    itr1 += 1
    for i in range(1, Nx):
        for j in range(1, Ny):
            u_jacobi[i, j] = 0.25 * (u_old1[i + 1, j] + u_old1[i - 1, j] + u_old1[i, j + 1] + u_old1[i, j - 1])

    error1 = np.sqrt(np.sum(np.abs(u_jacobi - u_old1)**2))
    error1_store.append(error1)
    itr1_store.append(itr1)
    
    end_time = time.time()
    cpu_time_taken = end_time - start_time
    cpu1.append(cpu_time_taken)
    u_old1 = np.copy(u_jacobi)
    print(itr1, error1)

# Gauss-Seidel
start_time = time.time()
while error2 > eps:
    #start_time = time.time()
    itr2 += 1
    for i in range(1, Nx):
        for j in range(1, Ny):
            u_Gauss[i, j] = 0.25 * (u_old2[i + 1, j] + u_Gauss[i - 1, j] + u_old2[i, j + 1] + u_Gauss[i, j - 1])

    error2 = np.sqrt(np.sum(np.abs(u_Gauss - u_old2)**2))
    error2_store.append(error2)
    itr2_store.append(itr2)
    
    end_time = time.time()
    cpu_time_taken = end_time - start_time
    cpu2.append(cpu_time_taken)
    #cpu2.append(time.time())
    u_old2 = np.copy(u_Gauss)
    print(itr2, error2)

# SOR
h = 1.5
start_time = time.time()
while error3 > eps:
    #start_time = time.time()
    itr3 += 1
    for i in range(1, Nx):
        for j in range(1, Ny):
            u_SOR[i, j] = (1 - h) * u_old3[i, j] + (h / 4) * (u_old3[i + 1, j] + u_SOR[i - 1, j] + u_old3[i, j + 1] + u_SOR[i, j - 1])

    error3 = np.sqrt(np.sum(np.abs(u_SOR - u_old3)**2))
    error3_store.append(error3)
    itr3_store.append(itr3)
    
    end_time = time.time()
    cpu_time_taken = end_time - start_time
    cpu3.append(cpu_time_taken)
    #cpu3.append(time.time())
    u_old3 = np.copy(u_SOR)
    print(itr3, error3)

# Plotting the contour
X, Y = np.meshgrid(x, y)
plt.figure(1)
plt.contourf(X, Y, u_jacobi)
plt.colorbar()
plt.axis('equal')
plt.xlabel('$X$', fontsize=12)
plt.ylabel('$Y$', fontsize=12)
plt.title('Solution of the Laplace equation', fontsize=14)

# Plotting the error vs iterations
plt.figure(2)
plt.plot(itr1_store, np.log10(error1_store))
#plt.hold(True)
plt.plot(itr2_store, np.log10(error2_store))
#plt.hold(True)
plt.plot(itr3_store, np.log10(error3_store))
plt.xlabel('$Iterations$', fontsize=12)
plt.ylabel('$Log(error)$', fontsize=12)
plt.title('Error vs Number of Iterations', fontsize=14)
plt.legend(['Jacobi', 'Gauss-seidel', 'SOR'], fontsize=12)

# Plotting the error vs total CPU time
plt.figure(3)
plt.plot(cpu1, np.log10(error1_store))
#plt.hold(True)
plt.plot(cpu2, np.log10(error2_store))
#plt.hold(True)
plt.plot(cpu3, np.log10(error3_store))
plt.xlabel('$CPU time$', fontsize=12)
plt.ylabel('$Log(error)$', fontsize=12)
plt.title('Error vs Total CPU time', fontsize=14)
plt.legend(['Jacobi', 'Gauss-seidel', 'SOR'], fontsize=12)

plt.show()
