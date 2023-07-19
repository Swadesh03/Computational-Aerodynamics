# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 08:49:34 2023

@author: NOTEBOOK
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the grid parameters
x_min = 0
x_max = 40
dx = 1
x = np.arange(x_min, x_max + dx, dx)

# Analytical solution
u_analytical = 0.5 * (1 + np.tanh(250 * ((x - 20) - 0.5 * 10)))
plt.plot(x, u_analytical, '-s', linewidth=2, label='Analytical Sol')

# Upwind 
# Initial condition
L = 40  # Length of domain
T = 10  # Final computational time
Nx = 41  # Number of spatial grid points
dx = L / (Nx - 1)  # Spatial grid size
c = 0.5  # Advection velocity
dt = 0.5  # Time step
Nt = int(np.ceil(T / dt))  # Number of time steps



# Solution matrix
u = np.zeros((Nx+1, Nt + 1))
x = np.arange(x_min, Nx+1, dx)
u[:, 0] = 0.5 * (1 + np.tanh(250 * (x - 20)))


# Dirichlet boundary condition
u[0, :] = u[0, 0]  # Set the boundary condition for the first row
u[-1, :] = u[-1, 0]  # Set the boundary condition for the last row


# Time marching using upwind scheme
for n in range(Nt):
    for i in range(1, Nx):
        u[i, n + 1] = u[i, n] - c * dt / dx * (u[i, n] - u[i - 1, n])

# Plot the numerical solution at the final time for upwind scheme
plt.plot(x, u[:, -1], '-s', linewidth=2, label='Upwind Scheme')

# Time marching using Lax scheme

for n in range(Nt):
    for i in range(1, Nx):
        u[i, n + 1] = 0.5 * (u[i + 1, n] + u[i - 1, n]) - (c / 2) * dt / dx * (u[i + 1, n] - u[i - 1, n])

# Plot the numerical solution at the final time for Lax scheme
plt.plot(x, u[:, -1], '-+', linewidth=2, label='Lax Scheme')

# Time marching using Leap Frog scheme

for n in range(0, Nt-1):
    for i in range(0, Nx):
        u[i, n + 1] = u[i, n - 1] - (c * dt / dx) * (u[i + 1, n] - u[i - 1, n])

# Plot the numerical solution at the final time for Leap Frog scheme
plt.plot(x, u[:, -1], '-x', linewidth=2, label='Leap Frog Scheme')

# Time marching using Lax-Wendroff scheme
u = np.zeros((Nx+1, Nt + 1))
x = np.arange(x_min, Nx+1, dx)
u[:, 0] = 0.5 * (1 + np.tanh(250 * (x - 20)))

# Dirichlet boundary condition
u[0, :] = u[0, 0]  # Set the boundary condition for the first row
u[-1, :] = u[-1, 0]  # Set the boundary condition for the last row

for n in range(Nt):
    for i in range(1, Nx):
        u[i, n + 1] = u[i, n] - (c * dt / (2 * dx)) * (u[i + 1, n] - u[i - 1, n]) + (
                c ** 2 * dt ** 2 / (2 * dx ** 2)) * (u[i + 1, n] - 2 * u[i, n] + u[i - 1, n])

# Plot the numerical solution at the final time for Lax-Wendroff scheme
plt.plot(x, u[:, -1], '-o', linewidth=2, label='Lax-Wendroff Scheme')

# Time marching using MacCormack scheme
u = np.zeros((Nx+1, Nt + 1))
x = np.arange(x_min, Nx+1, dx)
u[:, 0] = 0.5 * (1 + np.tanh(250 * (x - 20)))

# Dirichlet boundary condition
u[0, :] = u[0, 0]  # Set the boundary condition for the first row
u[-1, :] = u[-1, 0]  # Set the boundary condition for the last row
up = np.copy(u)
for n in range(Nt):
    for i in range(1, Nx):
        up[i, n + 1] = u[i, n] - c * dt / dx * (u[i + 1, n] - u[i, n])
        u[i, n + 1] = 0.5 * (up[i, n + 1] + u[i, n] - c * dt / dx * (up[i, n + 1] - up[i - 1, n + 1]))

# Plot the numerical solution at the final time for MacCormack scheme
plt.plot(x, u[:, -1], '-d', linewidth=2, label='MacCormack Scheme')

# Plot the analytical solution
#plt.plot(x, u_analytical, label='Analytical Solution', linewidth=2)

plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.title('1D linear advection equation (comparison between analytical and all schemes at dt = 0.5)',
          fontsize=14)
plt.legend()
plt.grid()
plt.show()
