/****
Author - Swadesh Suman
Code - Navier Stokes 2D Channel Flow
Time - forward stepping
Spatial - backward discretization - first order, Centered discretization - second order
Location - McGill University, Canada
****/

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <numeric>

// Constant variable declaration ---
const int nx = 41;  // no of points in x
const int ny = 41;  // no of points in y
const int nt = 100; // no of time steps
const int nit = 50; // no of iteration to find the optimum value of pressure
const double dx = 2.0 / (nx - 1);   // calculation of dx
const double dy = 2.0 / (ny - 1);   // calculation of dy
const double rho = 1.0; // density of fluid
const double nu = 0.1;  // viscosity of fluid
const double F = 1.0;   // source term in the x-momentum eqn to mimic the effect of pressure-driven flow
const double dt = 0.001; // each time step length - 0.001s

/*  function -- build_up_b
    return b vector -- which can be later used in poisson equation
    b is just a part of poisson equation, seperated from poisson eqn
    to keep it short.
    input - b, rho, dt, u, v, dx, dy
    output - b vector
*/

std::vector<std::vector<double>> build_up_b(std::vector<std::vector<double>> b,
                                           double rho, double dt,
                                           const std::vector<std::vector<double>>& u,
                                           const std::vector<std::vector<double>>& v,
                                           double dx, double dy) {
    std::vector<std::vector<double>> result = b;
    // updating the b
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            result[i][j] = (rho * (1.0 / dt *
                ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dx) + (v[i + 1][j] - v[i - 1][j]) / (2.0 * dy)) -
                ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dx)) * ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dx)) -
                2.0 * ((u[i + 1][j] - u[i - 1][j]) / (2.0 * dy) * (v[i][j + 1] - v[i][j - 1]) / (2.0 * dx)) -
                ((v[i + 1][j] - v[i - 1][j]) / (2.0 * dy)) * ((v[i + 1][j] - v[i - 1][j]) / (2.0 * dy))));
        }
        // Periodic BC Pressure @ x = 2
        result[i][nx - 1] = rho * (1.0 / dt * ((u[i][0] - u[i][nx - 2]) / (2.0 * dx) +
                                           (v[i + 1][nx - 1] - v[i - 1][nx - 1]) / (2.0 * dy)) -
                               ((u[i][0] - u[i][nx - 2]) / (2.0 * dx)) * ((u[i][0] - u[i][nx - 2]) / (2.0 * dx)) -
                               2.0 * ((u[i + 1][nx - 1] - u[i - 1][nx - 1]) / (2.0 * dy) * (v[i][0] - v[i][nx - 2]) / (2.0 * dx)) -
                               ((v[i + 1][nx - 1] - v[i - 1][nx - 1]) / (2.0 * dy)) * ((v[i + 1][nx - 1] - v[i - 1][nx - 1]) / (2.0 * dy)));

        // Periodic BC Pressure @ x = 0
        result[i][0] = rho * (1.0 / dt * ((u[i][1] - u[i][nx - 1]) / (2.0 * dx) +
                                      (v[i + 1][0] - v[i - 1][0]) / (2.0 * dy)) -
                          ((u[i][1] - u[i][nx - 1]) / (2.0 * dx)) * ((u[i][1] - u[i][nx - 1]) / (2.0 * dx)) -
                          2.0 * ((u[i + 1][0] - u[i - 1][0]) / (2.0 * dy) * (v[i][1] - v[i][nx - 1]) / (2.0 * dx)) -
                          ((v[i + 1][0] - v[i - 1][0]) / (2.0 * dy)) * ((v[i + 1][0] - v[i - 1][0]) / (2.0 * dy)));
    }

    return result;
}

/*  function -- pressure poisson
    solves the poisson eqn to find out the pressure
    input - p, dx, dy, b, nit
    output - pressure vector
*/

std::vector<std::vector<double>> pressure_poisson(std::vector<std::vector<double>>& p,
                                                  double dx, double dy,
                                                  const std::vector<std::vector<double>>& b,
                                                  int nit) {
    std::vector<std::vector<double>> pn = p;

    for (int q = 0; q < nit; ++q) {
        pn = p;

        for (int i = 1; i < ny - 1; ++i) {
            for (int j = 1; j < nx - 1; ++j) {
                p[i][j] = ((pn[i][j + 1] + pn[i][j - 1]) * dy * dy +
                           (pn[i + 1][j] + pn[i - 1][j]) * dx * dx) /
                          (2.0 * (dx * dx + dy * dy)) -
                          (dx * dx * dy * dy) / (2.0 * (dx * dx + dy * dy)) * b[i][j];
            }
            // Periodic BC Pressure @ x = 2
            p[i][nx - 1] = (((pn[i][0] + pn[i][nx - 2]) * dy * dy +
                              (pn[i + 1][nx - 1] + pn[i - 1][nx - 1]) * dx * dx) /
                             (2.0 * (dx * dx + dy * dy)) -
                             dx * dx * dy * dy / (2.0 * (dx * dx + dy * dy)) * b[i][nx - 1]);

            // Periodic BC Pressure @ x = 0
            p[i][0] = (((pn[i][1] + pn[i][nx - 1]) * dy * dy +
                         (pn[i + 1][0] + pn[i - 1][0]) * dx * dx) /
                        (2.0 * (dx * dx + dy * dy)) -
                        dx * dx * dy * dy / (2.0 * (dx * dx + dy * dy)) * b[i][0]);
        }

        // Wall boundary conditions, pressure
        for (int j = 0; j < nx; ++j) {
            p[ny - 1][j] = p[ny - 2][j];  // dp/dy = 0 at y = 2
            p[0][j] = p[1][j];            // dp/dy = 0 at y = 0
        }
    }

    return p;
}

/*  function -- update_uv_boundary_conditions
    solves the poisson eqn to find out the pressure
    input - un, vn, p, dt, dx, dy, rho, nu, F
    output - updated value of u and v
*/

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
update_uv_boundary_conditions(std::vector<std::vector<double>>& un,
                              std::vector<std::vector<double>>& vn,
                              const std::vector<std::vector<double>>& p,
                              double dt, double dx, double dy, double rho, double nu, double F) {
    int ny = un.size();
    int nx = un[0].size();

    // definition of u and v
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));

    // calculate u and v from the momentum equation
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            u[i][j] = (un[i][j] -
                        un[i][j] * dt / dx *
                        (un[i][j] - un[i][j - 1]) -
                        vn[i][j] * dt / dy *
                        (un[i][j] - un[i - 1][j]) -
                        dt / (2.0 * rho * dx) *
                        (p[i][j + 1] - p[i][j - 1]) +
                        nu * (dt / (dx * dx) *
                        (un[i][j + 1] - 2.0 * un[i][j] + un[i][j - 1]) +
                        dt / (dy * dy) *
                        (un[i + 1][j] - 2.0 * un[i][j] + un[i - 1][j])) +
                        F * dt);

            v[i][j] = (vn[i][j] -
                        un[i][j] * dt / dx *
                        (vn[i][j] - vn[i][j - 1]) -
                        vn[i][j] * dt / dy *
                        (vn[i][j] - vn[i - 1][j]) -
                        dt / (2.0 * rho * dy) *
                        (p[i + 1][j] - p[i - 1][j]) +
                        nu * (dt / (dx * dx) *
                        (vn[i][j + 1] - 2.0 * vn[i][j] + vn[i][j - 1]) +
                        dt / (dy * dy) *
                        (vn[i + 1][j] - 2.0 * vn[i][j] + vn[i - 1][j])));
        }
    }

    // Periodic BC u @ x = 2
    for (int i = 1; i < ny - 1; ++i) {
        u[i][nx - 1] = (un[i][nx - 1] - un[i][nx - 1] * dt / dx *
                        (un[i][nx - 1] - un[i][nx - 2]) -
                        vn[i][nx - 1] * dt / dy *
                        (un[i][nx - 1] - un[i - 1][nx - 1]) -
                        dt / (2.0 * rho * dx) *
                        (p[i][0] - p[i][nx - 2]) +
                        nu * (dt / (dx * dx) *
                        (un[i][0] - 2.0 * un[i][nx - 1] + un[i][nx - 2]) +
                        dt / (dy * dy) *
                        (un[i + 1][nx - 1] - 2.0 * un[i][nx - 1] + un[i - 1][nx - 1])) + F * dt);
    }

    // Periodic BC u @ x = 0
    for (int i = 1; i < ny - 1; ++i) {
        u[i][0] = (un[i][0] - un[i][0] * dt / dx *
                   (un[i][0] - un[i][nx - 1]) -
                   vn[i][0] * dt / dy *
                   (un[i][0] - un[i - 1][0]) -
                   dt / (2.0 * rho * dx) *
                   (p[i][1] - p[i][nx - 1]) +
                   nu * (dt / (dx * dx) *
                   (un[i][1] - 2.0 * un[i][0] + un[i][nx - 1]) +
                   dt / (dy * dy) *
                   (un[i + 1][0] - 2.0 * un[i][0] + un[i - 1][0])) + F * dt);
    }

    // Periodic BC v @ x = 2
    for (int i = 1; i < ny - 1; ++i) {
        v[i][nx - 1] = (vn[i][nx - 1] - un[i][nx - 1] * dt / dx *
                        (vn[i][nx - 1] - vn[i][nx - 2]) -
                        vn[i][nx - 1] * dt / dy *
                        (vn[i][nx - 1] - vn[i - 1][nx - 1]) -
                        dt / (2.0 * rho * dy) *
                        (p[i + 1][nx - 1] - p[i - 1][nx - 1]) +
                        nu * (dt / (dx * dx) *
                        (vn[i][0] - 2.0 * vn[i][nx - 1] + vn[i][nx - 2]) +
                        dt / (dy * dy) *
                        (vn[i + 1][nx - 1] - 2.0 * vn[i][nx - 1] + vn[i - 1][nx - 1])));
    }

    // Periodic BC v @ x = 0
    for (int i = 1; i < ny - 1; ++i) {
        v[i][0] = (vn[i][0] - un[i][0] * dt / dx *
                   (vn[i][0] - vn[i][nx - 1]) -
                   vn[i][0] * dt / dy *
                   (vn[i][0] - vn[i - 1][0]) -
                   dt / (2.0 * rho * dy) *
                   (p[i + 1][0] - p[i - 1][0]) +
                   nu * (dt / (dx * dx) *
                   (vn[i][1] - 2.0 * vn[i][0] + vn[i][nx - 1]) +
                   dt / (dy * dy) *
                   (vn[i + 1][0] - 2.0 * vn[i][0] + vn[i - 1][0])));
    }

    // Wall BC: u,v = 0 @ y = 0,2
    for (int j = 0; j < nx; ++j) {
        u[0][j] = 0;
        u[ny - 1][j] = 0;
        v[0][j] = 0;
        v[ny - 1][j] = 0;
    }

    return std::make_pair(u, v);
}

/*  function -- calculate_udiff
    udiff is the error at every step
    input - u, un
    output - udiff
*/

double calculate_udiff(const std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& un) {
    double sum_u = 0.0;
    double sum_un = 0.0;
    double udiff = 0.0;

    for (size_t i = 0; i < u.size(); ++i) {
        for (size_t j = 0; j < u[0].size(); ++j) {
            sum_u += u[i][j];
            sum_un += un[i][j];
        }
    }
    udiff = (sum_u - sum_un) / sum_u;

    return udiff;
}

/*  function -- linspace
    just a replication of linspace command in python
    output - vector of numPoints
*/

std::vector<double> linspace(double start, double end, int numPoints) {
    std::vector<double> result(numPoints);
    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
        result[i] = start + i * step;
    }
    return result;
}


int main() {

    // Variable Initialization
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> un(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> vn(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> pn(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx, 0.0));

    // x & y variable declaration
    std::vector<double> x = linspace(0.0, 2.0, nx);
    std::vector<double> y = linspace(0.0, 2.0, ny);

    double udiff = 1.0;
    int stepcount = 1;

    // Main loop to solve the equation
    while (udiff > 0.001) {

        un = u;  // reassigning the un to already calculated u
        vn = v;  // reassigning the vn to already calculated v

        // using the build_up_b function to calculate b
        b = build_up_b(b,rho, dt, u, v, dx, dy);

        // using the pressure_poisson equation to calculate p
        p = pressure_poisson(p, dx, dy, b, nit);


        // Update the boundary conditions using the update_uv_boundary_conditions function
        std::tie(u, v) = update_uv_boundary_conditions(un, vn, p, dt, dx, dy, rho, nu, F);

        // calculating udiff using calculate_udiff function
        udiff = calculate_udiff(u,un);
        printf("%f \n", udiff);
        stepcount++;  // calculating stepcount to find out number of iterations
    }

    // writing the x, y, u, v file to text files
    // These files can then later be used for post-processing in python
    std::ofstream file;

    // writing x text file
    file.open("x.txt");
    for (const double &val : x) {
        file << val << "\n";
    }
    file.close();

    // writing y text file
    file.open("y.txt");
    for (const double &val : y) {
        file << val << "\n";
    }
    file.close();

    // writing u and v text files
    std::ofstream u_file("u.txt");
    std::ofstream v_file("v.txt");


    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            u_file << u[i][j] << " ";
            v_file << v[i][j] << " ";
        }
        u_file << "\n";
        v_file << "\n";
    }

    u_file.close();
    v_file.close();

    return 0;
}
