
/****
Author - Swadesh Suman
Code - Navier Stokes 2D Channel Flow
Time - forward stepping
Spatial - backward discretization - first order, Centered discretization - second order
Location - McGill University, Canada
****/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


// Function for calculating the second term in possion's equation
std::vector<std::vector<double>> build_up_b(double rho, double dt, double dx, double dy,
                                           const std::vector<std::vector<double>>& u,
                                           const std::vector<std::vector<double>>& v) {
    int nx = u.size();
    int ny = u[0].size();

    std::vector<std::vector<double>> b(nx, std::vector<double>(ny, 0.0));

    // calculation of b term
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            b[i][j] = (rho * (1 / dt * ((u[i][j + 1] - u[i][j - 1]) / (2 * dx) +
                                        (v[i + 1][j] - v[i - 1][j]) / (2 * dy)) -
                              ((u[i][j + 1] - u[i][j - 1]) / (2 * dx)) * ((u[i][j + 1] - u[i][j - 1]) / (2 * dx)) -
                              2 * ((u[i + 1][j] - u[i - 1][j]) / (2 * dy) * (v[i][j + 1] - v[i][j - 1]) / (2 * dx)) -
                              ((v[i + 1][j] - v[i - 1][j]) / (2 * dy)) * ((v[i + 1][j] - v[i - 1][j]) / (2 * dy))));
        }
    }

    // Periodic BC Pressure @ x = 2
    for (int i = 1; i < nx - 1; ++i) {
        b[i][ny - 1] = (rho * (1 / dt * ((u[i][0] - u[i][ny - 2]) / (2 * dx) +
                                         (v[i + 1][ny - 1] - v[i - 1][ny - 1]) / (2 * dy)) -
                               ((u[i][0] - u[i][ny - 2]) / (2 * dx)) * ((u[i][0] - u[i][ny - 2]) / (2 * dx)) -
                               2 * ((u[i + 1][ny - 1] - u[i - 1][ny - 1]) / (2 * dy) * (v[i][0] - v[i][ny - 2]) / (2 * dx)) -
                               ((v[i + 1][ny - 1] - v[i - 1][ny - 1]) / (2 * dy)) * ((v[i + 1][ny - 1] - v[i - 1][ny - 1]) / (2 * dy))));
    }

    // Periodic BC Pressure @ x = 0
    for (int i = 1; i < nx - 1; ++i) {
        b[i][0] = (rho * (1 / dt * ((u[i][1] - u[i][ny - 1]) / (2 * dx) +
                                    (v[i + 1][0] - v[i - 1][0]) / (2 * dy)) -
                          ((u[i][1] - u[i][ny - 1]) / (2 * dx)) * ((u[i][1] - u[i][ny - 1]) / (2 * dx)) -
                          2 * ((u[i + 1][0] - u[i - 1][0]) / (2 * dy) * (v[i][1] - v[i][ny - 1]) / (2 * dx)) -
                          ((v[i + 1][0] - v[i - 1][0]) / (2 * dy)) * ((v[i + 1][0] - v[i - 1][0]) / (2 * dy))));
    }

    return b;
}

// Function to calculate the pressure in the poisson equation
std::vector<std::vector<double>> pressure_poisson_periodic(const std::vector<std::vector<double>>& p, const double dx, const double dy, const int nit, std::vector<std::vector<double>> b) {
    int nx = p.size();
    int ny = p[0].size();

    std::vector<std::vector<double>> pn(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> result(nx, std::vector<double>(ny, 0.0));

    for (int q = 0; q < nit; ++q) {
        pn = result;

        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                result[i][j] = (((pn[i][j + 1] + pn[i][j - 1]) * dy * dy +
                                 (pn[i + 1][j] + pn[i - 1][j]) * dx * dx) /
                                (2 * (dx * dx + dy * dy)) -
                                dx * dx * dy * dy / (2 * (dx * dx + dy * dy)) * b[i][j]);
            }

            // Periodic BC Pressure @ x = 2
            result[i][ny - 1] = (((pn[i][0] + pn[i][ny - 2]) * dy * dy +
                                  (pn[i + 1][ny - 1] + pn[i - 1][ny - 1]) * dx * dx) /
                                 (2 * (dx * dx + dy * dy)) -
                                 dx * dx * dy * dy / (2 * (dx * dx + dy * dy)) * b[i][ny - 1]);

            // Periodic BC Pressure @ x = 0
            result[i][0] = (((pn[i][1] + pn[i][ny - 1]) * dy * dy +
                             (pn[i + 1][0] + pn[i - 1][0]) * dx * dx) /
                            (2 * (dx * dx + dy * dy)) -
                            dx * dx * dy * dy / (2 * (dx * dx + dy * dy)) * b[i][0]);
        }

        // Wall boundary conditions, pressure
        for (int j = 0; j < ny; ++j) {
            result[nx - 1][j] = result[nx - 2][j];  // dp/dy = 0 at y = 2
            result[0][j] = result[1][j];  // dp/dy = 0 at y = 0
        }
    }

    return result;
}

void saveDataToFile(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream file(filename);
    for (const auto& row : data) {
        for (double value : row) {
            file << value << " ";
        }
        file << std::endl;
    }
    file.close();
}

//*********************************************************************///
int main() {
    // Variable declarations
    int nx = 41;  // points in x-direction
    int ny = 41;  // points in y-direction
    int nt = 10;    // no of time steps
    int nit = 50;   // no of pseudo time step for converging poisson equation
    double c = 1;  // constant value
    double dx = 2.0 / (nx - 1);  // grid spacing in x-direction
    double dy = 2.0 / (ny - 1);  // grid spacing in y -direction

    // Create a 2D vector initialized with zeros
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0)); // final u - velocity values
    std::vector<std::vector<double>> un(ny, std::vector<double>(nx, 0.0));  // initial u - velocity values
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));   // final v - velocity values
    std::vector<std::vector<double>> vn(ny, std::vector<double>(nx, 0.0)); // initial v - velocity values
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx, 0.0)); // final pressure values
    std::vector<std::vector<double>> pn(ny, std::vector<double>(nx, 0.0));  // initial pressure values
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx, 0.0));   // initial b term values



    // Grid points generation for the structured grid
    std::vector<double> x(nx);
    std::vector<double> y(ny);

    double start_x = 0.0;
    double end_x = 2.0;
    double start_y = 0.0;
    double end_y = 2.0;

    // x and y values in the grid
    for (int i = 0; i < nx; ++i) {
        x[i] = start_x + i * (end_x - start_x) / (nx - 1);
    }

    for (int i = 0; i < ny; ++i) {
        y[i] = start_y + i * (end_y - start_y) / (ny - 1);
    }

    //
    double rho = 1; //density
    double nu = 0.1;    // viscosity
    double F = 1;  // constant term in x-momentum equation to act as source term
    double dt = 0.01;  // time stepping


    double udiff = 1;  // initial error
    int stepcount = 0; // to count the number of iterations

    // Main simulation loop
    while (udiff > 0.05) {

        b = build_up_b(rho, dt, dx, dy, u, v);
        p = pressure_poisson_periodic(p, dx, dy, nit, b);

        // updating the u and v velocity values
        for (int i = 1; i < u.size() - 1; ++i) {
            for (int j = 1; j < u[i].size() - 1; ++j) {
                u[i][j] = (un[i][j] -
                           un[i][j] * dt / dx * (un[i][j] - un[i][j - 1]) -
                           vn[i][j] * dt / dy * (un[i][j] - un[i - 1][j]) -
                           dt / (2 * rho * dx) * (p[i][j + 1] - p[i][j - 1]) +
                           nu * (dt / dx * dx * (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) +
                                 dt / dy * dy * (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j])) + F * dt);

                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx * (vn[i][j] - vn[i][j - 1]) -
                           vn[i][j] * dt / dy * (vn[i][j] - vn[i - 1][j]) -
                           dt / (2 * rho * dy) * (p[i + 1][j] - p[i - 1][j]) +
                           nu * (dt / dx * dx * (vn[i][j + 1] - 2 * vn[i][j] + vn[i][j - 1]) +
                                 dt / dy * dy * (vn[i + 1][j] - 2 * vn[i][j] + vn[i - 1][j])));
            }
        }

        //  updating the u and v velocity values at the outlet
        int m = -1;
        // Periodic BC u @ x = 2
        u[u.size() - 1][m] = (un[u.size() - 1][m] - un[u.size() - 1][m] * dt / dx *
                              (un[u.size() - 1][m] - un[u.size() - 1][m - 1]) -
                               vn[u.size() - 1][m] * dt / dy *
                              (un[u.size() - 1][m] - un[u.size() - 2][m]) -
                               dt / (2 * rho * dx) *
                              (p[u.size() - 1][0] - p[u.size() - 1][m - 1]) +
                               nu * (dt / (dx * dx) *
                              (un[u.size() - 1][0] - 2 * un[u.size() - 1][m] + un[u.size() - 1][m - 1]) +
                               dt / (dy * dy) *
                              (un[1][m] - 2 * un[u.size() - 1][m] + un[u.size() - 2][m])) + F * dt);
        // Periodic BC v @ x = 2
        v[v.size() - 1][m] = (vn[v.size() - 1][m] - un[v.size() - 1][m] * dt / dx *
                              (vn[v.size() - 1][m] - vn[v.size() - 1][m - 1]) -
                               vn[v.size() - 1][m] * dt / dy *
                              (vn[v.size() - 1][m] - vn[v.size() - 2][m]) -
                               dt / (2 * rho * dy) *
                              (p[1][m] - p[u.size() - 1][m]) +
                               nu * (dt / (dx * dx) *
                              (vn[v.size() - 1][0] - 2 * vn[v.size() - 1][m] + vn[v.size() - 1][m - 1]) +
                               dt / (dy * dy) *
                              (vn[1][m] - 2 * vn[v.size() - 1][m] + vn[v.size() - 2][m])));

        //  updating the u and v velocity values at the inlet
        int k = 0;
        // Periodic BC u @ x = 0
        u[0][k] = (un[0][k] - un[0][k] * dt / dx *
                   (un[0][k] - un[0][k - 1]) -
                    vn[0][k] * dt / dy *
                   (un[0][k] - un[u.size() - 1][k]) -
                    dt / (2 * rho * dx) *
                   (p[0][1] - p[0][k - 1]) +
                    nu * (dt / (dx * dx) *
                   (un[0][1] - 2 * un[0][k] + un[0][k - 1]) +
                    dt / (dy * dy) *
                   (un[1][k] - 2 * un[0][k] + un[u.size() - 1][k])) + F * dt);


        // Periodic BC v @ x = 0
        v[0][k] = (vn[0][k] - un[0][k] * dt / dx *
                   (vn[0][k] - vn[0][k - 1]) -
                    vn[0][k] * dt / dy *
                   (vn[0][k] - vn[u.size() - 1][k]) -
                    dt / (2 * rho * dy) *
                   (p[1][k] - p[u.size() - 1][k]) +
                    nu * (dt / (dx * dx) *
                   (vn[0][1] - 2 * vn[0][k] + vn[0][k - 1]) +
                    dt / (dy * dy) *
                   (vn[1][k] - 2 * vn[0][k] + vn[u.size() - 1][k])));

        //  updating the u and v velocity values at the upper and lower wall
        for (int i = 0; i < u.size(); ++i) {
            u[0][i] = 0;
            u[u.size() - 1][i] = 0;
            v[0][i] = 0;
            v[v.size() - 1][i] = 0;
        }

        double sumU = 0.0;
        double sumUN = 0.0;

        // Calculating u difference as an error
        for (int i = 0; i < u.size(); ++i) {
            for (int j = 0; j < u[i].size(); ++j) {
                sumU += u[i][j];
                sumUN += un[i][j];
            }
        }
        udiff = (sumU - sumUN) / sumU;
        /*
        if (std::abs(sumU) > 1e-10) {  // Avoid division by zero

        } else {
            udiff = 0.0; // Handle the case when sumU is close to zero
        }
        */
        stepcount++;
        printf("%f \n", udiff);

    }

    // Save initial values to a file
    saveDataToFile(u, "values_u.txt");
    saveDataToFile(v, "values_v.txt");


    return 0;
}
    /*
    double** u = new double*[ny];
    double** un = new double*[ny];
    double** v = new double*[ny];
    double** vn = new double*[ny];
    double** p = new double*[ny];
    double** pn = new double*[ny];
    double** b = new double*[ny];

    for (int i = 0; i < ny; ++i) {
        u[i] = new double[nx];
        un[i] = new double[nx];
        v[i] = new double[nx];
        vn[i] = new double[nx];
        p[i] = new double[nx];
        pn[i] = new double[nx];
        b[i] = new double[nx];
    }

    */
    //u , v = solve_navier_stokes(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, double rho, double dt, double dx, double dy, double nu, double F, int nit);
    /*
    // Save initial and final values to a text file
    std::ofstream file("result.txt");
    // Save initial values
    file << "Initial Values:\n";
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            file << u[i][j] << " " << v[i][j] << " " << p[i][j] << "\n";
        }
    }
    // Save final values
    file << "\nFinal Values:\n";
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            file << un[i][j] << " " << vn[i][j] << " " << pn[i][j] << "\n";
        }
    }

    // Close the file
    file.close();
    */
/*
    // Initialize the pressure matrix
    std::vector<std::vector<double>> pressure(nx, std::vector<double>(ny, 0.0));

    // Call the function and get the result
    std::vector<std::vector<double>> result = pressure_poisson_periodic(pressure, dx, dy, nit);

    // Initialize the pressure matrix
    std::vector<std::vector<double>> pressure(nx, std::vector<double>(ny, 0.0));

    // Call the function
    pressure_poisson_periodic(pressure, dx, dy, nit);
std::vector<std::vector<double>> result = build_up_b(rho, dt, dx, dy, u, v);


std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> solve_navier_stokes(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, double rho, double dt, double dx, double dy, double nu, double F, int nit) {
    int stepcount = 0;
    double udiff = 1.0;

    while (udiff > 0.001) {
        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;

        std::vector<std::vector<double>> b = build_up_b(rho, dt, dx, dy, u, v);
        std::vector<std::vector<double>> p = pressure_poisson_periodic(p, dx, dy, nit, b);

        for (int i = 1; i < u.size() - 1; ++i) {
            for (int j = 1; j < u[i].size() - 1; ++j) {
                u[i][j] = (un[i][j] -
                           un[i][j] * dt / dx * (un[i][j] - un[i][j - 1]) -
                           vn[i][j] * dt / dy * (un[i][j] - un[i - 1][j]) -
                           dt / (2 * rho * dx) * (p[i][j + 1] - p[i][j - 1]) +
                           nu * (dt / dx * dx * (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) +
                                 dt / dy * dy * (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j])) + F * dt);

                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx * (vn[i][j] - vn[i][j - 1]) -
                           vn[i][j] * dt / dy * (vn[i][j] - vn[i - 1][j]) -
                           dt / (2 * rho * dy) * (p[i + 1][j] - p[i - 1][j]) +
                           nu * (dt / dx * dx * (vn[i][j + 1] - 2 * vn[i][j] + vn[i][j - 1]) +
                                 dt / dy * dy * (vn[i + 1][j] - 2 * vn[i][j] + vn[i - 1][j])));
            }
        }

        int j = 2;
        // Periodic BC u @ x = 2
        u[u.size() - 1][j] = (un[u.size() - 1][j] - un[u.size() - 1][j] * dt / dx *
                              (un[u.size() - 1][j] - un[u.size() - 1][j - 1]) -
                               vn[u.size() - 1][j] * dt / dy *
                              (un[u.size() - 1][j] - un[u.size() - 2][j]) -
                               dt / (2 * rho * dx) *
                              (p[u.size() - 1][0] - p[u.size() - 1][j - 1]) +
                               nu * (dt / (dx * dx) *
                              (un[u.size() - 1][0] - 2 * un[u.size() - 1][j] + un[u.size() - 1][j - 1]) +
                               dt / (dy * dy) *
                              (un[1][j] - 2 * un[u.size() - 1][j] + un[u.size() - 2][j])) + F * dt);
        // Periodic BC v @ x = 2
        v[v.size() - 1][j] = (vn[v.size() - 1][j] - un[v.size() - 1][j] * dt / dx *
                              (vn[v.size() - 1][j] - vn[v.size() - 1][j - 1]) -
                               vn[v.size() - 1][j] * dt / dy *
                              (vn[v.size() - 1][j] - vn[v.size() - 2][j]) -
                               dt / (2 * rho * dy) *
                              (p[1][j] - p[u.size() - 1][j]) +
                               nu * (dt / (dx * dx) *
                              (vn[v.size() - 1][0] - 2 * vn[v.size() - 1][j] + vn[v.size() - 1][j - 1]) +
                               dt / (dy * dy) *
                              (vn[1][j] - 2 * vn[v.size() - 1][j] + vn[v.size() - 2][j])));

        int k = 0;
        // Periodic BC u @ x = 0
        u[0][k] = (un[0][k] - un[0][k] * dt / dx *
                   (un[0][k] - un[0][k - 1]) -
                    vn[0][k] * dt / dy *
                   (un[0][k] - un[u.size() - 1][k]) -
                    dt / (2 * rho * dx) *
                   (p[0][1] - p[0][k - 1]) +
                    nu * (dt / (dx * dx) *
                   (un[0][1] - 2 * un[0][k] + un[0][j - 1]) +
                    dt / (dy * dy) *
                   (un[1][k] - 2 * un[0][k] + un[u.size() - 1][k])) + F * dt);


        // Periodic BC v @ x = 0
        v[0][k] = (vn[0][k] - un[0][k] * dt / dx *
                   (vn[0][k] - vn[0][k - 1]) -
                    vn[0][k] * dt / dy *
                   (vn[0][k] - vn[u.size() - 1][k]) -
                    dt / (2 * rho * dy) *
                   (p[1][k] - p[u.size() - 1][k]) +
                    nu * (dt / (dx * dx) *
                   (vn[0][1] - 2 * vn[0][k] + vn[0][k - 1]) +
                    dt / (dy * dy) *
                   (vn[1][k] - 2 * vn[0][k] + vn[u.size() - 1][k])));

        // Wall BC: u, v = 0 @ y = 0, 2
        for (int i = 0; i < u.size(); ++i) {
            u[0][i] = 0;
            u[u.size() - 1][i] = 0;
            v[0][i] = 0;
            v[v.size() - 1][i] = 0;
        }

        double sumU = 0.0;
        double sumUN = 0.0;

        // Calculate udiff...
        for (int i = 0; i < u.size(); ++i) {
            for (int j = 0; j < u[i].size(); ++j) {
                sumU += u[i][j];
                sumUN += un[i][j];
            }
        }

        if (std::abs(sumU) > 1e-10) {  // Avoid division by zero
            udiff = (sumU - sumUN) / sumU;
        } else {
            udiff = 0.0; // Handle the case when sumU is close to zero
        }

                stepcount++;
    }
    return std::make_pair(u, v);
    //return u,v;
}
*/
