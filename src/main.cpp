#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// hdf5 library
#include <H5Cpp.h>

#include "../include/NS.hpp"
#include "../include/misc.hpp"
#include "../include/object.hpp"
#include "../include/write.hpp"

using namespace Eigen;

#define EPS 1e-8
// #define interpolate_boundary 1

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

#define WRITE_ANIM() ((t - EPS < (plot_count + 1) * plot_dt) && (t + prm.dt - EPS > (plot_count + 1) * plot_dt))

int main(void) {
  const string filename_input = "config/input.txt";  // name of the input file to read the parameters
  const string space = "    ";                       // space to print
  const uint per = 10;                               // progress percentage interval to print (each count%)
  uint count = per;
  Prm prm;
  double T;             // final time of the simulation
  bool animation;       // flag to enable or disable the animation
  bool vorticity;       // flag to enable or disable plotting the vorticity
  double plot_dt;       // frequency to write the animation file
  uint plot_count = 0;  // counter to write the animation file
  string object_type;   // type of object to simulate
  int64_t total_files = 0, total_advection = 0, total_diffusion = 0,
          total_others = 0, total_build_poisson = 0, total_solve_laplace = 0,
          total_constr_obstacle = 0, total_boundary = 0;  // variables to measure the time
  chrono::steady_clock::time_point begin, end;            // variables to measure the time

  // ------------- File input setup ----------------
  START_TIMER();
  ifstream file_input;
  file_input.open(filename_input);
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> prm.LX;
    file_input >> tmp >> prm.LY;
    file_input >> tmp >> prm.nx;
    file_input >> tmp >> prm.ny;
    file_input >> tmp >> prm.dt;
    file_input >> tmp >> T;
    file_input >> tmp >> prm.Re;
    file_input >> tmp >> prm.L;
    file_input >> tmp >> prm.nu;
    file_input >> tmp >> vorticity;
    file_input >> tmp >> prm.obstacle_ON;
    file_input >> tmp >> object_type;
    file_input >> tmp >> plot_dt;
    file_input >> tmp >> animation;
  }
  file_input.close();

  prm.NX = prm.nx + 2;
  prm.NY = prm.ny + 2;
  prm.dx = prm.LX / prm.nx;
  prm.dy = prm.LY / prm.ny;
  prm.NXNY = (uint)(prm.NX * prm.NY);
  prm.nxny = (uint)(prm.nx * prm.ny);
  prm.U = prm.Re * prm.nu / prm.L;

  // ------------- Print plot setup -----------------
  cout << "dx:            " << space << prm.dx << endl;
  cout << "dy:            " << space << prm.dy << endl;
  cout << "dt:            " << space << prm.dt << endl;
  cout << "T:             " << space << T << endl;
  cout << "Re:            " << space << prm.Re << endl;
  cout << "U:             " << space << prm.U << endl;
  cout << "Obstacle?      " << space << (prm.obstacle_ON ? "yes, " + object_type : "no") << endl;
  cout << "Vorticity?     " << space << (vorticity ? "yes" : "no") << endl;
  cout << "Plot animation?" << space << (animation ? "yes" : "no") << endl;
  // ------------------------------------------------

  double t = 0.0;
  uint numSteps = 0;
  double fraction_completed = T / 100.;  // fraction of the integration time to print

  // allocate memory
  double* u = new double[prm.NXNY];      // x component of velocity
  double* ustar = new double[prm.NXNY];  // x component of velocity
  double* v = new double[prm.NXNY];      // y component of velocity
  double* vstar = new double[prm.NXNY];  // y component of velocity
  double* w = new double[prm.NXNY];      // vorticity (only 1 component, because it is 2D)
  double* p = new double[prm.NXNY];      // pressure
  double* adv_u = new double[prm.NXNY];  // advection term of u
  double* adv_v = new double[prm.NXNY];  // advection term of v
  // initialize to 0
  for (int i = 0; i < prm.NXNY; i++) {
    u[i] = 0;
    ustar[i] = 0;
    v[i] = 0;
    vstar[i] = 0;
    w[i] = 0;
    p[i] = 0;
    adv_u[i] = 0;
    adv_v[i] = 0;
  }
  // double* u = (double*)calloc(prm.NXNY, sizeof(double));      // x component of velocity, initialized to 0 by calloc
  // double* ustar = (double*)calloc(prm.NXNY, sizeof(double));  // x component of velocity, initialized to 0 by calloc
  // double* v = (double*)calloc(prm.NXNY, sizeof(double));      // y component of velocity, initialized to 0 by calloc
  // double* vstar = (double*)calloc(prm.NXNY, sizeof(double));  // y component of velocity, initialized to 0 by calloc
  // double* w = (double*)calloc(prm.NXNY, sizeof(double));      // vorticity (only 1 component, because it is 2D), initialized to 0 by calloc
  // double* p = (double*)calloc(prm.NXNY, sizeof(double));      // pressure, initialized to 0 by calloc
  // // DO NOT ACCESS THE GHOST POINTS OF THESE ARRAYS
  // double* adv_u = (double*)malloc(prm.NXNY * sizeof(double));  // advection term of u
  // double* adv_v = (double*)malloc(prm.NXNY * sizeof(double));  // advection term of v

  // if (animation) {
  // file_output_w.open(filename_solution_w);
  // file_output_u.open(filename_solution_u);
  // write_sol_w(file_output_w, w, 0, prm);
  // write_sol(file_output_u, u, v, 0, prm);
  // }
  saveSetupToHDF5(prm, T, vorticity, animation);
  saveDataToHDF5(plot_count, u, v, w, p, prm.NX, prm.NY, t);
  plot_count++;
  END_TIMER();
  ADD_TIME_TO(total_files);
  // -------------------------------------------

  // ------------- Obstacle setup ----------------
  START_TIMER();
  Object* obstacle;
  if (object_type == "circle") {
    obstacle = new Circle(prm.LX / 8, prm.LY / 2, prm.L / 2, prm);  // L is the diameter of the circle
  } else if (object_type == "rectangle") {
    obstacle = new Rectangle(prm.LX / 8, prm.LY / 2, prm.L, prm.L, prm);  // L is the side of the square
  } else if (object_type == "mountain") {
    obstacle = new Mountain(0.7, 0.9, 0.2, 3.1, prm);
  } else if (object_type == "airfoil") {
    obstacle = new Airfoil(0.17814, -0.0756, -0.21096, 0.17058, -0.0609, 0.2, prm.LX / 8, prm.LY / 2, prm);
  } else {
    cout << "Object type not recognized" << endl;
    return 1;
  }
  END_TIMER();
  ADD_TIME_TO(total_constr_obstacle);
  // -------------------------------------------

  // ------------- Build Poission matrix -----------
  START_TIMER();
  vector<Trip> coeffs;
  buildPoissonMatrix(coeffs, prm);

  SpMat A(prm.nx * prm.ny, prm.nx * prm.ny);
  VectorXd div(prm.nx * prm.ny);
  A.setFromTriplets(coeffs.begin(), coeffs.end());

  SimplicialLDLT<SpMat> chol;
  chol.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (chol.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  VectorXd p_solved;
  END_TIMER();
  ADD_TIME_TO(total_build_poisson);
  // -------------------------------------------

  numSteps++;
  double max_u_v = 0, lap, rhs_u, rhs_v;

  while (t < T - EPS) {
    if (prm.dt < EPS) {
      cout << "Time step too small. Exiting." << endl;
      return 1;
    }
    // adaptative time step
    for (int i = 0; i < prm.NXNY; i++) {
      max_u_v = max(max_u_v, max(abs(u[i]), abs(v[i])));
    }
    if (max_u_v > 1e-3) {
      // update CFL number
      prm.dt = prm.dx * prm.dy / (2 * max_u_v * (prm.dx + prm.dy));  // the factor 2 is to be conservative
    }

    // ---------- advection term ------------
    START_TIMER();
    Semilag2(u, v, u, adv_u, prm, *obstacle);
    Semilag2(u, v, v, adv_v, prm, *obstacle);
    END_TIMER();
    ADD_TIME_TO(total_advection);
    // --------------------------------------

    // ---------- diffusion term ------------
    START_TIMER();
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        if (prm.obstacle_ON && obstacle->IsInside[i * prm.NY + j]) {
          Ustar(i, j) = 0;
          Vstar(i, j) = 0;
          continue;
        }
        lap = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) / (prm.dx * prm.dx) +
              (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) / (prm.dy * prm.dy);
        rhs_u = ADV_U(i, j) + prm.dt * lap / prm.Re;
        lap = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) / (prm.dx * prm.dx) +
              (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) / (prm.dy * prm.dy);
        rhs_v = ADV_V(i, j) + prm.dt * lap / prm.Re;
        Ustar(i, j) = rhs_u;
        Vstar(i, j) = rhs_v;
      }
    }
    END_TIMER();
    ADD_TIME_TO(total_diffusion);
    // --------------------------------------

    // update velocity ghost points
    START_TIMER();
    BC_velocity(ustar, vstar, prm);
    END_TIMER();
    ADD_TIME_TO(total_boundary);

    // ---------- pressure term ------------
    START_TIMER();
    // compute the (minus) divergence of the velocity field (we omit multiplication by dt, because later on we will divide by dt) and store it in div
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        if (prm.obstacle_ON) {
          DIV(i - 1, j - 1) = obstacle->IsInside[i * prm.NY + j] ? 0 : -(Ustar(i + 1, j) - Ustar(i - 1, j)) / (2 * prm.dx) - (Vstar(i, j + 1) - Vstar(i, j - 1)) / (2 * prm.dy);
        } else {
          DIV(i - 1, j - 1) = -(Ustar(i + 1, j) - Ustar(i - 1, j)) / (2 * prm.dx) -
                              (Vstar(i, j + 1) - Vstar(i, j - 1)) / (2 * prm.dy);
        }
      }
    }
    END_TIMER();
    ADD_TIME_TO(total_others);

    START_TIMER();
    // solve the pressure Poisson equation

    p_solved = chol.solve(div);
    if (chol.info() != Success) {
      cout << "Cholesky solve failed" << endl;
      return 1;
    }
    // convert p_solved back to pointer
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        P(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
      }
    }

    END_TIMER();
    ADD_TIME_TO(total_solve_laplace);

    // update pressure ghost points
    START_TIMER();
    BC_pressure(p, prm);
    END_TIMER();
    ADD_TIME_TO(total_boundary);

    START_TIMER();
    // ---------- projection step ------------
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        U(i, j) = Ustar(i, j) - (P(i + 1, j) - P(i - 1, j)) / (2 * prm.dx);
        V(i, j) = Vstar(i, j) - (P(i, j + 1) - P(i, j - 1)) / (2 * prm.dy);
      }
    }
    BC_velocity(u, v, prm);
    END_TIMER();
    ADD_TIME_TO(total_others);
    // --------------------------------------

    t += prm.dt;

    // ---------- printing to file ------------
    START_TIMER();
    if (WRITE_ANIM()) {
      // we cannot change the value of u and v (in order not to produce divergence with the pressure at next step)
      // linear interpolation on the outer discrete boundary of the object to smooth the boundary jump values.
      memcpy(ustar, u, prm.NXNY * sizeof(double));
      memcpy(vstar, v, prm.NXNY * sizeof(double));
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          if (prm.obstacle_ON && obstacle->IsInside[i * prm.NY + j]) {
            Ustar(i, j) = 0;
            Vstar(i, j) = 0;
          }
#ifdef interpolate_boundary
          if (prm.obstacle_ON && obstacle->IsInterface[i * prm.NY + j]) {
            Ustar(i, j) = interpolate(i, j, ustar, prm, *obstacle);
            Vstar(i, j) = interpolate(i, j, vstar, prm, *obstacle);
          }
#endif
        }
      }
      set_vorticity(ustar, vstar, w, prm);
      saveDataToHDF5(plot_count, ustar, vstar, w, p, prm.NX, prm.NY, t);
      // write_sol_w(file_output_w, w, t, prm);
      // write_sol(file_output_u, ustar, vstar, t, prm);
      plot_count++;
    }

    if (t > fraction_completed * count - EPS) {
      cout << count << "%"
           << " dt: " << prm.dt << endl;
      count += per;
    }
    END_TIMER();
    ADD_TIME_TO(total_files);
    // --------------------------------------

    numSteps++;
  }

  // ---------- auxiliar file to pass to the animation script ------------
  START_TIMER();
  if (animation) {
    // if (prm.obstacle_ON)
    write_tmp_file(object_type, obstacle->data);
    // else
    //   write_tmp_file(prm.LX, prm.LY, vorticity);
  }
  END_TIMER();
  ADD_TIME_TO(total_files);
  // --------------------------------------------------------------------

  print("Total time for files:                   " + space, total_files);
  print("Total time for construction of obstacle:" + space, total_constr_obstacle);
  print("Total time for build laplace:           " + space, total_build_poisson);
  print("Total time for advection:               " + space, total_advection);
  print("Total time for diffusion:               " + space, total_diffusion);
  print("Total time for solve laplace:           " + space, total_solve_laplace);
  print("Total time for boundary:                " + space, total_boundary);
  print("Total time for rest:                    " + space, total_others);

  // free memory
  delete[] u;
  delete[] ustar;
  delete[] v;
  delete[] vstar;
  delete[] w;
  delete[] p;
  delete[] adv_u;
  delete[] adv_v;
  delete obstacle;

  return 0;
}
