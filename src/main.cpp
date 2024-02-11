#include <chrono>
#include <cstring>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../include/NS.hpp"
#include "../include/misc.hpp"
#include "../include/write.hpp"

using namespace Eigen;

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

#define WRITE_ANIM() ((animation) && (numSteps % freq_write_anim == 0))

int main(void) {
  const string filename_solution = "data/solution.txt";  // name of the output file to write the solution
  const string filename_input = "data/input.txt";        // name of the input file to read the parameters
  const string space = "    ";                           // space to print
  Prm prm;
  double T;                                                                                     // final time of the simulation
  bool animation;                                                                               // flag to enable or disable the animation
  uint freq_write_anim;                                                                         // frequency to write the animation file
  int64_t total_files = 0, total_others = 0, total_build_laplace = 0, total_solve_laplace = 0;  // variables to measure the time
  chrono::steady_clock::time_point begin, end;                                                  // variables to measure the time

  // ------------- File input setup ----------------
  START_TIMER();
  ifstream file_input;
  ofstream file_output_anim;
  file_input.open("data/input.txt");
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> prm.LX;
    file_input >> tmp >> prm.LY;
    file_input >> tmp >> prm.nx;
    file_input >> tmp >> prm.ny;
    file_input >> tmp >> prm.dt;
    file_input >> tmp >> T;
    file_input >> tmp >> prm.Re;
    file_input >> tmp >> animation;
    file_input >> tmp >> freq_write_anim;
  }
  file_input.close();

  prm.NX = prm.nx + 2;
  prm.NY = prm.ny + 2;
  prm.dx = prm.LX / prm.nx;
  prm.dy = prm.LY / prm.ny;
  prm.size_ghost = (uint)(prm.NX * prm.NY);
  prm.size = (uint)(prm.nx * prm.ny);

  // ------------- Plot setup ----------------
  cout << "dx:            " << space << prm.dx << endl;
  cout << "dy:            " << space << prm.dy << endl;
  cout << "dt:            " << space << prm.dt << endl;
  cout << "T:             " << space << T << endl;
  cout << "Re:            " << space << prm.Re << endl;
  cout << "Plot animation?" << space << (animation ? "yes" : "no") << endl;

  double t = 0.0;
  uint numSteps = 0;

  // allocate memory
  double* u = (double*)calloc(prm.size_ghost, sizeof(double));      // x component of velocity, initialized to 0 by calloc
  double* u1 = (double*)calloc(prm.size_ghost, sizeof(double));     // x component of velocity, initialized to 0 by calloc
  double* v = (double*)calloc(prm.size_ghost, sizeof(double));      // y component of velocity, initialized to 0 by calloc
  double* v1 = (double*)calloc(prm.size_ghost, sizeof(double));     // y component of velocity, initialized to 0 by calloc
  double* p = (double*)calloc(prm.size_ghost, sizeof(double));      // pressure, initialized to 0 by calloc
  double* adv_u = (double*)calloc(prm.size_ghost, sizeof(double));  // advection term of u
  double* adv_v = (double*)calloc(prm.size_ghost, sizeof(double));  // advection term of v
  double* div = (double*)malloc(prm.size * sizeof(double));         // auxiliary variable

  if (animation) {
    file_output_anim.open(filename_solution);
    write_sol(file_output_anim, u, v, 0, prm);
  }
  END_TIMER();
  ADD_TIME_TO(total_files);

  START_TIMER();
  double lap;
  vector<Trip> coeffs;
  buildLaplaceMatrix(coeffs, prm);

  SpMat A(prm.nx * prm.ny, prm.nx * prm.ny);
  A.setFromTriplets(coeffs.begin(), coeffs.end());

  SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
  VectorXd p_solved;
  END_TIMER();
  ADD_TIME_TO(total_build_laplace);

  while (t < T) {
    START_TIMER();
    // ---------- advection term ------------
    Semilag2(u, v, u, adv_u, prm);
    Semilag2(u, v, v, adv_v, prm);

    // ---------- diffusion term ------------
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        lap = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) / (prm.dx * prm.dx) +
              (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) / (prm.dy * prm.dy);
        u1[i * prm.NY + j] = U(i, j) + prm.dt * lap / prm.Re;
        lap = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) / (prm.dx * prm.dx) +
              (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) / (prm.dy * prm.dy);
        v1[i * prm.NY + j] = V(i, j) + prm.dt * lap / prm.Re;
      }
    }

    // update velocity ghost points
    BC_velocity(u1, v1, prm);

    // ---------- pressure term ------------
    // compute the divergence of the velocity field (we omit multiplication by dt, because later on we will divide by dt) and store it in div
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        div[(i - 1) * prm.ny + (j - 1)] = (u1[(i + 1) * prm.NY + j] - u1[(i - 1) * prm.NY + j]) / (2 * prm.dx) +
                                          (v1[i * prm.NY + j + 1] - v1[i * prm.NY + j - 1]) / (2 * prm.dy);
      }
    }
    END_TIMER();
    ADD_TIME_TO(total_others);

    START_TIMER();
    // solve the pressure Poisson equation
    p_solved = chol.solve(VectorXd::Map(div, prm.nx * prm.ny));

    // convert p_solved back to pointer
    memcpy(p, p_solved.data(), prm.size * sizeof(double));
    END_TIMER();
    ADD_TIME_TO(total_solve_laplace);

    START_TIMER();
    // update pressure ghost points
    BC_pressure(p, prm);

    // ---------- projection step ------------
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        U(i, j) = u1[i * prm.NY + j] - (P(i + 1, j) - P(i - 1, j)) / (2 * prm.dx);
        V(i, j) = v1[i * prm.NY + j] - (P(i, j + 1) - P(i, j - 1)) / (2 * prm.dy);
      }
    }

    t += prm.dt;
    END_TIMER();
    ADD_TIME_TO(total_others);

    START_TIMER();
    if (WRITE_ANIM()) {
    }
    END_TIMER();
    ADD_TIME_TO(total_files);

    numSteps++;
  }

  START_TIMER();
  if (animation) write_tmp_file(prm.LX, prm.LY);
  END_TIMER();
  ADD_TIME_TO(total_files);

  print("Total time for files: ", total_files);
  print("Total time for build laplace: ", total_build_laplace);
  print("Total time for solve laplace: ", total_solve_laplace);
  print("Total time for rest: ", total_others);

  free(u);
  free(v);
  free(u1);
  free(v1);
  free(p);
  free(adv_u);
  free(adv_v);
  free(div);

  return 0;
}
