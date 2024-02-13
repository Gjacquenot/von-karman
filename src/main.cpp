#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../include/NS.hpp"
#include "../include/misc.hpp"
#include "../include/write.hpp"

using namespace Eigen;

#define EPS 1e-8

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

#define WRITE_ANIM() ((animation) && (numSteps % freq_write_anim == 0))

int main(void) {
  const string filename_solution = "data/u_solution.txt";  // name of the output file to write the solution
  const string filename_input = "data/input.txt";          // name of the input file to read the parameters
  const string space = "    ";                             // space to print
  const uint per = 10;                                     // progress percentage interval to print (each count%)
  uint count = per;
  Prm prm;
  double T;              // final time of the simulation
  bool animation;        // flag to enable or disable the animation
  uint freq_write_anim;  // frequency to write the animation file
  int64_t total_files = 0, total_advection = 0, total_diffusion = 0,
          total_others = 0, total_build_laplace = 0, total_solve_laplace = 0;  // variables to measure the time
  chrono::steady_clock::time_point begin, end;                                 // variables to measure the time

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
    file_input >> tmp >> prm.L;
    file_input >> tmp >> prm.nu;
    file_input >> tmp >> animation;
    file_input >> tmp >> freq_write_anim;
  }
  file_input.close();

  prm.NX = prm.nx + 2;
  prm.NY = prm.ny + 2;
  prm.dx = prm.LX / prm.nx;
  prm.dy = prm.LY / prm.ny;
  prm.NXNY = (uint)(prm.NX * prm.NY);
  prm.nxny = (uint)(prm.nx * prm.ny);
  prm.U = prm.Re * prm.nu / prm.L;

  // ------------- Plot setup ----------------
  cout << "dx:            " << space << prm.dx << endl;
  cout << "dy:            " << space << prm.dy << endl;
  cout << "dt:            " << space << prm.dt << endl;
  cout << "T:             " << space << T << endl;
  cout << "Re:            " << space << prm.Re << endl;
  cout << "U:             " << space << prm.U << endl;
  cout << "Plot animation?" << space << (animation ? "yes" : "no") << endl;

  double t = 0.0;
  uint numSteps = 0;
  double fraction_completed = T / 100.;  // fraction of the integration time to print

  // allocate memory
  double* u = (double*)calloc(prm.NXNY, sizeof(double));   // x component of velocity, initialized to 0 by calloc
  double* u1 = (double*)calloc(prm.NXNY, sizeof(double));  // x component of velocity, initialized to 0 by calloc
  double* v = (double*)calloc(prm.NXNY, sizeof(double));   // y component of velocity, initialized to 0 by calloc
  double* v1 = (double*)calloc(prm.NXNY, sizeof(double));  // y component of velocity, initialized to 0 by calloc
  double* p = (double*)calloc(prm.NXNY, sizeof(double));   // pressure, initialized to 0 by calloc

  // DO NOT ACCESS THE GHOST POINTS OF THESE ARRAYS
  double* adv_u = (double*)malloc(prm.NXNY * sizeof(double));  // advection term of u
  double* adv_v = (double*)malloc(prm.NXNY * sizeof(double));  // advection term of v

  if (animation) {
    file_output_anim.open(filename_solution);
    write_sol(file_output_anim, u, v, 0, prm, false);
  }
  END_TIMER();
  ADD_TIME_TO(total_files);

  START_TIMER();
  double lap, mean;
  vector<Trip> coeffs;
  buildLaplaceMatrix(coeffs, prm);

  SpMat A(prm.nx * prm.ny, prm.nx * prm.ny);
  VectorXd div(prm.nx * prm.ny);
  // fill b with 1,2,3,4,5,6
  // for (int i = 0; i < prm.nx * prm.ny; i++) b(i) = i + 1;
  A.setFromTriplets(coeffs.begin(), coeffs.end());
  // plot the matrix A (it does a nice plotof the matrix, don't worry about it for now)
  // file_output_anim << A << endl;
  // file_output_anim << "RHS" << endl;
  // file_output_anim << b << endl;

  // SimplicialCholesky<SpMat> chol(A);  // performs a Cholesky factorization of A
  SimplicialLDLT<SpMat> chol;  // the matrix has to be symmetric for this to work as expected
  // SparseLU<SpMat> chol;
  // ConjugateGradient<SpMat, Lower | Upper, IncompleteLUT<double>> chol;
  // chol.preconditioner().setDroptol(0.001);
  chol.compute(A);  // performs a Cholesky factorization of A
  if (chol.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }
  // print L
  // file_output_anim << "L" << endl;
  // file_output_anim << chol.matrixL() << endl;
  // // print U
  // file_output_anim << "U" << endl;
  // file_output_anim << chol.matrixU() << endl;

  // multiply L and U to get A
  // file_output_anim << "L * U" << endl;
  // SpMat result(prm.nx * prm.ny, prm.nx * prm.ny);
  // result = chol.matrixL() * chol.matrixU();
  // file_output_anim << result << endl;
  // for (int k = 0; k < A.outerSize(); ++k) {
  //   for (SpMat::InnerIterator it(A, k); it; ++it) {
  //     if (it.row() >= it.col()) {
  //       file_output_anim << it.value() << " ";
  //     }
  //   }
  //   file_output_anim << endl;
  // }

  // convert chol.matrixL() to a SpMat
  // convert chol.matrixL() to a dense matrix

  // complete form info in firefox
  // VectorXd XX = chol.solve(b);  // use the factorization to solve for the given right hand side
  // if (chol.info() != Success) {
  //   cout << "Cholesky solve failed" << endl;
  //   return 1;
  // }
  // // print x
  // file_output_anim << "Solution" << endl;
  // for (int i = 0; i < prm.nx * prm.ny; i++) {
  //   file_output_anim << XX(i) << " ";
  // }
  VectorXd p_solved;
  END_TIMER();
  ADD_TIME_TO(total_build_laplace);
  numSteps++;
  uint ITER = 1;
  while (t < T - EPS) {
    // ---------- advection term ------------
    START_TIMER();
    Semilag2(u, v, u, adv_u, prm);
    Semilag2(u, v, v, adv_v, prm);
    END_TIMER();
    ADD_TIME_TO(total_advection);
    // if (numSteps == ITER) {
    //   file_output_anim << "Advection" << endl;
    //   write_sol(file_output_anim, adv_u, adv_v, 500000000000000, prm, true);
    // }
    // ---------- diffusion term ------------
    START_TIMER();
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        lap = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) / (prm.dx * prm.dx) +
              (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) / (prm.dy * prm.dy);
        U1(i, j) = ADV_U(i, j) + prm.dt * lap / prm.Re;
        lap = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) / (prm.dx * prm.dx) +
              (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) / (prm.dy * prm.dy);
        V1(i, j) = ADV_V(i, j) + prm.dt * lap / prm.Re;
      }
    }

    // update velocity ghost points
    BC_velocity(u1, v1, prm);
    END_TIMER();
    ADD_TIME_TO(total_diffusion);

    // if (numSteps == ITER) {
    //   file_output_anim << "Diffusion" << endl;
    //   write_sol(file_output_anim, u1, v1, 500000000000000, prm, true);
    // }

    // ---------- pressure term ------------
    START_TIMER();
    // compute the (minus) divergence of the velocity field (we omit multiplication by dt, because later on we will divide by dt) and store it in div
    mean = 0;
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        DIV(i - 1, j - 1) = -(U1(i + 1, j) - U1(i - 1, j)) / (2 * prm.dx) -
                            (V1(i, j + 1) - V1(i, j - 1)) / (2 * prm.dy);
        mean += DIV(i - 1, j - 1);
      }
    }
    mean /= prm.nx * prm.ny;
    // if (numSteps == ITER) {
    //   file_output_anim << "Divergence" << endl;
    //   for (int i = 0; i < prm.nx; i++) {
    //     for (int j = 0; j < prm.ny; j++) {
    //       file_output_anim << div[i * prm.ny + j] << " ";
    //     }
    //     file_output_anim << endl;
    //   }
    // }
    // we remove the mean to ensure that the Poisson equation has a solution.
    for (int i = 0; i < prm.nx * prm.ny; i++) div[i] -= mean;

    // print div
    // if (numSteps == ITER) {
    // file_output_anim << "(minus) Divergence" << endl;
    // for (int i = 0; i < prm.nx; i++) {
    //   for (int j = 0; j < prm.ny; j++) {
    //     file_output_anim << DIV(i, j) << " ";
    //   }
    //   file_output_anim << endl;
    // }

    END_TIMER();
    ADD_TIME_TO(total_others);

    START_TIMER();
    // solve the pressure Poisson equation

    p_solved = chol.solve(div);
    if (chol.info() != Success) {
      cout << "Cholesky solve failed" << endl;
      return 1;
    }
    // print p_solved
    // file_output_anim << "Solution" << endl;
    // file_output_anim << p_solved << endl;

    // print p_solved
    // file_output_anim << "Pressure" << endl;
    // convert p_solved back to pointer
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        P(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
      }
    }
    // write_sol(file_output_anim, p, p, 500000000000000, prm, false);

    END_TIMER();
    ADD_TIME_TO(total_solve_laplace);

    START_TIMER();
    // update pressure ghost points
    BC_pressure(p, prm);
    // if (numSteps == ITER) {
    //   file_output_anim << "Diffusion2" << endl;
    //   write_sol(file_output_anim, u1, v1, 500000000000000, prm, true);
    // }
    // ---------- projection step ------------
    memcpy(u, u1, prm.NXNY * sizeof(double));
    memcpy(v, v1, prm.NXNY * sizeof(double));
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        U(i, j) -= (P(i + 1, j) - P(i - 1, j)) / (2 * prm.dx);
        V(i, j) -= (P(i, j + 1) - P(i, j - 1)) / (2 * prm.dy);
      }
    }
    // if (numSteps == ITER) {
    // file_output_anim << "Final" << endl;
    // write_sol(file_output_anim, u, v, 500000000000000, prm, true);
    // }
    t += prm.dt;
    END_TIMER();
    ADD_TIME_TO(total_others);

    START_TIMER();
    if (WRITE_ANIM()) {
      write_sol(file_output_anim, u, v, t, prm, false);
    }

    if (t > fraction_completed * count - EPS) {
      cout << count << "%" << endl;
      count += per;
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
  print("Total time for advection: ", total_advection);
  print("Total time for diffusion: ", total_diffusion);
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

  return 0;
}
