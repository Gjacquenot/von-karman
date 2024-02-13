#include <chrono>
#include <cmath>
#include <cstring>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace Eigen;
using namespace std;
typedef Eigen::SparseMatrix<double> SpMat;  // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;        // declares a triplet type of double

int main(void) {
  vector<Trip> coeffs;

  int nx = 8;
  int ny = 4;
  int n = nx * ny;
  SpMat A(n, n);
  VectorXd b(n);
  for (int i = 0; i < n; i++) b(i) = 0;
  b(ny - 1) = -0.000256;
  b(n - 1) = 0.000256;
  double dx = 1. / nx;
  double dy = 1. / ny;
  double dx2 = dx * dx;
  double dy2 = dy * dy;
  double diagX, diagY;

  coeffs.reserve(5 * n + 2);
  for (int i = 0; i < n; i++) {
    // diagonal
    diagX = -2. / dx2;
    diagY = -2. / dy2;
    if (i % ny == 0 || i % ny == ny - 1) diagY = -1. / dy2;
    // for full neumann BC
    if (i < ny || i >= n - ny) diagX = -1. / dx2;
    // for mixed dirichlet-neumann BC (our case)
    // if (i < ny) diagX = -1. / dx2;
    // if (i >= dim - ny) diagX = -3. / dx2;
    coeffs.push_back(Trip(i, i, -diagX - diagY));

    // Dxx part (secondary diagonals)
    if (i < n - ny) {
      coeffs.push_back(Trip(i, i + ny, -1. / dx2));  // upper diagonal
      coeffs.push_back(Trip(i + ny, i, -1. / dx2));  // lower diagonal
    }

    // Dyy part (secondary diagonals)
    if ((i + 1) % ny != 0) {
      coeffs.push_back(Trip(i, i + 1, -1. / dy2));  // upper diagonal
      coeffs.push_back(Trip(i + 1, i, -1. / dy2));  // lower diagonal
    }
  }
  coeffs.push_back(Trip(0, 1, 1. / dy2));
  coeffs.push_back(Trip(1, 0, 1. / dy2));

  A.setFromTriplets(coeffs.begin(), coeffs.end());
  // plot the matrix A (it does a nice plotof the matrix, don't worry about it for now)
  cout << A << endl;
  // file_output_anim << "RHS" << endl;
  cout << b << endl;

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
  // for (int i = 0; i < n; i++) {
  //   file_output_anim << XX(i) << " ";
  // }
  VectorXd p_solved;
  p_solved = chol.solve(b);
  if (chol.info() != Success) {
    cout << "Cholesky solve failed" << endl;
    return 1;
  }
  cout << "Solution:" << endl;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      cout << p_solved(i * ny + j) << " ";
    }
    cout << endl;
  }

  // cout << "#iterations:     " << chol.iterations() << std::endl;
  // cout << "estimated error: " << chol.error() << std::endl;

  return 0;
}
