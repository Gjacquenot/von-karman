#include "../include/NS.hpp"

#include <cstring>

// 1st order semi-Lagrangian advection
void Semilag(double* u, double* v, double* q, Prm prm, int sign) {
  int sign_u, sign_v;
  double a, b;
  double* aux = (double*)calloc(prm.NXNY, sizeof(double));
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      // if (U(i, j) > 0) {
      //   a = 1 - U(i, j) * prm.dt / prm.dx;
      //   if (V(i, j) > 0) {
      //     b = 1 - V(i, j) * prm.dt / prm.dy;
      //     C00 = a * b * q[i * prm.NY + j];
      //     C01 = a * (1 - b) * q[i * prm.NY + j - 1];
      //     C10 = (1 - a) * b * q[(i - 1) * prm.NY + j];
      //     C11 = (1 - a) * (1 - b) * q[(i - 1) * prm.NY + j - 1];
      //   } else {
      //     b = 1 + V(i, j) * prm.dt / prm.dy;
      //     C00 = a * b * q[i * prm.NY + j];
      //     C01 = a * (1 - b) * q[i * prm.NY + j + 1];
      //     C10 = (1 - a) * b * q[(i - 1) * prm.NY + j];
      //     C11 = (1 - a) * (1 - b) * q[(i - 1) * prm.NY + j + 1];
      //   }
      // } else {
      //   a = 1 + U(i, j) * prm.dt / prm.dx;
      //   if (V(i, j) > 0) {
      //     b = 1 - V(i, j) * prm.dt / prm.dy;
      //     C00 = a * b * q[i * prm.NY + j];
      //     C01 = a * (1 - b) * q[i * prm.NY + j - 1];
      //     C10 = (1 - a) * b * q[(i + 1) * prm.NY + j];
      //     C11 = (1 - a) * (1 - b) * q[(i + 1) * prm.NY + j - 1];
      //   } else {
      //     b = 1 + V(i, j) * prm.dt / prm.dy;
      //     C00 = a * b * q[i * prm.NY + j];
      //     C01 = a * (1 - b) * q[i * prm.NY + j + 1];
      //     C10 = (1 - a) * b * q[(i + 1) * prm.NY + j];
      //     C11 = (1 - a) * (1 - b) * q[(i + 1) * prm.NY + j + 1];
      //   }
      // }
      // ADVq[i * prm.NY + j] = C00 + C01 + C10 + C11;

      if (sign * U(i, j) > 0) {
        a = 1 - sign * U(i, j) * prm.dt / prm.dx;
        sign_u = 1;
      } else {
        a = 1 + sign * U(i, j) * prm.dt / prm.dx;
        sign_u = -1;
      }
      if (sign * V(i, j) > 0) {
        b = 1 - sign * V(i, j) * prm.dt / prm.dy;
        sign_v = 1;
      } else {
        b = 1 + sign * V(i, j) * prm.dt / prm.dy;
        sign_v = -1;
      }
      aux[i * prm.NY + j] = a * b * q[i * prm.NY + j] +
                            (1 - a) * b * q[(i - sign_u) * prm.NY + j] +
                            a * (1 - b) * q[i * prm.NY + j - sign_v] +
                            (1 - a) * (1 - b) * q[(i - sign_u) * prm.NY + j - sign_v];
    }
  }
  memcpy(q, aux, prm.NXNY * sizeof(double));
  free(aux);
}

// 2nd order semi-Lagrangian advection
void Semilag2(double* u, double* v, double* q0, double* q1, Prm prm) {
  memcpy(q1, q0, prm.NXNY * sizeof(double));
  Semilag(u, v, q1, prm, 1);
  Semilag(u, v, q1, prm, -1);
  for (int i = 0; i < prm.NX * prm.NY; i++)
    q1[i] = q0[i] + (q0[i] - q1[i]) / 2;
  Semilag(u, v, q1, prm, 1);
}

void BC_velocity(double* u, double* v, Prm prm) {
  // mine
  for (int j = 1; j < prm.NY - 1; j++) {
    // left boundary: u = U, v = 0, \partial_x p = 0
    U(0, j) = 2 * prm.U - U(1, j);
    V(0, j) = -V(1, j);
    // right boundary: \partial_x u = 0, \partial_x v = 0, p = 0
    U(prm.NX - 1, j) = U(prm.NX - 2, j);
    V(prm.NX - 1, j) = V(prm.NX - 2, j);
  }
  for (int i = 1; i < prm.NX - 1; i++) {
    // bottom boundary: \partial_y u = 0, v = 0, \partial_y p = 0
    U(i, 0) = U(i, 1);
    V(i, 0) = -V(i, 1);
    // top boundary: \partial_y u = 0, v = 0, \partial_y p = 0
    U(i, prm.NY - 1) = U(i, prm.NY - 2);
    V(i, prm.NY - 1) = -V(i, prm.NY - 2);
  }

  // CFD lecture example
  // for (int j = 1; j < prm.NY - 1; j++) {
  //   // left boundary: u = U, v = 0, \partial_x p = 0
  //   U(0, j) = -U(1, j);
  //   V(0, j) = -V(1, j);
  //   // right boundary: \partial_x u = 0, \partial_x v = 0, p = 0
  //   U(prm.NX - 1, j) = -U(prm.NX - 2, j);
  //   V(prm.NX - 1, j) = -V(prm.NX - 2, j);
  // }
  // for (int i = 1; i < prm.NX - 1; i++) {
  //   // bottom boundary: \partial_y u = 0, v = 0, \partial_y p = 0
  //   U(i, 0) = -U(i, 1);
  //   V(i, 0) = -V(i, 1);
  //   // top boundary: \partial_y u = 0, v = 0, \partial_y p = 0
  //   U(i, prm.NY - 1) = 2 - U(i, prm.NY - 2);
  //   V(i, prm.NY - 1) = -V(i, prm.NY - 2);
  // }
}

void BC_pressure(double* p, Prm prm) {
  for (int j = 1; j < prm.NY - 1; j++) {
    // left boundary: u = U, v = 0, \partial_x p = 0
    P(0, j) = P(1, j);
    // right boundary: \partial_x u = 0, \partial_x v = 0, p = 0
    // P(prm.NX - 1, j) = P(prm.NX - 2, j);  // CFD lecture example
    P(prm.NX - 1, j) = -P(prm.NX - 2, j);  // mine
  }
  for (int i = 1; i < prm.NX - 1; i++) {
    // bottom boundary: \partial_y u = 0, v = 0, \partial_y p = 0
    P(i, 0) = P(i, 1);
    // top boundary: \partial_y u = 0, v = 0, \partial_y p = 0
    P(i, prm.NY - 1) = P(i, prm.NY - 2);
  }
}

void buildLaplaceMatrix(vector<Trip>& coeffs, Prm prm) {
  // Remeber we are taking the following BC for the pressure:
  // on the top, bottom and left boundaries, the normal derivative of the pressure is 0, i.e.:  P(i, NY - 1) = P(i, NY - 2)
  //                                                                                            P(i, 0) = P(i, 1);
  //                                                                                            P(0, j) = P(1, j);
  // respectively. On the right boundary, the pressure is 0, i.e.: P(NX - 1, j) = -P(NX - 2, j).
  // Also, note that the coefficients are ordered by column, i.e. the first NY elements are the first column, the next NY elements are the second column, and so on.
  // The matrix has size (nx * ny) x (nx * ny), where nx = NX - 2 and ny = NY - 2 (i.e. the number of points in the domain minus the ghost points).
  // Matrix for the 2nd derivative in x (except for the division by dx^2)
  //  ------ ny ------
  // ||-----------------------------------------------------------------------------------------||
  // || -1           |  1           |                                                           || |
  // ||   -1         |    1         |                                                           || |
  // ||      .       |      .       |                                                           || ny
  // ||        .     |        .     |                                                           || |
  // ||          .   |          .   |                                                           || |
  // ||           -1 |            1 |                                                           || |
  // ||-----------------------------------------------------------------------------------------||
  // || 1            | -2           | 1            |                                            ||
  // ||   1          |   -2         |   1          |                                            ||
  // ||      .       |      .       |      .       |                                            ||
  // ||        .     |        .     |        .     |                                            ||
  // ||          .   |          .   |          .   |                                            ||
  // ||            1 |           -2 |            1 |                                            ||
  // ||-----------------------------------------------------------------------------------------||
  // ||              | 1            | -2           | 1            |                             ||
  // ||              |   1          |   -2         |   1          |                             ||
  // ||              |      .       |      .       |      .       |                             ||
  // ||              |        .     |        .     |        .     |                             ||
  // ||              |          .   |          .   |          .   |                             ||
  // ||              |            1 |           -2 |            1 |                             ||
  // ||-----------------------------------------------------------------------------------------||
  // ||                             |  .           |  .           |  .           |              ||
  // ||                             |    .         |    .         |    .         |              ||
  // ||                             |      .       |      .       |      .       |              ||
  // ||                             |        .     |        .     |        .     |              ||
  // ||                             |          .   |          .   |          .   |              ||
  // ||                             |            . |            . |            . |              ||
  // ||-----------------------------------------------------------------------------------------||
  // ||                                            | 1            | -2           | 1            ||
  // ||                                            |   1          |   -2         |   1          ||
  // ||                                            |      .       |      .       |      .       ||
  // ||                                            |        .     |        .     |        .     ||
  // ||                                            |          .   |          .   |          .   ||
  // ||                                            |            1 |           -2 |            1 ||
  // ||-----------------------------------------------------------------------------------------||
  // ||                                                           | 1            | -3           ||
  // ||                                                           |   1          |   -3         ||
  // ||                                                           |      .       |      .       ||
  // ||                                                           |        .     |        .     ||
  // ||                                                           |          .   |          .   ||
  // ||                                                           |            1 |           -3 ||
  // ||-----------------------------------------------------------------------------------------||

  // Matrix for the 2nd derivative in y (except for the division by dy^2)
  //  --------- ny --------
  // ||------------------------------------------------------------------------------------||
  // || -1  1              |                                                               || |
  // ||  1 -2  1           |                                                               || |
  // ||    1 -2  1         |                                                               || |
  // ||      .  .  .       |                                                               || ny
  // ||        .  .  .     |                                                               || |
  // ||          .  .  .   |                                                               || |
  // ||            1 -2  1 |                                                               || |
  // ||               1 -1 |                                                               || |
  // ||------------------------------------------------------------------------------------||
  // ||                    | -1  1              |                                          ||
  // ||                    |  1 -2  1           |                                          ||
  // ||                    |    1 -2  1         |                                          ||
  // ||                    |      .  .  .       |                                          ||
  // ||                    |        .  .  .     |                                          ||
  // ||                    |          .  .  .   |                                          ||
  // ||                    |            1 -2  1 |                                          ||
  // ||                    |               1 -1 |                                          ||
  // ||------------------------------------------------------------------------------------||
  // ||                                          |   .                |                    ||
  // ||                                          |     .              |                    ||
  // ||                                          |       .            |                    ||
  // ||                                          |         .          |                    ||
  // ||                                          |           .        |                    ||
  // ||                                          |             .      |                    ||
  // ||                                          |               .    |                    ||
  // ||                                          |                 .  |                    ||
  // ||------------------------------------------------------------------------------------||
  // ||                                                               | -1  1              ||
  // ||                                                               |  1 -2  1           ||
  // ||                                                               |    1 -2  1         ||
  // ||                                                               |      .  .  .       ||
  // ||                                                               |        .  .  .     ||
  // ||                                                               |          .  .  .   ||
  // ||                                                               |            1 -2  1 ||
  // ||                                                               |               1 -1 ||
  // ||------------------------------------------------------------------------------------||

  // number of non-zero coefficients in the Laplace matrix
  // dx -> 3 * ny * (nx - 2) + 2 * (2 * ny) = ny * (3 * nx - 2)
  // dy -> (3 * ny - 2) * nx
  // common (diag) -> nx * ny
  // total = dx + dy - common = 5 * nx * ny - 2 * ny - 2 * nx

  // we will store the minus of the coefficients, because we will use the cholensky decomposition, which requires a positive definite matrix

  // coeffs.reserve((uint)(prm.nx * prm.ny));
  // coeffs.push_back(Trip(0, 0, 32));
  // coeffs.push_back(Trip(0, 1, 9));
  // coeffs.push_back(Trip(0, 2, 16));
  // coeffs.push_back(Trip(0, 3, 23));
  // coeffs.push_back(Trip(0, 4, 30));
  // coeffs.push_back(Trip(0, 5, 37));
  // coeffs.push_back(Trip(1, 0, 9));
  // coeffs.push_back(Trip(1, 1, 46));
  // coeffs.push_back(Trip(1, 2, 23));
  // coeffs.push_back(Trip(1, 3, 30));
  // coeffs.push_back(Trip(1, 4, 37));
  // coeffs.push_back(Trip(1, 5, 44));
  // coeffs.push_back(Trip(2, 0, 16));
  // coeffs.push_back(Trip(2, 1, 23));
  // coeffs.push_back(Trip(2, 2, 60));
  // coeffs.push_back(Trip(2, 3, 37));
  // coeffs.push_back(Trip(2, 4, 44));
  // coeffs.push_back(Trip(2, 5, 51));
  // coeffs.push_back(Trip(3, 0, 23));
  // coeffs.push_back(Trip(3, 1, 30));
  // coeffs.push_back(Trip(3, 2, 37));
  // coeffs.push_back(Trip(3, 3, 74));
  // coeffs.push_back(Trip(3, 4, 51));
  // coeffs.push_back(Trip(3, 5, 58));
  // coeffs.push_back(Trip(4, 0, 30));
  // coeffs.push_back(Trip(4, 1, 37));
  // coeffs.push_back(Trip(4, 2, 44));
  // coeffs.push_back(Trip(4, 3, 51));
  // coeffs.push_back(Trip(4, 4, 88));
  // coeffs.push_back(Trip(4, 5, 65));
  // coeffs.push_back(Trip(5, 0, 37));
  // coeffs.push_back(Trip(5, 1, 44));
  // coeffs.push_back(Trip(5, 2, 51));
  // coeffs.push_back(Trip(5, 3, 58));
  // coeffs.push_back(Trip(5, 4, 65));
  // coeffs.push_back(Trip(5, 5, 102));

  coeffs.reserve((uint)(5 * prm.nx * prm.ny - 2 * (prm.ny + prm.nx)));
  int dim = prm.nx * prm.ny;
  double dx2 = prm.dx * prm.dx;
  double dy2 = prm.dy * prm.dy;
  double diagX, diagY;

  for (int i = 0; i < dim; i++) {
    // diagonal
    diagX = -2. / dx2;
    diagY = -2. / dy2;
    if (i % prm.ny == 0 || i % prm.ny == prm.ny - 1) diagY = -1. / dy2;
    // for full neumann BC (test case, remebr that in this case the matrix is singular and the cholensky decomposition will fail, you have to change for example the component (0, 1) to 0, and, to make the matrix symmetric, the component (1, 0) to 0)
    // if (i < prm.ny || i >= dim - prm.ny) diagX = -1. / dx2;
    // for mixed dirichlet-neumann BC (our case)
    if (i < prm.ny) diagX = -1. / dx2;
    if (i >= dim - prm.ny) diagX = -3. / dx2;
    coeffs.push_back(Trip(i, i, -diagX - diagY));

    // Dxx part (secondary diagonals)
    if (i < dim - prm.ny) {
      coeffs.push_back(Trip(i, i + prm.ny, -1. / dx2));  // upper diagonal
      coeffs.push_back(Trip(i + prm.ny, i, -1. / dx2));  // lower diagonal
    }

    // Dyy part (secondary diagonals)
    if ((i + 1) % prm.ny != 0) {
      coeffs.push_back(Trip(i, i + 1, -1. / dy2));  // upper diagonal
      coeffs.push_back(Trip(i + 1, i, -1. / dy2));  // lower diagonal
    }
  }
  // coeffs.push_back(Trip(0, 1, 1. / dy2));
  // coeffs.push_back(Trip(1, 0, 1. / dy2));
}