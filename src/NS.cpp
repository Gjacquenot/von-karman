#include "../include/NS.hpp"

#include <cstring>

// 1st order semi-Lagrangian advection
void Semilag(double* u, double* v, double* q, Prm prm, int sign, Object& obstacle) {
  int sign_u, sign_v;
  double a, b;
  double* aux = (double*)calloc(prm.NXNY, sizeof(double));
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      if (prm.obstacle_ON && obstacle.IsInside[i * prm.NY + j]) {
        continue;
      }
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
void Semilag2(double* u, double* v, double* q0, double* q1, Prm prm, Object& obstacle) {
  memcpy(q1, q0, prm.NXNY * sizeof(double));
  Semilag(u, v, q1, prm, 1, obstacle);
  Semilag(u, v, q1, prm, -1, obstacle);
  for (int i = 0; i < prm.NX * prm.NY; i++)
    q1[i] = q0[i] + (q0[i] - q1[i]) / 2;
  Semilag(u, v, q1, prm, 1, obstacle);
}

// void upwind(double* u, double* v, double* q, Prm prm, int sign, Object& obstacle) {
//   int sign_u, sign_v;
//   double* aux = (double*)calloc(prm.NXNY, sizeof(double));
//   double q_x, q_y;
//   for (int i = 1; i < prm.NX - 1; i++) {
//     for (int j = 1; j < prm.NY - 1; j++) {
//       if (prm.obstacle_ON && obstacle.IsInside[i * prm.NY + j]) {
//         continue;
//       }
//       if (sign * U(i, j) > 0) {
//         q_x = q[i * prm.NY + j] - q[(i - 1) * prm.NY + j];
//       } else {
//         q_x = q[(i + 1) * prm.NY + j] - q[i * prm.NY + j];
//       }
//       if (sign * V(i, j) > 0) {
//         q_y = q[i * prm.NY + j] - q[i * prm.NY + j - 1];
//       } else {
//         q_y = q[i * prm.NY + j + 1] - q[i * prm.NY + j];
//       }
//       aux[i * prm.NY + j] = q[i * prm.NY + j] - sign * prm.dt * (U(i, j) * q_x / prm.dx + V(i, j) * q_y / prm.dy);
//     }
//   }

//   free(aux);
// }

// void upwind2(double* u, double* v, double* q0, double* q1, Prm prm, Object& obstacle) {
//   memcpy(q1, q0, prm.NXNY * sizeof(double));
//   upwind(u, v, q1, prm, 1, obstacle);
//   upwind(u, v, q1, prm, -1, obstacle);
//   for (int i = 0; i < prm.NX * prm.NY; i++)
//     q1[i] = q0[i] + (q0[i] - q1[i]) / 2;
//   upwind(u, v, q1, prm, 1, obstacle);
// }

void BC_velocity(double* u, double* v, Prm prm, Object& obstacle) {
  // outer boundaries
  for (int i = 0; i < prm.NX; i++) {
    // bottom boundary: \partial_y u = 0, v = 0
    // U(i, 0) = -U(i, 1);
    U(i, 0) = U(i, 1);
    V(i, 0) = -V(i, 1);
    // top boundary: \partial_y u = 0, v = 0
    // U(i, prm.NY - 1) = -U(i, prm.NY - 2);
    U(i, prm.NY - 1) = U(i, prm.NY - 2);
    V(i, prm.NY - 1) = -V(i, prm.NY - 2);
  }

  for (int j = 0; j < prm.NY; j++) {
    // left boundary: u = U, v = 0
    // first order extrapolation
    U(0, j) = 2 - U(1, j);
    V(0, j) = -V(1, j);
    // right boundary:
    U(prm.NX - 1, j) = U(prm.NX - 2, j);
    V(prm.NX - 1, j) = V(prm.NX - 2, j);
    // U(prm.NX - 1, j) = 2 * U(prm.NX - 2, j) - U(prm.NX - 3, j);
    // V(prm.NX - 1, j) = 2 * V(prm.NX - 2, j) - V(prm.NX - 3, j);
  }

  // for (int j = 1; j < prm.NY - 1; j++) {
  //   // left boundary: u = U, v = 0
  //   U(0, j) = -U(1, j);
  //   V(0, j) = -V(1, j);
  //   // right boundary:
  //   U(prm.NX - 1, j) = -U(prm.NX - 2, j);
  //   V(prm.NX - 1, j) = -V(prm.NX - 2, j);
  // }
  // for (int i = 1; i < prm.NX - 1; i++) {
  //   // bottom boundary: \partial_y u = 0, v = 0
  //   U(i, 0) = U(i, 1);
  //   V(i, 0) = -V(i, 1);
  //   // top boundary: \partial_y u = 0, v = 0
  //   U(i, prm.NY - 1) = 2 - U(i, prm.NY - 2);
  //   V(i, prm.NY - 1) = -V(i, prm.NY - 2);
  // }

  // if (prm.obstacle_ON) {
  //   // immersed boundary
  //   Eigen::MatrixXd A(6, 6);
  //   // A will be the matrix of the system Ax = b and it has the form
  //   // | 1 x_1 y_1 x_1^2 x_1*y_1 y_1^2 |
  //   // | 1 x_2 y_2 x_2^2 x_2*y_2 y_2^2 |
  //   // | 1 x_3 y_3 x_3^2 x_3*y_3 y_3^2 |
  //   // | 1 x_4 y_4 x_4^2 x_4*y_4 y_4^2 |
  //   // | 1 x_5 y_5 x_5^2 x_5*y_5 y_5^2 |
  //   // | 1 x_6 y_6 x_6^2 x_6*y_6 y_6^2 |
  //   // where x_i and y_i are the x and y coordinates of the interpolating points
  //   Eigen::VectorXd b(6), sol(6);
  //   // b will be the vector of the system Ax = b and it has the form
  //   // | phi_1 |
  //   // | phi_2 |
  //   // | phi_3 |
  //   // | phi_4 |
  //   // | phi_5 |
  //   // | phi_6 |
  //   // where phi_i is the value of each velocity component (separately) at the interpolating points
  //   double phi_m;  // value of the velocity component at the mirror point
  //   for (int i = 0; i < obstacle.count_ghost; i++) {
  //     A << 1, x(p(0).i), y(p(0).j), x(p(0).i) * x(p(0).i), x(p(0).i) * y(p(0).j), y(p(0).j) * y(p(0).j),
  //         1, x(p(1).i), y(p(1).j), x(p(1).i) * x(p(1).i), x(p(1).i) * y(p(1).j), y(p(1).j) * y(p(1).j),
  //         1, x(p(2).i), y(p(2).j), x(p(2).i) * x(p(2).i), x(p(2).i) * y(p(2).j), y(p(2).j) * y(p(2).j),
  //         1, x(p(3).i), y(p(3).j), x(p(3).i) * x(p(3).i), x(p(3).i) * y(p(3).j), y(p(3).j) * y(p(3).j),
  //         1, x(p(4).i), y(p(4).j), x(p(4).i) * x(p(4).i), x(p(4).i) * y(p(4).j), y(p(4).j) * y(p(4).j),
  //         1, x(p(5).i), y(p(5).j), x(p(5).i) * x(p(5).i), x(p(5).i) * y(p(5).j), y(p(5).j) * y(p(5).j);
  //     b << U(p(0).i, p(0).j), U(p(1).i, p(1).j), U(p(2).i, p(2).j), U(p(3).i, p(3).j), U(p(4).i, p(4).j), U(p(5).i, p(5).j);
  //     sol = A.colPivHouseholderQr().solve(b);
  //     phi_m = sol(0) + sol(1) * m_x + sol(2) * m_y + sol(3) * m_x * m_x + sol(4) * m_x * m_y + sol(5) * m_y * m_y;
  //     U(obstacle.GhostPoints[i].i, obstacle.GhostPoints[i].j) = -phi_m;  // u = 0 at the boundary

  //     b << V(p(0).i, p(0).j), V(p(1).i, p(1).j), V(p(2).i, p(2).j), V(p(3).i, p(3).j), V(p(4).i, p(4).j), V(p(5).i, p(5).j);
  //     sol = A.colPivHouseholderQr().solve(b);
  //     phi_m = sol(0) + sol(1) * m_x + sol(2) * m_y + sol(3) * m_x * m_x + sol(4) * m_x * m_y + sol(5) * m_y * m_y;
  //     V(obstacle.GhostPoints[i].i, obstacle.GhostPoints[i].j) = -phi_m;  // v = 0 at the boundary
  //   }
  // }

  // immersed boundary

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

void BC_pressure(double* p, Prm prm, Object& obstacle) {
  for (int i = 0; i < prm.NX; i++) {
    // bottom boundary: \partial_y p = 0
    P(i, 0) = P(i, 1);
    // top boundary: \partial_y p = 0
    P(i, prm.NY - 1) = P(i, prm.NY - 2);
  }
  for (int j = 0; j < prm.NY; j++) {
    // left boundary: \partial_x p = 0
    P(0, j) = P(1, j);
    // right boundary: p = 0
    // P(prm.NX - 1, j) = P(prm.NX - 2, j);  // CFD lecture example
    P(prm.NX - 1, j) = -P(prm.NX - 2, j);  // mine
  }
  // if (prm.obstacle_ON) {
  //   // immersed boundary
  //   Eigen::MatrixXd A(6, 6);
  //   // A will be the matrix of the system Ax = b and it has the form
  //   // | 1 x_1 y_1 x_1^2 x_1*y_1 y_1^2 |
  //   // | 1 x_2 y_2 x_2^2 x_2*y_2 y_2^2 |
  //   // | 1 x_3 y_3 x_3^2 x_3*y_3 y_3^2 |
  //   // | 1 x_4 y_4 x_4^2 x_4*y_4 y_4^2 |
  //   // | 1 x_5 y_5 x_5^2 x_5*y_5 y_5^2 |
  //   // | 1 x_6 y_6 x_6^2 x_6*y_6 y_6^2 |
  //   // where x_i and y_i are the x and y coordinates of the interpolating points
  //   Eigen::VectorXd b(6), sol(6);
  //   // b will be the vector of the system Ax = b and it has the form
  //   // | phi_1 |
  //   // | phi_2 |
  //   // | phi_3 |
  //   // | phi_4 |
  //   // | phi_5 |
  //   // | phi_6 |
  //   // where phi_i is the value of the pressure at the interpolating points
  //   double P_m;  // value of the pressure at the mirror point
  //   for (int i = 0; i < obstacle.count_ghost; i++) {
  //     A << 1, x(p(0).i), y(p(0).j), x(p(0).i) * x(p(0).i), x(p(0).i) * y(p(0).j), y(p(0).j) * y(p(0).j),
  //         1, x(p(1).i), y(p(1).j), x(p(1).i) * x(p(1).i), x(p(1).i) * y(p(1).j), y(p(1).j) * y(p(1).j),
  //         1, x(p(2).i), y(p(2).j), x(p(2).i) * x(p(2).i), x(p(2).i) * y(p(2).j), y(p(2).j) * y(p(2).j),
  //         1, x(p(3).i), y(p(3).j), x(p(3).i) * x(p(3).i), x(p(3).i) * y(p(3).j), y(p(3).j) * y(p(3).j),
  //         1, x(p(4).i), y(p(4).j), x(p(4).i) * x(p(4).i), x(p(4).i) * y(p(4).j), y(p(4).j) * y(p(4).j),
  //         1, x(p(5).i), y(p(5).j), x(p(5).i) * x(p(5).i), x(p(5).i) * y(p(5).j), y(p(5).j) * y(p(5).j);
  //     b << P(p(0).i, p(0).j), P(p(1).i, p(1).j), P(p(2).i, p(2).j), P(p(3).i, p(3).j), P(p(4).i, p(4).j), P(p(5).i, p(5).j);
  //     sol = A.colPivHouseholderQr().solve(b);
  //     P_m = sol(0) + sol(1) * m_x + sol(2) * m_y + sol(3) * m_x * m_x + sol(4) * m_x * m_y + sol(5) * m_y * m_y;
  //     P(obstacle.GhostPoints[i].i, obstacle.GhostPoints[i].j) = P_m;  // partial_n p = 0 at the boundary
  //   }
  // }
}

double interpolate(double x, double y, double* phi, Prm prm) {
  // find the indices of the cell containing the point (x, y)
  // find I, J such that x(I) < x < x(I+1) and y(J) < y < y(J+1)

  int I = 0, J = 0;  // minimum value for I and J (the maximum value is I+2 and J+2)
  while (x(I + 1) < x + TOL) I++;
  while (y(J + 1) < y + TOL) J++;

  // interpolate (bi-linear interpolation)
  double phi_00 = phi[I * prm.NY + J] * (x(I + 1) - x) * (y(J + 1) - y);
  double phi_01 = phi[I * prm.NY + J + 1] * (x(I + 1) - x) * (y - y(J));
  double phi_10 = phi[(I + 1) * prm.NY + J] * (x - x(I)) * (y(J + 1) - y);
  double phi_11 = phi[(I + 1) * prm.NY + J + 1] * (x - x(I)) * (y - y(J));
  return (phi_00 + phi_01 + phi_10 + phi_11) / (prm.dx * prm.dy);
}

void set_vorticity(double* u, double* v, double* w, Prm prm) {
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      W(i, j) = (V(i + 1, j) - V(i - 1, j)) / (2 * prm.dx) - (U(i, j + 1) - U(i, j - 1)) / (2 * prm.dy);
    }
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
  double dx_2 = 1. / (prm.dx * prm.dx);
  double dy_2 = 1. / (prm.dy * prm.dy);
  double diagX, diagY;

  for (int i = 0; i < dim; i++) {
    // diagonal
    diagX = -2. * dx_2;
    diagY = -2. * dy_2;
    if (i % prm.ny == 0 || i % prm.ny == prm.ny - 1) diagY = -dy_2;
    // for full neumann BC (test case, remebr that in this case the matrix is singular and the cholensky decomposition will fail, you have to change for example the component (0, 1) to 0, and, to make the matrix symmetric, the component (1, 0) to 0)
    // if (i < prm.ny || i >= dim - prm.ny) diagX = -dx_2;
    // for mixed dirichlet-neumann BC (our case)
    if (i < prm.ny) diagX = -dx_2;
    if (i >= dim - prm.ny) diagX = -3. * dx_2;
    coeffs.push_back(Trip(i, i, -diagX - diagY));

    // Dxx part (secondary diagonals)
    if (i < dim - prm.ny) {
      coeffs.push_back(Trip(i, i + prm.ny, -dx_2));  // upper diagonal
      coeffs.push_back(Trip(i + prm.ny, i, -dx_2));  // lower diagonal
    }

    // Dyy part (secondary diagonals)
    if ((i + 1) % prm.ny != 0) {
      coeffs.push_back(Trip(i, i + 1, -dy_2));  // upper diagonal
      coeffs.push_back(Trip(i + 1, i, -dy_2));  // lower diagonal
    }
  }
  // coeffs.push_back(Trip(0, 1, dy_2));
  // coeffs.push_back(Trip(1, 0, dy_2));
}