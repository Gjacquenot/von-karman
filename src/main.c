#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define U(i, j) u[(i) * prm.NY + (j)]
#define V(i, j) v[(i) * prm.NY + (j)]
#define P(i, j) p[(i) * prm.NY + (j)]
#define x(i) (-0.5 * prm.dx + (i) * prm.dx)
#define y(j) (-0.5 * prm.dy + (j) * prm.dy)

#define TOL 1.e-8
#define SIGN(x) ((x) > TOL ? 1 : ((x) < -TOL ? -1 : 0))
#define ABS(x) ((x) > 0 ? (x) : -(x))
#define PI 3.14159265358979323846

typedef struct prm {
  int NX;
  int NY;
  int nx;
  int ny;
  double LX;
  double LY;
  double dt;
  double dx;
  double dy;
  double U;
  double L;
  double nu;
  double Re;
} Prm;

// void Semilag_slow(double* u, double* v, double* q, Prm prm, int sign) {
//   double x, y, x0, y0;
//   int i0, j0;
//   double* ADVq = calloc(prm.NX * prm.NY, sizeof(double));
//   double* Mx2 = calloc(prm.NX * prm.NY, sizeof(double));
//   double* My2 = calloc(prm.NX * prm.NY, sizeof(double));
//   double* Mx1 = calloc(prm.NX * prm.NY, sizeof(double));
//   double* My1 = calloc(prm.NX * prm.NY, sizeof(double));
//   double u_abs, v_abs;
//   double Cc, Ce, Cmx, Cmy;
//   for (int i = 1; i < prm.NX - 1; i++) {
//     for (int j = 1; j < prm.NY - 1; j++) {
//       Mx2[i * prm.NY + j] = SIGN(SIGN(sign * U(i, j)) + 1);
//       Mx1[i * prm.NY + j] = 1. - Mx2[i * prm.NY + j];
//       My2[i * prm.NY + j] = SIGN(SIGN(sign * V(i, j)) + 1);
//       My1[i * prm.NY + j] = 1. - My2[i * prm.NY + j];

//       u_abs = ABS(U(i, j));
//       v_abs = ABS(V(i, j));

//       Cc = (prm.dx - u_abs * prm.dt) * (prm.dy - v_abs * prm.dt) / (prm.dx * prm.dy);
//       Ce = prm.dt * prm.dt * u_abs * v_abs / (prm.dx * prm.dy);
//       Cmx = (prm.dx - u_abs * prm.dt) * (v_abs * prm.dt) / (prm.dx * prm.dy);
//       Cmy = (u_abs * prm.dt) * (prm.dy - v_abs * prm.dt) / (prm.dx * prm.dy);

//       ADVq[i * prm.NY + j] = Cc * q[i * prm.NY + j] +
//                              Ce * (Mx1[i * prm.NY + j] * My1[i * prm.NY + j] * q[(i + 1) * prm.NY + j + 1] +
//                                    Mx2[i * prm.NY + j] * My1[i * prm.NY + j] * q[(i - 1) * prm.NY + j + 1] +
//                                    Mx1[i * prm.NY + j] * My2[i * prm.NY + j] * q[(i + 1) * prm.NY + j - 1] +
//                                    Mx2[i * prm.NY + j] * My2[i * prm.NY + j] * q[(i - 1) * prm.NY + j - 1]) +
//                              Cmx * (My2[i * prm.NY + j] * q[i * prm.NY + j - 1] +
//                                     My1[i * prm.NY + j] * q[i * prm.NY + j + 1]) +
//                              Cmy * (Mx2[i * prm.NY + j] * q[(i - 1) * prm.NY + j] +
//                                     Mx1[i * prm.NY + j] * q[(i + 1) * prm.NY + j]);
//     }
//   }
//   memcpy(q, ADVq, prm.NX * prm.NY * sizeof(double));
//   free(ADVq);

//   // for (int i = 1; i < prm.NX - 1; i++) {
//   //   for (int j = 1; j < prm.NY - 1; j++) {
//   //     x = -0.5 * prm.dx + i * prm.dx;
//   //     y = -0.5 * prm.dy + j * prm.dy;
//   //     if (U(i, j) > 0)
//   //       x0 = x - U(i, j) * prm.dt;
//   //     else
//   //       x0 = x + U(i, j) * prm.dt;
//   //     if (V(i, j) > 0)
//   //       y0 = y - V(i, j) * prm.dt;
//   //     else
//   //       y0 = y + V(i, j) * prm.dt;

//   //     i0 = (int)(x0 / prm.dx + 0.49999999);
//   //     j0 = (int)(y0 / prm.dy + 0.49999999);
//   //     // printf("x = %lf, x0 = %lf, y = %lf, y0 = %lf, i0 = %d, j0 = %d\n", x, x0, y, y0, i0, j0);

//   //     double C00, C01, C10, C11;
//   //     C00 = (x(i0 + 1) - x0) * (y(j0 + 1) - y0);
//   //     C10 = (x0 - x(i0)) * (y(j0 + 1) - y0);
//   //     C01 = (x(i0 + 1) - x0) * (y0 - y(j0));
//   //     C11 = (x0 - x(i0)) * (y0 - y(j0));

//   //     C00 /= prm.dx * prm.dy;
//   //     C10 /= prm.dx * prm.dy;
//   //     C01 /= prm.dx * prm.dy;
//   //     C11 /= prm.dx * prm.dy;

//   //     ADVq[i * prm.NY + j] =
//   //         C00 * q[i0 * prm.NY + j0] +
//   //         C10 * q[(i0 + 1) * prm.NY + j0] +
//   //         C01 * q[i0 * prm.NY + j0 + 1] +
//   //         C11 * q[(i0 + 1) * prm.NY + j0 + 1];
//   //   }
//   // }
//   // memcpy(q, ADVq, prm.NX * prm.NY * sizeof(double));
//   // free(ADVq);
// }

// 1st order semi-Lagrangian advection
void Semilag(double* u, double* v, double* q, Prm prm, int sign) {
  int sign_u, sign_v;
  double a, b;
  double* aux = (double*)malloc(prm.NX * prm.NY * sizeof(double));
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
  memcpy(q, aux, prm.NX * prm.NY * sizeof(double));
  free(aux);
}

// 2nd order semi-Lagrangian advection
void Semilag2(double* u, double* v, double* q, Prm prm) {
  double* q0 = (double*)malloc(prm.NX * prm.NY * sizeof(double));
  memcpy(q0, q, prm.NX * prm.NY * sizeof(double));
  Semilag(u, v, q, prm, 1);
  Semilag(u, v, q, prm, -1);
  for (int i = 0; i < prm.NX * prm.NY; i++)
    q[i] = q0[i] + (q0[i] - q[i]) / 2;
  Semilag(u, v, q, prm, 1);
}

void init(Prm* prm) {
  prm->NX = 100;                        // total grid points (including ghost points)
  prm->NY = 100;                        // total grid points (including ghost points)
  prm->nx = prm->NX - 2;                // real grid points
  prm->ny = prm->NY - 2;                // real grid points
  prm->LX = 10.0;                       // Length of domain (x direction)
  prm->LY = 5.0;                        // Length of domain (y direction)
  prm->dx = prm->LX / prm->nx;          // grid spacing (x direction)
  prm->dy = prm->LY / prm->ny;          // grid spacing (y direction)
  prm->dt = 0.001;                      // time step
  prm->U = 1.0;                         // enterning speed
  prm->L = 1.0;                         // characteristic length of the object
  prm->nu = 1.51e-5;                    // kinematic viscosity of air at atmospheric conditions
  prm->Re = prm->U * prm->L / prm->nu;  // Reynolds number
}

// we set the boundary conditions using the ghost points and linear interpolation
void boundary_conditions(double* u, double* v, double* p, Prm prm) {
  for (int j = 1; j < prm.NY - 1; j++) {
    // left boundary: u = U, v = 0, \partial_x p = 0
    U(0, j) = 2 * prm.U - U(1, j);
    V(0, j) = -V(1, j);
    P(0, j) = P(1, j);
    // right boundary: \partial_x u = 0, \partial_x v = 0, p = 0
    U(prm.NX - 1, j) = U(prm.NX - 2, j);
    V(prm.NX - 1, j) = V(prm.NX - 2, j);
    P(prm.NX - 1, j) = -P(prm.NX - 2, j);
  }
  for (int i = 1; i < prm.NX - 1; i++) {
    // bottom boundary: \partial_y u = 0, v = 0, \partial_y p = 0
    U(i, 0) = U(i, 1);
    V(i, 0) = -V(i, 1);
    P(i, 0) = P(i, 1);
    // top boundary: \partial_y u = 0, v = 0, \partial_y p = 0
    U(i, prm.NY - 1) = U(i, prm.NY - 2);
    V(i, prm.NY - 1) = -V(i, prm.NY - 2);
    P(i, prm.NY - 1) = P(i, prm.NY - 2);
  }
}

int main(void) {
  Prm prm;
  init(&prm);

  double t = 0.0;

  // allocate memory
  double* u = (double*)calloc(prm.NX * prm.NY, sizeof(double));  // x component of velocity, initialized to 0 by calloc
  double* v = (double*)calloc(prm.NX * prm.NY, sizeof(double));  // y component of velocity, initialized to 0 by calloc
  double* p = (double*)calloc(prm.NX * prm.NY, sizeof(double));  // pressure, initialized to 0 by calloc

  // // set initial condition
  // for (int i = 1; i < prm.NX - 1; i++) {
  //   for (int j = 1; j < prm.NY - 1; j++) {
  //     U(i, j) = 0;
  //     V(i, j) = -sin(2 * PI * x(i)) * sin(PI * y(j)) * sin(PI * y(j));
  //     r = sqrt((x(i) - X0) * (x(i) - X0) + (y(j) - Y0) * (y(j) - Y0));
  //     if (r > alpha) r = alpha;
  //     r /= alpha;
  //     q[i * prm.NY + j] = 0.25 * (1 + cos(PI * r));
  //   }
  // }

  // count time
  clock_t start, end;
  start = clock();
  while (t < 2) {
    // advection term
    Semilag2(u, v, u, prm);
    Semilag2(u, v, v, prm);
    t += prm.dt;
  }
  end = clock();
  printf("time = %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

  free(u);
  free(v);

  return 0;
}
