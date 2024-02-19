#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define U(i, j) u[(i) * prm.NY + (j)]
#define V(i, j) v[(i) * prm.NY + (j)]
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
  double Re;
} Prm;

// 1st order semi-Lagrangian advection
void Semilag(double* u, double* v, double* q, double* aux, Prm prm, int sign) {
  int sign_u, sign_v;
  double a, b;
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
}

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

// 2nd order semi-Lagrangian advection
void Semilag2(double* u, double* v, double* q, double* q0_aux, double* aux, Prm prm) {
  memcpy(q0_aux, q, prm.NX * prm.NY * sizeof(double));
  Semilag(u, v, q, aux, prm, 1);
  Semilag(u, v, q, aux, prm, -1);
  for (int i = 0; i < prm.NX * prm.NY; i++)
    q[i] = q0_aux[i] + (q0_aux[i] - q[i]) / 2;
  Semilag(u, v, q, aux, prm, 1);
}

int try() {
  return 2;
}

int main(void) {
  Prm prm;
  prm.NX = 102;
  prm.NY = 100;
  prm.nx = prm.NX - 2;
  prm.ny = prm.NY - 2;
  prm.LX = 1.0;
  double aspect_ratio = 1.0;
  prm.LY = prm.LX / aspect_ratio;
  prm.dx = prm.LX / prm.nx;
  prm.dy = prm.LY / prm.ny;
  prm.dt = 0.001;
  prm.Re = 100.0;
  double t = 0.0;
  // printf("dx = %lf, dy = %lf\n", prm.dx, prm.dy);
  // Allocate arrays
  double* u = (double*)calloc(prm.NX * prm.NY, sizeof(double));
  double* v = (double*)calloc(prm.NX * prm.NY, sizeof(double));
  double* aux = (double*)calloc(prm.NX * prm.NY, sizeof(double));
  double* q = (double*)calloc(prm.NX * prm.NY, sizeof(double));
  double* q_mine = (double*)calloc(prm.NX * prm.NY, sizeof(double));
  double* q0 = malloc(prm.NX * prm.NY * sizeof(double));
  double* ADVq = (double*)calloc(prm.NX * prm.NY, sizeof(double));
  double r, alpha = 0.2, X0 = 0.5, Y0 = 0.75;

  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      U(i, j) = sin(PI * x(i)) * sin(PI * x(i)) * sin(2 * PI * y(j));
      V(i, j) = -sin(2 * PI * x(i)) * sin(PI * y(j)) * sin(PI * y(j));
      r = sqrt((x(i) - X0) * (x(i) - X0) + (y(j) - Y0) * (y(j) - Y0));
      if (r > alpha) r = alpha;
      r /= alpha;
      q[i * prm.NY + j] = 0.25 * (1 + cos(PI * r));
    }
  }

  memcpy(q_mine, q, prm.NX * prm.NY * sizeof(double));

  // count time
  clock_t start, end;
  start = clock();
  while (t < 2) {
    Semilag2(u, v, q, q0, aux, prm);
    t += prm.dt;
  }
  end = clock();
  printf("time = %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

  double max_diff = 0.0, s;
  printf("t = %lf\n", t);
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      s = fabs(q[i * prm.NY + j] - q_mine[i * prm.NY + j]);
      if (max_diff < s) max_diff = s;
    }
  }
  printf("max_diff = %g\n", max_diff);

  free(u);
  free(v);
  free(q);
  free(ADVq);

  return 0;
}
