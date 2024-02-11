#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>

#define U(i, j) u[(i) * prm.NY + (j)]
#define V(i, j) v[(i) * prm.NY + (j)]
#define P(i, j) p[(i) * prm.NY + (j)]
#define x(i) (-0.5 * prm.dx + (i) * prm.dx)
#define y(j) (-0.5 * prm.dy + (j) * prm.dy)

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
  size_t size_ghost;
  size_t size;
} Prm;

#endif  // MISC_HPP