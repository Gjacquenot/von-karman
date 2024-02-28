#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>

#define U(i, j) u[(i) * prm.NY + (j)]
#define V(i, j) v[(i) * prm.NY + (j)]
#define W(i, j) w[(i) * prm.NY + (j)]
#define ADV_U(i, j) adv_u[(i) * prm.NY + (j)]
#define ADV_V(i, j) adv_v[(i) * prm.NY + (j)]
#define Ustar(i, j) ustar[(i) * prm.NY + (j)]
#define Vstar(i, j) vstar[(i) * prm.NY + (j)]
#define P(i, j) p[(i) * prm.NY + (j)]
#define DIV(i, j) div((i) * prm.ny + (j))
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
  size_t NXNY;
  size_t nxny;
  bool obstacle_ON;
} Prm;

#endif  // MISC_HPP