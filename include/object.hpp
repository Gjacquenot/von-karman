#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <cmath>
#include <iostream>

#include "misc.hpp"

using namespace std;

typedef struct point {
  double x;
  double y;
} Point;

typedef struct grid {
  int i;
  int j;
} Grid;

#define TOL 1e-8

class Object {
 public:
  bool* IsInside;       // 1 if the point is inside the object, 0 otherwise.
  bool* IsGhost;        // 1 if the point is a ghost point, 0 otherwise.
  Grid* GhostPoints;    // Array of points that are just inside the object.
  Point* MirrorPoints;  // Array of points that are mirrored by the object from the ghost points.
  // create a matrix of 6 interpolating points for each ghost point
  Grid* InterpolatingPoints;  // Array of points that are used to interpolate any quantity (velocity u, velocity v or pressure) at the mirror points.
  int count_ghost;

  void init(Prm prm) {
    IsInside = new bool[prm.NX * prm.NY];
    IsGhost = new bool[prm.NX * prm.NY];
    int count_gh = 0;
    for (int i = 0; i < prm.NX; i++) {
      for (int j = 0; j < prm.NY; j++) {
        // set the IsInside array
        IsInside[i * prm.NY + j] = is_inside(i, j, prm);

        // set the IsGhost array
        if (IsInside[i * prm.NY + j] && (!is_inside(i - 1, j, prm) || !is_inside(i + 1, j, prm) || !is_inside(i, j - 1, prm) || !is_inside(i, j + 1, prm))) {
          IsGhost[i * prm.NY + j] = true;
          count_gh++;
        } else {
          IsGhost[i * prm.NY + j] = false;
        }
      }
    }
    this->count_ghost = count_gh;
    GhostPoints = new Grid[count_gh];
    MirrorPoints = new Point[count_gh];
    InterpolatingPoints = new Grid[count_gh * 6];
    int count = 0;
    for (int i = 0; i < prm.NX; i++) {
      for (int j = 0; j < prm.NY && count < count_gh; j++) {
        if (IsGhost[i * prm.NY + j]) {
          GhostPoints[count] = {i, j};
          MirrorPoints[count] = mirror_point(x(i), y(j));
          set_interpolating_points(i, j, MirrorPoints[count], count, prm);
          count++;
        }
      }
    }
  }

  // checks if a point is inside the boundary
  virtual bool is_inside(int i, int j, Prm prm) = 0;

  virtual Point closest_boundary_point(double x, double y) = 0;

  Point mirror_point(double x, double y) {
    Point p = closest_boundary_point(x, y);
    return {2 * p.x - x, 2 * p.y - y};
  }

  int sign_nx(double x, double mirror_x) {
    if (mirror_x - x > TOL) {
      return 1;
    } else if (mirror_x - x < -TOL) {
      return -1;
    } else {
      return 0;
    }
  }

  int sign_ny(double y, double mirror_y) {
    if (mirror_y - y > TOL) {
      return 1;
    } else if (mirror_y - y < -TOL) {
      return -1;
    } else {
      return 0;
    }
  }

  void set_interpolating_points(int i, int j, Point mirror, int count, Prm prm) {
    // based on algorithm of the paper: High order ghost-cell immersed boundary method for generalized boundary conditions (by Mehrdad Yousefzadeh, Ilenia Battiato), p. 589-590

    // find I, J such that x(I) < mirror.x < x(I+1) and y(J) < mirror.y < y(J+1)
    int I = i - 2, J = j - 2;  // minimum value for I and J (the maximum value is I+2 and J+2)
    while (x(I + 1) < mirror.x + TOL) I++;
    while (y(J + 1) < mirror.y + TOL) J++;

    int x1, x2, x3, x4, x5, x6_1, x6_2, x6;
    int y1, y2, y3, y4, y5, y6_1, y6_2, y6;
    double d_1, d_2;
    int sgn_nx = sign_nx(x(i), mirror.x);
    int sgn_ny = sign_ny(y(j), mirror.y);
    if (sgn_nx > 0) {
      x1 = I + 1;
      x2 = I + 2;
      x5 = I + 1;
      x6_1 = I + 1;
      x6_2 = I + 2;
    } else {  // sgn_nx <= 0
      I++;
      x1 = I - 1;
      x2 = I - 2;
      x5 = I - 1;
      x6_1 = I - 1;
      x6_2 = I - 2;
    }
    x3 = I;
    x4 = I;

    if (sgn_ny > 0) {
      y3 = J + 1;
      y4 = J + 2;
      y5 = J + 1;
      y6_1 = J + 2;
      y6_2 = J + 1;
    } else {
      J++;
      y3 = J - 1;
      y4 = J - 2;
      y5 = J - 1;
      y6_1 = J - 2;
      y6_2 = J - 1;
    }
    y1 = J;
    y2 = J;

    // choose from x6_1 and x6_2
    d_1 = (x(x6_1) - mirror.x) * (x(x6_1) - mirror.x) + (y(y6_1) - mirror.y) * (y(y6_1) - mirror.y);
    d_2 = (x(x6_2) - mirror.x) * (x(x6_2) - mirror.x) + (y(y6_2) - mirror.y) * (y(y6_2) - mirror.y);

    if (d_1 < d_2) {
      x6 = x6_1;
      y6 = y6_1;
    } else {
      x6 = x6_2;
      y6 = y6_2;
    }

    InterpolatingPoints[count * 6] = {x1, y1};
    InterpolatingPoints[count * 6 + 1] = {x2, y2};
    InterpolatingPoints[count * 6 + 2] = {x3, y3};
    InterpolatingPoints[count * 6 + 3] = {x4, y4};
    InterpolatingPoints[count * 6 + 4] = {x5, y5};
    InterpolatingPoints[count * 6 + 5] = {x6, y6};
  }
};

class Circle : public Object {
 public:
  double x0, y0, r;
  Circle(double x0_, double y0_, double r_, Prm prm) : x0(x0_), y0(y0_), r(r_) {
    init(prm);
  }
  bool is_inside(int i, int j, Prm prm) override {
    return (x(i) - x0) * (x(i) - x0) + (y(j) - y0) * (y(j) - y0) < r * r;
  }

  Point closest_boundary_point(double x, double y) override {
    double d = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
    return {x0 + r * (x - x0) / d, y0 + r * (y - y0) / d};
  }
};

class Rectangle : public Object {
 public:
  double x0, y0, Lx, Ly;
  Rectangle(double x0_, double y0_, double Lx_, double Ly_, Prm prm) : x0(x0_), y0(y0_), Lx(Lx_), Ly(Ly_) {
    init(prm);
  }

  bool is_inside(int i, int j, Prm prm) override {
    return (x(i) - x0) * (x(i) - x0) < Lx * Lx / 4 && (y(j) - y0) * (y(j) - y0) < Ly * Ly / 4;
  }

  Point closest_boundary_point(double x, double y) override {
    double dx1 = x - (x0 - Lx / 2);
    double dx2 = (x0 + Lx / 2) - x;
    double dx = (dx1 < dx2) ? dx1 : dx2;
    double dy1 = y - (y0 - Ly / 2);
    double dy2 = (y0 + Ly / 2) - y;
    double dy = (dy1 < dy2) ? dy1 : dy2;
    if (dx < dy) {
      return {(x > x0) ? x0 + Lx / 2 : x0 - Lx / 2, y};
    } else {
      return {x, (y > y0) ? y0 + Ly / 2 : y0 - Ly / 2};
    }
  }
};

#endif  // OBJECT_HPP