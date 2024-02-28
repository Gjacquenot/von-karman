#ifndef NS_HPP
#define NS_HPP

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "../include/object.hpp"
#include "misc.hpp"
#define p(jj) (obstacle.InterpolatingPoints[6 * i + (jj)])
#define m_x (obstacle.MirrorPoints[i].x)
#define m_y (obstacle.MirrorPoints[i].y)

typedef Eigen::SparseMatrix<double> SpMat;  // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;        // declares a triplet type of double

using namespace std;

// @brief 1st order semi-Lagrangian for advection
// @param u velocity in x direction
// @param v velocity in y direction
// @param q scalar field to be advected
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param sign sign of the direction of the velocity fields (to be easily inverted). If unsure, use 1.
// @param obstacle object representing the obstacle in the domain
void Semilag(double* u, double* v, double* q, Prm prm, int sign, Object& obstacle);

// @brief 2nd order semi-Lagrangian for advection
// @param u velocity in x direction
// @param v velocity in y direction
// @param q0 initial scalar field to be advected
// @param q1 final scalar field to be advected
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void Semilag2(double* u, double* v, double* q0, double* q1, Prm prm, Object& obstacle);

void upwind(double* u, double* v, double* q, Prm prm, int sign, Object& obstacle);
void upwind2(double* u, double* v, double* q0, double* q1, Prm prm, Object& obstacle);

// @brief Set the boundary conditions for the velocity fields, but on the ghost points (i.e. fictitious points outside the domain) by linear interpolation
// @param u velocity in x direction
// @param v velocity in y direction
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void BC_velocity(double* u, double* v, Prm prm, Object& obstacle);
// @brief Set the boundary conditions for the pressure field, but on the ghost points (i.e. fictitious points outside the domain) by linear interpolation
// @param p pressure
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void BC_pressure(double* p, Prm prm, Object& obstacle);
double interpolate(double x, double y, double* phi, Prm prm);
void set_vorticity(double* u, double* v, double* w, Prm prm);
// @brief Build the Laplace matrix for the pressure Poisson equation laplace(p) = div(u)
// @param coeffs vector of triplets to store the coefficients of the Laplace matrix. Triplets consist of (i, j, value), where i and j are the indices of the matrix and value is the value of the coefficient
// @param prm parameters of the simulation (dx, dy, dt, etc.)
void buildLaplaceMatrix(vector<Trip>& coeffs, Prm prm);

#endif  // NS_HPP