Lx                      Length of the domain in the x-direction. Positive real number.
Ly                      Length of the domain in the y-direction. Positive real number.
nx                      Number of grid points in the x-direction (excluding ghost cells). Positive integer.
ny                      Number of grid points in the y-direction (excluding ghost cells). Positive integer.
dt                      Initial time step (if all goes well, it should be adaptative). Positive real number.
Tfinal                  Time at which the simulation should stop. Positive real number.
Re                      Reynolds number. Positive real number.
nu                      Kinematic viscosity of the fluid. Positive real number.
vorticity               Whether to plot the vorticity or not. Boolean: 1 (yes) or 0 (no).
obstacle                Whether to include an obstacle in the domain or not. Boolean: 1 (yes) or 0 (no).
object                  Type of obstacle. Default strings: "circle", "circle_fin", "rectangle" or "airfoil".
plot_dt                 Time step for the animation (print solution each animation_dt time-steps). Positive integer.
animation               Whether to animate the solution or not. Boolean: 1 (yes) or 0 (no).
