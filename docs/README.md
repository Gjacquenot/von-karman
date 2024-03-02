# Numerical simulation of the Von Kármán vortex street

This project is dedicated to the numerical simulation of the Kármán vortex street in two dimensions, utilizing finite difference, the Chorin method and a Semi-Lagrangian solver. Users can explore various boundary conditions, initial conditions, and visualize the evolution and final states of simulations. The

## Description

The Kármán vortex street is a repeating pattern of swirling vortices caused by the unsteady separation of flow of a fluid around blunt bodies, such as a cylinder or a sphere, at certain ranges of Reynolds numbers. This phenomenon is characterized by its strikingly regular formation of alternating vortices downstream of the obstacle, which can lead to oscillating forces on the object and is studied extensively in fluid dynamics for its implications in engineering and natural systems.

## Compilation and execution of the code

The code is written in `C++` and the animation is done in `Python`. You will need the library `Eigen` to compile it, as it is used for the linear algebra operations (with sparse matrices). You can check its webpage from [here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

In order to properly compile the code, you will to ensure that the `eigen3` folder is correctly linked in the `Makefile` file. Either you have the `eigen3` folder in the path `/usr/include/`, being this folder in the path of the system, or you will need to specify the library path in the `Makefile` in the `INCLUDES` variable. Keep in mind that you may also need to change `python` to `python3` in the `run.sh` file if you are using `python3` instead of `python`.

To compile the code, first download it:

```
git clone git@github.com:victorballester7/von-karman.git
cd von-karman
```

Then, you can compile and run the code (having previously check your input parameters in the file `config/input.txt`) with the following command:

```
./run.sh
```

If activated in the file `config/input.txt` (check the readme of the input file in `docs/input_README.txt`), the code will also animate the solution.

## Other information

The code is mainly in the `src` folder, and the headers are in the `include` folder. In the `data` folder you can find the input file and the output files. Finally, in the `docs` folder you can find the documentation of the input file and other papers with information about the Von Kármán vortex street and numerical methods used in the code.
