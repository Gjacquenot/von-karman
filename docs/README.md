# Numerical simulation of the Von Kármán vortex street

This project is dedicated to the numerical simulation of the Kármán vortex street in two dimensions, utilizing finite difference, the Chorin method and a Semi-Lagrangian solver. Users can explore various boundary conditions, initial conditions, and visualize the evolution and final states of simulations. 

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

If you get an error, you may need to activate the execution permissions of the file `run.sh` with the following command:

```
chmod +x run.sh
```

If activated in the file `config/input.txt` (check the readme of the input file in `docs/input_README.txt`), the code will also animate the solution.

## Other information

The code is mainly in the `src` folder, and the headers are in the `include` folder. In the `config` folder you can find the input file and the output files. Finally, in the `docs` folder you can find the documentation of the input file and other papers with information about the Von Kármán vortex street and numerical methods used in the code.

### Instructions to add a new object shape

In order to add a new object shape, you will need to create a new class in the `object.hpp` file inside the `include` folder. This class will need to inherit from the `Object` class and implement the `is_inside`and `closest_boundary_point` methods (you can check the other classes in the file to see the general idea in how to implement them). Then, in the `src/main.cpp` you will need to add an equivalent `else if` condition for the detection of the object. Finally, you will need to add the object in the python script `src/plot_functions.py` in order to plot the object in the animation.

Finally, create a new input file in the `config` folder with the parameters of the new object and run the code (remember to change the obstacle section in `config/input.txt` with the name of your new shape).

## Results

The following is an example of the evolution of the Von Kármán vortex street, with different objects and Reynolds numbers. We plot the vorticity field in order to visualize more clearly the formation of the vortices.

<!-- add video in data/videos/circle_Re=500.0.mp4 -->

https://github.com/victorballester7/von-karman/assets/78110444/ef5da9ef-3f67-44c4-bcba-cd4eb07d0270

Circle at Re = 500

https://github.com/victorballester7/von-karman/assets/78110444/2d75f60b-4d38-431c-b721-f691c689a034

Circle with a fin at Re = 500

https://github.com/victorballester7/von-karman/assets/78110444/77fc098b-49ff-458b-99cc-7e5267056435

Airfoil at Re = 5000
