# Numerical simulation of the Von Kármán vortex street

## Compilation of the code

The code is written in C++. You will need the library `Eigen` to compile it, as it is used for the linear algebra operations (with sparse matrices). You can check its webpage from [here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

To compile the code, first download it:

```
git clone git@github.com:victorballester7/von-karman.git
cd von-karman
```

Then, you can compile and run the code (having previously check your input parameters in the file `data/input.txt`) with the following command:

```
./run.sh
```

If activated in the file `data/input.txt` (check the readme of the input file in `docs/input_README.txt`), the code will also animate the solution.
