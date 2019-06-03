# NumericPDE

A repository containing practical implementations of the finite-difference method to partial differential equations. 
This repository contains some projects from a module I completed at York, as well as some of my own work. The code written in this repository should be able to run in both MATLAB and Octave.
For the projects from my numerical methods for PDEs module at York, there are pdf documents accompanying the code to explain the problems I was trying to solve.

### Project 1
Parabolic Heat Equation implemented solved using the Crank-Nicolson method. Implemented using the double-sweep method.

- PDE problem and implementation [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/pde_project_1.pdf)
- Files required
   * Crank-Nicolson implementation [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/project_c_n.m) 
   * Plotting results [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/plot_u_on_single_figure_1.m)

### Project 2
PDE problem and implementation of code [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/pde_project_2.pdf)
1. Solve a generalised PDE using MATLAB's built-in `pdepe` function.
   * Specifying BC [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/pde_1_bc.m)
   * Specifying IC [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/pde_1_ic.m)
   * Specifying PDE [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/pde_1_function.m)
   * Plotting results [here](https://github.com/thomasarmstrong98/numericPDE/blob/master/plot_u_on_single_figure.m)
2. Solving a two-dimensional heat equation with the ADI method.
   * ADI implementation [here](ADI_method.m)
