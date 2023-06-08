# Parameter Estimation in Dynamical Systems

This repository contains code to estimate the parameters of a dynamical system
using gradient descent. The `config/` directory contains a set of configuration
files that are used by the program files located in the `code/` directory. There
is also a report located in the `report/` directory that provides an overview of
the problem and a summary of the results.

To run the code, first create a symbolic link in the `config/` directory for 
either macOS or Linux. If you are on macOS there are configuration files for 
Homebrew and MacPorts. Choose the file that aligns with your package manager. If
using macOS, you may also have to modify the compiler to align with your
system's configuration. For example, the Homebrew file sets `g++-12` as the
compiler, but you may have to change this to `g++-10` or `g++-11`, depending on
your own system's preferences.

Navigate to the `code/` directory, and type `make` in the terminal to compile
the code. The code is structured as follows:

* `dopri54`: a general-purpose ODE solver that implements the Dormand-Prince 
  5(4) method
* `lotka_volterra`: a class that specifically implements the right-hand side of 
  the Lotka-Volterra equations
* `lv_dop54`: a subclass of `dopri54` and `lotka_volterra` that solves the
  Lotka-Volterra equations using the Dormand-Prince algorithm
* `lv_sensitivity`: a class that specifically implements the right-hand side of
  the Lotka-Volterra equations with the sensitivity equations included
* `lvs_dop54`: a subclass of `dopri54` and `lv_sensitivity` that solves the
  Lotka-Volterra equations with sensitivity using the Dormand-Prince algorithm
* `param_est.cc`: front-end script that performs the parameter estimation

