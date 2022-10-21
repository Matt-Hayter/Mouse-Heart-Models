Matt Hayter, University of Manchester
05/2022

A simplification of the Zhang et al. (2020) atrial myocyte model, for steady stimulus pacing.

*****Running Sumlations*****

1. Prepare the initial cell-state variable file for the atrial myocyte within the input folder. For measurements, 
   these should typically be stabalised output states from a previous simulation (under the same experimental conditions).

2. Configure the atrial cell and simulation within the "configuration.hpp" file.

3. Build an executable binary and run the program.

4. Program will output data, final state-variables, measurements and the applied configuration to the output_folder.
   Some measurements will also be outputted to the terminal, including the simulation's current beat number.

*****Key Alterations in Simplified Model*****

- Numerical integration method converted from using BDFs and Newton methods to the forward Euler method (validated).

- Computational implementation restructured and simplified, in C++ only. No dependencies on external libraries/external
  ODE solvers.

- Separated from the generic cell model, from which the Zhang et al. atrial cell model was derived. Configuration
  file designed to be self-contained.

*****Other notes*****

- Cell parameters are defined within the cell's class declaration (or configuration file).

- No measurements are taken in the first beat.

- Program compiled and tested on Windows and Linux OS, with GNU's g++ compiler.