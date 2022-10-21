# A Simplified Atrial Myocyte Model Electrophysiologically Coupled to a SAN Pacemaker Cell Model

- Atrial cell: A simplified atrial cell model (Euler based, 2022) based upon the H. Zhang et al. (2020) atrial cell model.
- SAN pacemaker cell: Restructured C++ version of Wei's SAN model (2012) based upon the S. Kharche (2011) SAN model.

## Configuring Simulations

1. Alter variables within atrial_cell_configuration.hpp and SAN_cell_configuration.hpp to desired cell properties.

2. Configure experiment as a whole (coupled ODE solver timesteps, output timeteps, temperatures...) within
   simulation_configuration.hpp. Both "output_dt" and "dt" are driving factors in processing times.

## Running Sumlations

1. Ensure initial cell-state variable files are present for both the atrial and SAN cells within the input folder.

2. Configure simulation as above. Name the current simulation, if desired.

3. Build an executable binary and run the program.

4. Program will output data, final states, measurements and the applied configurations to the output folder.
   Some measurements will also be outputted to the terminal, which informs the user of the current beat as
   the simulation progresses.

## Other notes

- Cell parameters are defined within each cell's class declaration.

- Follow the function's called in main.cpp for an overview of the program's structure. Both the SAN and atrial
  cells are derived from an abstract base class "base_cell". All respective files for each cell contain the cell
  type name as a prefix. The "Simulation" class contains members that apply to the simulation as a whole.

- No measurements are taken in the first beat.

- Program compiled and tested on Windows and Linux OS, with GNU's g++ compiler.

Matt Hayter, University of Manchester
