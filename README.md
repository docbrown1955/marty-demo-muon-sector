# MARTY Demo - Muon Sector

Muon sector example (self-energy + anomalous magnetic moment) using MARTY in a simple QED model.

To download and install MARTY follow the instructions [here]().

## Summary

The example presents the complete calculation, including the model building, of:
 - The muon self-energy (divergent contribution to an off-shell one-loop amplitude)
 - The muon anomalous magnetic moment (g-2)µ at the one-loop level

## Execute the MARTY program

Just type

  make
  ./main

In case the compilation does not work, just change the compiler in the Makefile to set it to your C++ compiler.

The model will be displayed and several results of calculations. The program will ask for input step-by-step just to pause the program, and GRAFED will be launched displaying the relevant Feynman diagrams for the unique vertex in the theory and the two calculations (self-energy and magnetic moment).
Execute the numerical example the check the numbers

Just type

cp example_demolib.cpp demolib/script
cd demolib
make
bin/example_demolib.x

This will execute the program and evaluate the quantities that have been calculated to compare them to their theoretical values.
