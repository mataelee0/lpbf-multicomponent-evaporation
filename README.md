# lpbf-multicomponent-evaporation
A multicomponent evaporation–thermo-fluid model incorporating non-ideal thermodynamic activities and composition feedback for laser powder bed fusion.

## Requirements
1. OpenFOAM Foundation version 10 is required.
   The OpenFOAM-10 source code is available at:
   https://openfoam.org/download/10/
2. The laserbeamFoam solver developed by Flint et al. must be installed.

## Solver compilation
The solver directory contains source code for compiling the
laserbeamFoam-based solver.

Before compilation, the include path to the `laserHeatSource` library
installed with laserbeamFoam must be specified in the solver build options.

Specifically, the path to `laserHeatSource/lnInclude` should be added
to the `EXE_INC` entry in the `Make/options` file.
The correct path can be identified from the local laserbeamFoam installation.

1. Navigate to the solver directory:
```bash
cd solver
```
2. Compile the solver:
```bash
./Allwmake
```

After successful compilation, the solver can be executed using:
```bash
   laserbeamFoam
```

## Case setup and execution
Simulation cases are provided in the following directories:
- run_stationary
- run_scanning

### Running a simulation
From within either the run_stationary or run_scanning directory,
execute the following commands:
1. Generate the computational mesh:
```bash
blockMesh
```
2. Decompose the domain for parallel execution:

The initial fields, including alpha.metal, are already initialized in the 0/ directory. Therefore, no additional field initialization is required.
```bash
decomposePar
```

### Solver execution
Serial execution:
   laserbeamFoam
Parallel execution (example using 32 processors):
```bash
mpirun -np 32 laserbeamFoam -parallel >& log
```
