# Compressible Navier–Stokes Solver with Chombo Signed Distance Function (SDF) Converted to AMReX EB2

This folder contains the code to run a compressible Navier–Stokes solver with embedded geometry created from the Chombo signed distance function (SDF).

## Build and Run

```bash
cd Exec/Pulse
make -j8
mpirun -np <num_ranks> <exec> inputs
```

This will read the plotfile written by the `ReadSDFandWritePlotfile` code, create an AMReX EB from the SDF, and run a simulation.
