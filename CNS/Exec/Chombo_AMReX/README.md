# Run AMReX flow solver with Chombo SDF geometry

## Installation
```
git clone --recursive https://github.com/amrex-suite/chombo_to_amrex.git
```

## Compilation
```
cd chombo_to_amrex/CNS/Exec/Chombo_AMReX
make -j8
```

## Run
In the inputs file give the path to the chombo hdf5 file.
```
eb2.chombo_sdf_hdf5_file=<path-to=chombo-hdf5-file>
```
Run the flow solver. Currently, only serial run is supported.
```
./CNS3d.gnu.x86-milan.MPI.ex inputs_simple_plane
```
