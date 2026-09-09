# Run AMReX flow solver with Chombo SDF geometry

## Installation
```
git clone https://github.com/nataraj2/amrex.git
git clone https://github.com/amrex-suite/chombo_to_amrex.git
cd amrex
git checkout chombo_SDF_to_amrex_EB2 
cd .. 
```

In the inputs file give the path to the chombo hdf5 file.
```
eb2.chombo_sdf_hdf5_file=<path-to=chombo-hdf5-file>
```

# Compilation
```
cd chombo_to_amrex/CNS/Exec/Chombo_AMReX
make -j8
```


## Run
```
./CNS3d.gnu.x86-milan.MPI.ex inputs_simple_plane
```
