#include <AMReX.H>
#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include <hdf5.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include <hdf5.h>
struct H5Box
{
    int lo_i, lo_j, lo_k;
    int hi_i, hi_j, hi_k;
};

struct H5IntVect
{
    int i, j, k;
};

static hid_t make_box_type()
{
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(H5Box));

    H5Tinsert(type, "lo_i", HOFFSET(H5Box,lo_i), H5T_NATIVE_INT);
    H5Tinsert(type, "lo_j", HOFFSET(H5Box,lo_j), H5T_NATIVE_INT);
    H5Tinsert(type, "lo_k", HOFFSET(H5Box,lo_k), H5T_NATIVE_INT);

    H5Tinsert(type, "hi_i", HOFFSET(H5Box,hi_i), H5T_NATIVE_INT);
    H5Tinsert(type, "hi_j", HOFFSET(H5Box,hi_j), H5T_NATIVE_INT);
    H5Tinsert(type, "hi_k", HOFFSET(H5Box,hi_k), H5T_NATIVE_INT);

    return type;
};

void 
read_level(
    hid_t file,
    int lev,
    amrex::MultiFab& mf);

