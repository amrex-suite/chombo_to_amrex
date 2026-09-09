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

#include "ReadChomboHDF5.H"

using namespace amrex;

void read_level(
    hid_t file,
    int lev,
    MultiFab& mf)
{
    std::string prefix = "level_" + std::to_string(lev);

    // ------------------------------------------------------------
    // Read boxes
    // ------------------------------------------------------------

    std::string boxes_name = prefix + "/boxes";

    hid_t boxes_dset = H5Dopen2(file, boxes_name.c_str(), H5P_DEFAULT);

    if (boxes_dset < 0) {
        Abort("Could not open boxes dataset");
    }

    hid_t boxes_space = H5Dget_space(boxes_dset);

    hsize_t dims[1];
    H5Sget_simple_extent_dims(boxes_space, dims, nullptr);

    int nboxes = static_cast<int>(dims[0]);

    std::vector<H5Box> h5boxes(nboxes);

    hid_t box_type = make_box_type();

    H5Dread(
        boxes_dset,
        box_type,
        H5S_ALL,
        H5S_ALL,
        H5P_DEFAULT,
        h5boxes.data());

    H5Tclose(box_type);
    H5Sclose(boxes_space);
    H5Dclose(boxes_dset);

    // ------------------------------------------------------------
    // Read offsets
    // ------------------------------------------------------------

    std::string offsets_name = prefix + "/data:offsets=0";

    hid_t offsets_dset =
        H5Dopen2(file, offsets_name.c_str(), H5P_DEFAULT);

    hid_t offsets_space = H5Dget_space(offsets_dset);

    hsize_t offset_dims[1];
    H5Sget_simple_extent_dims(
        offsets_space, offset_dims, nullptr);

    int noffsets = static_cast<int>(offset_dims[0]);

    std::vector<long long> offsets(noffsets);

    H5Dread(
        offsets_dset,
        H5T_NATIVE_LLONG,
        H5S_ALL,
        H5S_ALL,
        H5P_DEFAULT,
        offsets.data());

    H5Sclose(offsets_space);
    H5Dclose(offsets_dset);

    AMREX_ALWAYS_ASSERT(noffsets == nboxes + 1);

    // ------------------------------------------------------------
    // Read complete SDF data array
    // ------------------------------------------------------------

    std::string data_name =
        prefix + "/data:datatype=0";

    hid_t data_dset =
        H5Dopen2(file, data_name.c_str(), H5P_DEFAULT);

    hid_t data_space = H5Dget_space(data_dset);

    hsize_t data_dims[1];
    H5Sget_simple_extent_dims(
        data_space, data_dims, nullptr);

    std::size_t ndata =
        static_cast<std::size_t>(data_dims[0]);

    std::vector<Real> data(ndata);

    H5Dread(
        data_dset,
        H5T_NATIVE_DOUBLE,
        H5S_ALL,
        H5S_ALL,
        H5P_DEFAULT,
        data.data());

    H5Sclose(data_space);
    H5Dclose(data_dset);

    // ------------------------------------------------------------
    // Copy each FArrayBox into the MultiFab
    // ------------------------------------------------------------

    for (int b = 0; b < nboxes; ++b)
    {
        Box const& bx = mf.boxArray()[b];

        long long begin = offsets[b];
        long long end   = offsets[b+1];

const int nx = bx.length(0);
const int ny = bx.length(1);
const int nz = bx.length(2);

const int ng = 0;

// Chombo FArrayBox contains ghost cells:
// (nx + 2*ng) x (ny + 2*ng) x (nz + 2*ng)
const long long chombo_size =
    static_cast<long long>(nx + 2*ng) *
    (ny + 2*ng) *
    (nz + 2*ng);

AMREX_ALWAYS_ASSERT(
    end - begin == chombo_size);

const Real* src =
    data.data() + begin;

auto& fab = mf[b];

Real* dst = fab.dataPtr(0);

// Copy only the valid 128^3 region.
for (int k = 0; k < nz; ++k)
{
    for (int j = 0; j < ny; ++j)
    {
        const Real* src_row =
            src +
            (k + ng) * (ny + 2*ng) * (nx + 2*ng) +
            (j + ng) * (nx + 2*ng) +
            ng;

        Real* dst_row =
            dst +
            k * ny * nx +
            j * nx;

        std::copy(
            src_row,
            src_row + nx,
            dst_row);
    }
}
    }
}
