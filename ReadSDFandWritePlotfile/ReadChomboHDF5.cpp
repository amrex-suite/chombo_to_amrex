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

        long long expected =
            static_cast<long long>(bx.numPts());

        AMREX_ALWAYS_ASSERT(
            end - begin == expected);

        const Real* src =
            data.data() + begin;

        auto& fab = mf[b];

        Real* dst = fab.dataPtr(0);

        std::copy(
            src,
            src + expected,
            dst);
    }
}
