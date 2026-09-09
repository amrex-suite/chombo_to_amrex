#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <hdf5.h>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <filesystem>

#include "ReadChomboHDF5.H"

using namespace amrex;


// ============================================================
// Read Chombo HDF5, create AMReX MultiFabs/Geometry,
// write multilevel AMReX plotfile, and return the data.
// ============================================================

std::pair<Vector<MultiFab>, Vector<Geometry>>
read_chombo_hdf5_and_plot(const std::string& hdf5_file)
{
    std::filesystem::path hdf5_path(hdf5_file);

    std::string base_string = hdf5_path.stem().string();

    std::string amrex_plt_file = "plt_" + base_string;
    // --------------------------------------------------------
    // Open HDF5 file
    // --------------------------------------------------------

    hid_t file =
        H5Fopen(
            hdf5_file.c_str(),
            H5F_ACC_RDONLY,
            H5P_DEFAULT);

    if (file < 0)
    {
        Abort("Could not open HDF5 file: " + hdf5_file);
    }

    // --------------------------------------------------------
    // Determine number of AMR levels
    // --------------------------------------------------------

    int nlev = 0;

    while (true)
    {
        std::string level_name =
            "level_" + std::to_string(nlev);

        htri_t exists =
            H5Lexists(
                file,
                level_name.c_str(),
                H5P_DEFAULT);

        if (exists <= 0)
        {
            break;
        }

        ++nlev;
    }

    if (nlev == 0)
    {
        H5Fclose(file);
        Abort("No AMR levels found in Chombo HDF5 file");
    }

    // --------------------------------------------------------
    // Allocate output containers
    // --------------------------------------------------------

    Vector<MultiFab> mf(nlev);
    Vector<Geometry> geom(nlev);
    Vector<IntVect> ref_ratio(nlev - 1);

    // ========================================================
    // Loop over AMR levels
    // ========================================================

    for (int lev = 0; lev < nlev; ++lev)
    {
        std::string prefix =
            "level_" + std::to_string(lev);

        // ----------------------------------------------------
        // Read dx
        // ----------------------------------------------------

        double dx;

        hid_t dx_attr =
            H5Aopen_by_name(
                file,
                prefix.c_str(),
                "dx",
                H5P_DEFAULT,
                H5P_DEFAULT);

        if (dx_attr < 0)
        {
            H5Fclose(file);
            Abort("Could not open dx attribute for " + prefix);
        }

        H5Aread(
            dx_attr,
            H5T_NATIVE_DOUBLE,
            &dx);

        H5Aclose(dx_attr);

        // ----------------------------------------------------
        // Read prob_domain
        // ----------------------------------------------------

        hid_t pd_attr =
            H5Aopen_by_name(
                file,
                prefix.c_str(),
                "prob_domain",
                H5P_DEFAULT,
                H5P_DEFAULT);

        if (pd_attr < 0)
        {
            H5Fclose(file);
            Abort(
                "Could not open prob_domain attribute for "
                + prefix);
        }

        H5Box pd;

        hid_t box_type = make_box_type();

        H5Aread(
            pd_attr,
            box_type,
            &pd);

        H5Tclose(box_type);
        H5Aclose(pd_attr);

        Box domain(
            IntVect(
                pd.lo_i,
                pd.lo_j,
                pd.lo_k),
            IntVect(
                pd.hi_i,
                pd.hi_j,
                pd.hi_k));

        // ----------------------------------------------------
        // Read boxes
        // ----------------------------------------------------

        std::string boxes_name =
            prefix + "/boxes";

        hid_t boxes_dset =
            H5Dopen2(
                file,
                boxes_name.c_str(),
                H5P_DEFAULT);

        if (boxes_dset < 0)
        {
            H5Fclose(file);
            Abort("Could not open " + boxes_name);
        }

        hid_t boxes_space =
            H5Dget_space(boxes_dset);

        hsize_t dims[1];

        H5Sget_simple_extent_dims(
            boxes_space,
            dims,
            nullptr);

        int nboxes =
            static_cast<int>(dims[0]);

        std::vector<H5Box> h5boxes(nboxes);

        box_type = make_box_type();

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

        // ----------------------------------------------------
        // Construct BoxArray
        // ----------------------------------------------------

        BoxList box_list;

        for (int b = 0; b < nboxes; ++b)
        {
            Box bx(
                IntVect(
                    h5boxes[b].lo_i,
                    h5boxes[b].lo_j,
                    h5boxes[b].lo_k),
                IntVect(
                    h5boxes[b].hi_i,
                    h5boxes[b].hi_j,
                    h5boxes[b].hi_k));

            box_list.push_back(bx);
        }

        BoxArray ba(box_list);

        // ----------------------------------------------------
        // DistributionMapping
        //
        // We don't reproduce the original Chombo MPI
        // distribution.
        // ----------------------------------------------------

        DistributionMapping dm(ba);

        // ----------------------------------------------------
        // Geometry
        // ----------------------------------------------------

        RealBox real_box(
            {0.0, 0.0, 0.0},
            {
                dx * domain.length(0),
                dx * domain.length(1),
                dx * domain.length(2)
            });

        Vector<int> is_per(3, 0);

        geom[lev] =
            Geometry(
                domain,
                &real_box,
                CoordSys::cartesian,
                is_per.data());

        // ----------------------------------------------------
        // Allocate MultiFab
        // ----------------------------------------------------

        mf[lev].define(
            ba,
            dm,
            1,      // one component: SDF
            0);     // no ghost cells

        mf[lev].setVal(0.0);

        // ----------------------------------------------------
        // Read offsets
        // ----------------------------------------------------

        std::string offsets_name =
            prefix + "/data:offsets=0";

        hid_t offsets_dset =
            H5Dopen2(
                file,
                offsets_name.c_str(),
                H5P_DEFAULT);

        if (offsets_dset < 0)
        {
            H5Fclose(file);
            Abort("Could not open " + offsets_name);
        }

        hid_t offsets_space =
            H5Dget_space(offsets_dset);

        hsize_t offset_dims[1];

        H5Sget_simple_extent_dims(
            offsets_space,
            offset_dims,
            nullptr);

        int noffsets =
            static_cast<int>(offset_dims[0]);

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

        AMREX_ALWAYS_ASSERT(
            noffsets == nboxes + 1);

        // ----------------------------------------------------
        // Read complete SDF data array
        // ----------------------------------------------------

        std::string data_name =
            prefix + "/data:datatype=0";

        hid_t data_dset =
            H5Dopen2(
                file,
                data_name.c_str(),
                H5P_DEFAULT);

        if (data_dset < 0)
        {
            H5Fclose(file);
            Abort("Could not open " + data_name);
        }

        hid_t data_space =
            H5Dget_space(data_dset);

        hsize_t data_dims[1];

        H5Sget_simple_extent_dims(
            data_space,
            data_dims,
            nullptr);

        std::size_t ndata =
            static_cast<std::size_t>(
                data_dims[0]);

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

        // ----------------------------------------------------
        // Copy each Chombo box's data into its AMReX FArrayBox
        // ----------------------------------------------------

        for (int b = 0; b < nboxes; ++b)
        {
            Box const& bx =
                mf[lev].boxArray()[b];

            long long begin = offsets[b];
            long long end   = offsets[b + 1];

            long long expected =
                static_cast<long long>(
                    bx.numPts());

            AMREX_ALWAYS_ASSERT(
                end - begin == expected);

            const Real* src =
                data.data() + begin;

            Real* dst =
                mf[lev][b].dataPtr(0);

            std::copy(
                src,
                src + expected,
                dst);
        }

        // ----------------------------------------------------
        // Refinement ratio
        // ----------------------------------------------------

        if (lev < nlev - 1)
        {
            ref_ratio[lev] =
                IntVect(4, 4, 4);
        }
    }

    // --------------------------------------------------------
    // Close HDF5 file
    // --------------------------------------------------------

    H5Fclose(file);

    // ========================================================
    // Write AMReX multilevel plotfile
    // ========================================================

    Vector<const MultiFab*> mf_ptrs(nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        mf_ptrs[lev] = &mf[lev];
    }

    Vector<std::string> varnames = {"SDF"};

    Vector<int> level_steps(nlev, 0);

    WriteMultiLevelPlotfile(
        amrex_plt_file,
        nlev,
        mf_ptrs,
        varnames,
        geom,
        0.0,
        level_steps,
        ref_ratio);

    // ========================================================
    // Return MultiFabs and Geometry
    // ========================================================

    return {
        std::move(mf),
        std::move(geom)
    };
}

Vector<MultiFab>
make_ghosted_multifab(
    const Vector<MultiFab>& mf,
    int nghost)
{
    const int nlevels = mf.size();

    Vector<MultiFab> mf_ghosted(nlevels);

    for (int lev = 0; lev < nlevels; ++lev)
    {
        const MultiFab& src = mf[lev];

        mf_ghosted[lev].define(
            src.boxArray(),
            src.DistributionMap(),
            src.nComp(),
            nghost);

        // Copy the valid cells.
        MultiFab::Copy(
            mf_ghosted[lev],
            src,
            0,              // srccomp
            0,              // destcomp
            src.nComp(),    // number of components
            0);             // nghost

        // Fill ghost cells from neighboring boxes on this level.
        mf_ghosted[lev].FillBoundary();
    }

    return mf_ghosted;
}

std::pair<amrex::Vector<amrex::MultiFab>,
          amrex::Vector<amrex::Geometry>>
read_and_write_plotfile(
    const std::string& input_plotfile,
    const std::string& output_plotfile,
    const std::string& varname)
{
    amrex::PlotFileData pf(input_plotfile);

    const int finest_level = pf.finestLevel();
    const int nlev = finest_level + 1;

    // Check that variable exists.
    const auto& names = pf.varNames();

    auto it = std::find(names.begin(), names.end(), varname);

    if (it == names.end()) {
        amrex::Abort("Variable not found in plotfile: " + varname);
    }

    // Read the MultiFab at every AMR level.
    amrex::Vector<amrex::MultiFab> mfs(nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        mfs[lev] = pf.get(lev, varname);

        amrex::Print() << "Read level " << lev
                       << ": " << mfs[lev].boxArray().size()
                       << " boxes\n";
    }

    // Geometries.
    amrex::Vector<amrex::Geometry> geom(nlev);

    amrex::Array<int, AMREX_SPACEDIM> is_periodic{
        AMREX_D_DECL(0, 0, 0)
    };

    for (int lev = 0; lev < nlev; ++lev)
    {
        geom[lev] = amrex::Geometry(
            pf.probDomain(lev),
            amrex::RealBox(pf.probLo(), pf.probHi()),
            pf.coordSys(),
            is_periodic
        );
    }

    // Refinement ratios.
    amrex::Vector<amrex::IntVect> ref_ratio(nlev);

    for (int lev = 0; lev < finest_level; ++lev)
    {
        ref_ratio[lev] = pf.refRatioVect(lev);
    }

    // Level steps.
    amrex::Vector<int> level_steps(nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        level_steps[lev] = pf.levelStep(lev);
    }

    // MultiFab pointers.
    amrex::Vector<const amrex::MultiFab*> mf_ptrs(nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        mf_ptrs[lev] = &mfs[lev];
    }

    // Variable names.
    amrex::Vector<std::string> varnames = {varname};

    // Write new plotfile.
    amrex::WriteMultiLevelPlotfile(
        output_plotfile,
        nlev,
        mf_ptrs,
        varnames,
        geom,
        pf.time(),
        level_steps,
        ref_ratio
    );

    amrex::Print() << "Wrote plotfile: "
                   << output_plotfile << "\n";

    return {std::move(mfs), std::move(geom)};
}

amrex::Vector<amrex::MultiFab>
copy_to_amrex_hierarchy(
    const amrex::Vector<amrex::MultiFab>& vec_mf_chombo,
    Amr& amr)
{
    amrex::Print() << "amr.finestLevel() = "
               << amr.finestLevel() << "\n";

amrex::Print() << "vec_mf_chombo.size() = "
               << vec_mf_chombo.size() << "\n";
    const int finest_level = 0;//amr.finestLevel();

    if (static_cast<int>(vec_mf_chombo.size()) < finest_level + 1) {
        amrex::Abort(
            "vec_mf_chombo does not contain enough levels for the AMReX hierarchy");
    }

    amrex::Vector<amrex::MultiFab> vec_mf_amrex(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // AMReX-created hierarchy at this level.
        AmrLevel& level = amr.getLevel(lev);

        const amrex::BoxArray& ba = level.boxArray();
        const amrex::DistributionMapping& dm = level.DistributionMap();

        const amrex::MultiFab& src = vec_mf_chombo[lev];

        // Create destination MultiFab on AMReX's BoxArray and
        // DistributionMapping. No ghost cells for now.
        vec_mf_amrex[lev].define(
            ba,
            dm,
            src.nComp(),
            0
        );

        // Copy valid data from the source layout to the
        // AMReX-created layout.
        vec_mf_amrex[lev].ParallelCopy(src);

        amrex::Print() << "Copied level " << lev
                       << ": "
                       << src.boxArray().size() << " source boxes -> "
                       << ba.size() << " AMReX boxes\n";
    }

    return vec_mf_amrex;
}

