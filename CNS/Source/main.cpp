#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Amr.H>
#include <AMReX_EB2.H>
#include "AMReX_EB2_DistributedGeometryShop.H"

#include <CNS.H>
#include "ReadChomboHDF5.H"

using namespace amrex;

amrex::LevelBld* getLevelBld ();
void initialize_EB2 (const Geometry& geom, const int required_level, const int max_level);

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
    const int finest_level = amr.finestLevel();

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

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    double timer_tot = amrex::second();
    double timer_init = 0.;
    double timer_advance = 0.;

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("stop_time",stop_time);
    }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

     // 1. Read SDF and Geometries from plotfile
    auto [vec_mf_chombo, vec_geom_chombo] = read_and_write_plotfile("plt_sphere_1lev", "plt_sphere_test", "SDF");

    {
        timer_init = amrex::second();

        Amr amr(getLevelBld());
        AmrLevel::SetEBSupportLevel(EBSupport::full);
        AmrLevel::SetEBMaxGrowCells(CNS::numGrow(),4,2);

        // vec_mf_chombo has already been read from the plotfile
        // and has no ghost cells.
        auto vec_mf_amrex = copy_to_amrex_hierarchy(vec_mf_chombo, amr);

        // Make a version of the sdf vector<multifab> with 1 ghost cell for 3d-interpolation to work
        auto mf_ghosted = make_ghosted_multifab(vec_mf_amrex, 1);

        // 1. Convert Vector<MultiFab> to Vector<std::shared_ptr<MultiFab>>

        amrex::Vector<amrex::MultiFab> mf_ptrs;

        for (auto& mf : mf_ghosted) {
            mf_ptrs.push_back(std::move(mf));
        }

        // 2. Create the Distributed SDF functor
        // Pass the vector of shared_ptrs
        amrex::EB2::DistributedSDF sdf_functor(std::move(mf_ptrs), std::move(vec_geom_chombo));

        // 3. Create the specialized shop
        auto gshop = amrex::EB2::makeDistributedShop(std::move(sdf_functor));

        // 4. Build the EB2 IndexSpace
        amrex::EB2::Build(gshop, amr.Geom(amr.maxLevel()), amr.maxLevel(), amr.maxLevel(), 4);

        //initialize_EB2(amr.Geom(amr.maxLevel()), amr.maxLevel(), amr.maxLevel());

        amr.init(strt_time,stop_time);

        timer_init = amrex::second() - timer_init;

        timer_advance = amrex::second();

        while ( amr.okToContinue() &&
                 (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time || stop_time < 0.0) )

        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
            //
            amr.coarseTimeStep(stop_time);
        }

        timer_advance = amrex::second() - timer_advance;

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }
    }

    timer_tot = amrex::second() - timer_tot;

    ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance},
                                      ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot     << "\n"
                   << "Run Time init         = " << timer_init    << "\n"
                   << "Run Time advance      = " << timer_advance << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
