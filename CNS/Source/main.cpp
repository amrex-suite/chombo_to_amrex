
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Amr.H>
#include <AMReX_EB2.H>

#include <CNS.H>
#include "ReadChomboHDF5.H"

using namespace amrex;

amrex::LevelBld* getLevelBld ();
void initialize_EB2 (const Geometry& geom, const int required_level, const int max_level);

void read_and_write_plotfile(
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

    amrex::Vector<amrex::Geometry> geom(nlev);

    Array<int,AMREX_SPACEDIM> is_periodic
    {
        AMREX_D_DECL(0,0,0)
    };

    for (int lev = 0; lev < nlev; ++lev)
    {
        Geometry geom_lev(pf.probDomain(lev),
                          amrex::RealBox(
                            pf.probLo(),
                            pf.probHi()),
                          pf.coordSys(),
                          is_periodic
                         );

        geom[lev] = geom_lev;
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


           read_and_write_plotfile(
            "plt_sphere",
            "plt_sphere_test",
            "SDF"
        );

    {
        timer_init = amrex::second();

        Amr amr(getLevelBld());
        AmrLevel::SetEBSupportLevel(EBSupport::full);
        AmrLevel::SetEBMaxGrowCells(CNS::numGrow(),4,2);

        initialize_EB2(amr.Geom(amr.maxLevel()), amr.maxLevel(), amr.maxLevel());

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
