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

    ParmParse pp("eb2");
    std::string chombo_sdf_hdf5_file = "";
    pp.query("chombo_sdf_hdf5_file", chombo_sdf_hdf5_file);

    {
        timer_init = amrex::second();

        Amr amr(getLevelBld());
        AmrLevel::SetEBSupportLevel(EBSupport::full);
        AmrLevel::SetEBMaxGrowCells(CNS::numGrow(),4,2);

       // 1. Read SDF and Geometries from plotfile
        //auto [vec_mf_chombo, vec_geom_chombo] = read_and_write_plotfile("plt_sphere_1lev", "plt_sphere_test", "SDF");

        auto [vec_mf_chombo, vec_geom_chombo] = read_chombo_hdf5_and_plot(chombo_sdf_hdf5_file);

        // Make a version of the sdf vector<multifab> with 1 ghost cell for 3d-interpolation to work
        auto vec_mf_chombo_ghosted = make_ghosted_multifab(vec_mf_chombo, 1);

        // 2. Create the Distributed SDF functor
        // Pass the vector of shared_ptrs
        amrex::EB2::DistributedSDF sdf_functor(vec_mf_chombo_ghosted, std::move(vec_geom_chombo));

        // 3. Create the specialized shop
        auto gshop = amrex::EB2::makeDistributedShop(std::move(sdf_functor));

        // 4. Build the EB2 IndexSpace
        amrex::EB2::Build(gshop, amr.Geom(amr.maxLevel()), amr.maxLevel(), amr.maxLevel(), 4, false);

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
