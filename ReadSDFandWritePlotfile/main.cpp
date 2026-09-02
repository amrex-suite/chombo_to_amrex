#include <hdf5.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "ReadChomboHDF5.H"

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    std::string chombo_hdf5_file;
    std::string amrex_plt_file;

    for (int i = 1; i < argc; ++i)
    {
	    std::string arg(argv[i]);

	    const std::string h5_prefix =
	        "--chombo-hdf5-file=";

	    const std::string plt_prefix =
	        "--amrex-plt-file=";

	    if (arg.rfind(h5_prefix, 0) == 0)
	    {
	        chombo_hdf5_file =
		    arg.substr(h5_prefix.size());
	    }
	    else if (arg.rfind(plt_prefix, 0) == 0)
	    {
	        amrex_plt_file =
		    arg.substr(plt_prefix.size());
	    }
	    else
	    {
	        amrex::Abort(
		    "Unknown command-line argument: " + arg);
	    }
    }

    if (chombo_hdf5_file.empty())
    {
	amrex::Abort(
	    "Missing required argument: "
	    "--chombo-hdf5-file=<filename>");
    }

    if (amrex_plt_file.empty())
    {
	amrex::Abort(
	    "Missing required argument: "
	    "--amrex-plt-file=<filename>");
    }

	hid_t file =
	    H5Fopen(
		chombo_hdf5_file.c_str(),
		H5F_ACC_RDONLY,
		H5P_DEFAULT);

	if (file < 0) {
	    amrex::Abort("Could not open HDF5 file");
	}

	// --------------------------------------------------------
	// Three AMR levels
	// --------------------------------------------------------

	const int nlev = 3;

	Vector<std::unique_ptr<MultiFab>> mf(nlev);
	Vector<Geometry> geom(nlev);
	Vector<IntVect> ref_ratio(nlev-1);

	// --------------------------------------------------------
	// Loop over levels
	// --------------------------------------------------------

	for (int lev = 0; lev < nlev; ++lev)
	{
	    std::string prefix =
		"level_" + std::to_string(lev);

	    // ====================================================
	    // Read dx
	    // ====================================================

	    double dx;

	    std::string dx_name = prefix + "/dx";

	    hid_t dx_attr =
		H5Aopen_by_name(
		    file,
		    prefix.c_str(),
		    "dx",
		    H5P_DEFAULT,
		    H5P_DEFAULT);

	    H5Aread(
		dx_attr,
		H5T_NATIVE_DOUBLE,
		&dx);

	    H5Aclose(dx_attr);

	    // ====================================================
	    // Read prob_domain
	    // ====================================================

	    hid_t pd_attr =
		H5Aopen_by_name(
		    file,
		    prefix.c_str(),
		    "prob_domain",
		    H5P_DEFAULT,
		    H5P_DEFAULT);

	    H5Box pd;

	    hid_t box_type = make_box_type();

	    H5Aread(
		pd_attr,
		box_type,
		&pd);

	    H5Tclose(box_type);
	    H5Aclose(pd_attr);

	    Box domain(
		IntVect(pd.lo_i, pd.lo_j, pd.lo_k),
		IntVect(pd.hi_i, pd.hi_j, pd.hi_k));

	    // ====================================================
	    // Read boxes
	    // ====================================================

	    std::string boxes_name =
		prefix + "/boxes";

	    hid_t boxes_dset =
		H5Dopen2(
		    file,
		    boxes_name.c_str(),
		    H5P_DEFAULT);

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

	    // ====================================================
	    // Construct BoxArray
	    // ====================================================

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

	    // ====================================================
	    // DistributionMapping
	    //
	    // We do NOT need to reproduce the original MPI
	    // distribution stored in "Processors".
	    // ====================================================

	    DistributionMapping dm(ba);

	    // ====================================================
	    // Geometry
	    // ====================================================

	    RealBox real_box(
		{0.0, 0.0, 0.0},
		{
		    dx * (domain.length(0)),
		    dx * (domain.length(1)),
		    dx * (domain.length(2))
		});

	    Vector<int> is_per(3, 0);

	    geom[lev] =
		Geometry(
		    domain,
		    &real_box,
		    CoordSys::cartesian,
		    is_per.data());

	    // ====================================================
	    // Allocate MultiFab
	    // ====================================================

	    mf[lev] =
		std::make_unique<MultiFab>(
		    ba,
		    dm,
		    1,   // one component
		    0);  // no ghost cells

	    mf[lev]->setVal(0.0);

	    // ====================================================
	    // Read SDF
	    // ====================================================

	    read_level(
		file,
		lev,
		*mf[lev]);

	    // ====================================================
	    // Refinement ratio
	    // ====================================================

	    if (lev < nlev-1)
	    {
		ref_ratio[lev] =
		    IntVect(4,4,4);
	    }
	}

	H5Fclose(file);

	// --------------------------------------------------------
	// Write AMReX plotfile
	// --------------------------------------------------------

	Vector<const MultiFab*> mf_ptrs(nlev);

	for (int lev = 0; lev < nlev; ++lev)
	{
	    mf_ptrs[lev] = mf[lev].get();
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

    amrex::Finalize();

    return 0;
}
