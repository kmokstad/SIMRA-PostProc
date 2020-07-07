// $Id$
//==============================================================================
//!
//! \file main_PostProc.C
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for postprocessing for SIMRA results.
//!
//==============================================================================

#include "IFEM.h"
#include "ASMs3DSimra.h"
#include "SIMenums.h"
#include "SIMSimraProject.h"
#include "Profiler.h"
#include "SIMargsBase.h"


/*!
  \brief Main program for the SIMRA post processor.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");
  char* infile = nullptr;

  IFEM::Init(argc,argv,"SIMRA-PostProc");
  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!infile)
      infile = argv[i];
    else if (!strcasecmp(argv[i],"-double"))
      ASMs3DSimra::useDouble = true;
   else
     std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    IFEM::cout <<"usage: "<< argv[0]
               <<" <inputfile> [-double]\n"
               <<"             [-hdf5]\n"
               <<"             [-vtf <format> [-nviz <nviz>]\n";
    return 1;
  }

  SIMSimraProject model;

  if (!model.read(argv[1]))
    return 2;

  if (!model.preprocess())
    return 3;

  if (!model.readResults())
    return 4;

  utl::profiler->stop("Initialization");

  Vector sol = model.getSolution();

  // Project the secondary solution
  size_t idx = 0;
  model.setMode(SIM::RECOVERY);
  Vectors projs;
  projs.resize(model.opt.project.size());
  TimeDomain time;
  time.t = model.getSolutionTime();
  std::vector<std::string> prefix;
  for (const SIMoptions::ProjectionMap::value_type& prj : model.opt.project) {
    if (prj.first <= SIMoptions::NONE)
      idx++; // No projection for this norm group
    else {
      Matrix stmp(projs[idx++]);
      if (!model.project(stmp,sol,prj.first,time))
        return 5;
    }
    prefix.push_back(prj.second);
  }

  Vectors gNorm;
  Matrix eNorm;
  if (!model.solutionNorms(time, Vectors(1,sol), projs, gNorm, &eNorm))
    return 7;

  model.printSolutionNorms(gNorm);

  if (model.opt.format > -1) {
    int nBlock = 0;
    int geoBlk = 0;
    model.writeGlvG(geoBlk, argv[1]);
    model.writeSolutionVectors(nBlock);
    model.writeGlvS2(sol, 1, nBlock, model.getSolutionTime());
    idx = 1;
    for (const SIMoptions::ProjectionMap::value_type& prj : model.opt.project) {
      if (prj.first > SIMoptions::NONE)
        if (!model.writeGlvP(projs[idx-1], 1, nBlock,
                             100+model.getProblem()->getNoFields(2)*idx, prj.second.c_str()))
          return 6;
      ++idx;
    }

    if (!model.writeGlvN(eNorm, 1, nBlock, prefix))
      return 8;

    model.writeGlvStep(1, model.getSolutionTime(), 0);
  }

  return 0;
}
