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
#include "DataExporter.h"
#include "HDF5Writer.h"
#include "SIMenums.h"
#include "SIMFSWallDistance.h"
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

  SIM3D::msgLevel = 1;

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
    IFEM::cout <<"usage: "<< argv[0] <<" <inputfile>\n"
               << "\t[-double] - Use double precision input data\n"
               << "\t[-vtf <format>] - Write VTF result file\n"
               <<"\t[-hdf5] - Write HDF5 result file\n";
    return 1;
  }

  SIMSimraProject model;

  if (!model.read(argv[1]))
    return 2;

  if (!model.preprocess())
    return 3;

  utl::profiler->stop("Initialization");

  int nBlock = 0;
  int geoBlk = 0;
  if (model.opt.format > -1)
    if (!model.writeGlvG(geoBlk, argv[1])) {
      std::cerr << "Error writing VTF file." << std::endl;
      return 4;
    }

  Vector sol;
  Vectors gNorm;
  Matrix eNorm;
  Vectors projs;

  std::unique_ptr<DataExporter> exporter;
  if (model.opt.dumpHDF5(infile)) {
    exporter.reset(new DataExporter(true, 1));
    exporter->registerWriter(new HDF5Writer(model.opt.hdf5, model.getProcessAdm()));
  }

  int iStep = 1;
  std::vector<std::string> prefix;
  while (model.readResults()) {

    IFEM::cout << "\n  step=" << iStep << "  time=" << model.getSolutionTime() << std::endl;
    if (exporter && iStep == 1)
      model.registerFields(*exporter, sol, projs, eNorm);

    if (iStep == 1) {
      Vector& dist = model.getDistance();
      if (!dist.empty()) {
        if (!model.orthogonalDistance()) {
          if (!SIMFSWallDistance::solveProblem(dist, argv[1])) {
            std::cerr << "Error solving terrain distance problem." << std::endl;
            return 6;
          }
        }
      }
    }

    // Project the secondary solution
    sol = model.getSolution();
    size_t idx = 0;
    model.setMode(SIM::RECOVERY);
    projs.resize(model.opt.project.size());
    TimeDomain time;
    time.t = model.getSolutionTime();
    for (const SIMoptions::ProjectionMap::value_type& prj : model.opt.project) {
      if (prj.first <= SIMoptions::NONE)
        idx++; // No projection for this norm group
      else {
        Matrix stmp(projs[idx++]);
        if (!model.project(stmp,sol,prj.first,time))
          return 5;
      }
      if (iStep == 1)
        prefix.push_back(prj.second);
    }
    if (exporter)
      exporter->setNormPrefixes(prefix);

    if (!model.solutionNorms(time, Vectors(1,sol), projs, gNorm, &eNorm)) {
      std::cerr << "Error calculating solution norms." << std::endl;
      return 6;
    }

    model.printSolutionNorms(gNorm);

    if (model.opt.format > -1) {
      if (!model.writeSolutionVectors(nBlock, sol)) {
        std::cerr << "Error writing solution fields to VTF" << std::endl;
        return 7;
      }

      idx = 1;
      for (const SIMoptions::ProjectionMap::value_type& prj : model.opt.project) {
        if (prj.first > SIMoptions::NONE)
          if (!model.writeGlvP(projs[idx-1], 1, nBlock,
                               100+model.getProblem()->getNoFields(2)*idx, prj.second.c_str())) {
            std::cerr << "Error writing projection to VTF" << std::endl;
            return 8;
          }
        ++idx;
      }

      if (!model.writeGlvN(eNorm, iStep, nBlock, prefix)) {
        std::cerr << "Error writing norms to VTF" << std::endl;;
        return 9;
      }

      if (iStep == 1 && !model.getDistance().empty())
        model.writeGlvS(model.getDistance(), "terrain distance", iStep, nBlock, 110);

      model.writeGlvStep(iStep, model.getSolutionTime(), 0);
    }
    if (exporter)
      exporter->dumpTimeLevel();

    ++iStep;
  }

  return 0;
}
