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
#include "SIMSimraProject.h"
#include "Profiler.h"


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

  SIMSimraProject model;

  if (!model.read(argv[1]))
    return 1;

  if (!model.preprocess())
    return 2;

  if (!model.readResults())
    return 3;

  int nBlock = 0;
  int geoBlk = 0;
  model.writeGlvG(geoBlk, argv[1]);
  model.writeGlvS1(model.getSolution(), 1, nBlock, model.getSolutionTime());
  model.writeGlvStep(1, model.getSolutionTime(), 0);

  return 0;
}
