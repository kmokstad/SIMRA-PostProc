// $Id$
//==============================================================================
//!
//! \file main_SolutionTransfer.C
//!
//! \date Mar 17 2021
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for solution transfer for SIMRA models.
//!
//==============================================================================

#include "ASMs3DSimra.h"
#include "SIMSimraTransfer.h"

#include "ASMbase.h"
#include "IFEM.h"
#include "MatVec.h"
#include "LogStream.h"
#include "Profiler.h"
#include "SIM3D.h"
#include "SIMadmin.h"
#include "SIMoptions.h"
#include "Vec3.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <strings.h>
#include <vector>



/*!
  \brief Main program for the SIMRA solution transfer program.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");
  char* infile = nullptr;
  bool selfTransfer = false;

  SIM3D::msgLevel = 1;

  IFEM::Init(argc,argv,"SIMRA-SolutionTransfer");
  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!infile)
      infile = argv[i];
    else if (!strcasecmp(argv[i],"-double"))
      ASMs3DSimra::useDouble = true;
    else if (!strcasecmp(argv[i],"-selfTransfer"))
      selfTransfer = true;
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

  SIMSimraTransfer input_model("simra_input");
  SIMSimraTransfer output_model("simra_output");

  if (!input_model.read(argv[1]) || !output_model.read(argv[1]))
    return 2;

  if (!input_model.preprocess() || !output_model.preprocess())
    return 3;

  if (input_model.hasNestingRegion()) {
      IFEM::cout << "\n>>> Determining nesting region";
      SIMSimraTransfer::NestedRegion region = input_model.determineNestingRegion(output_model);
      if (std::find(region[0].begin(), region[0].end(), -1) != region[0].end() ||
          std::find(region[1].begin(), region[1].end(), -1) != region[1].end()) {
        std::cerr << "*** Failed to detect nesting region." << std::endl;
        return 4;
      }
      IFEM::cout << "\n\tNesting region (" << region[0][0] << "," << region[0][1] << "," << region[0][2] << ')'
                 << " -> (" << region[1][0] << "," << region[1][1] << "," << region[1][2] << ")" << std::endl;

      size_t nx, ny, nz;
      static_cast<const ASMs3DSimra*>(input_model.getPatch(1))->getNoStructNodes(nx,ny,nz);

      if (region[0][2] != 1 || region[1][2] != static_cast<int>(nz)) {
        std::cerr << "*** Nesting region must span the entire z direction.\n";
        return 4;
      }

      output_model.setParams(region, input_model);
  }

  if (selfTransfer && input_model.getNoNodes() != output_model.getNoNodes()) {
    std::cerr << "** selfTransfer requested but models have different sizes, disabling.." << std::endl;
    selfTransfer = false;
  }

  Vec3 in_ll = input_model.getNodeCoord(1);
  Vec3 in_tr = input_model.getNodeCoord(input_model.getNoNodes());
  Vec3 out_ll = output_model.getNodeCoord(1);
  Vec3 out_tr = output_model.getNodeCoord(output_model.getNoNodes());

  IFEM::cout << "\nInput model:  ("
             << in_ll.x << ", " << in_ll.y << ", " << in_ll.z
             << ") -> ("
             << in_tr.x << ", " << in_tr.y << ", " << in_tr.z << ")\n"
             << "Output model: ("
             << out_ll.x << ", " << out_ll.y << ", " << out_ll.z
             << ") -> ("
             << out_tr.x << ", " << out_tr.y << ", " << out_tr.z << ")" << std::endl;

  utl::profiler->stop("Initialization");

  IFEM::cout << "\n>>> Reading input data from file '" << input_model.getResultFile() << "'" << std::endl;
  if (!input_model.readResults())
    return 5;

  Matrix inputs = input_model.getSolutions();

  if (selfTransfer) {
    output_model.setSolutions(input_model.getSolutions());
    IFEM::cout << "\tUsing original solution vectors since selfTransfer was requested." << std::endl;
  } else {
    IFEM::cout << "\tTransferring nodal data using linear interpolation" << std::endl;
    Vector output_vector;
    if (!output_model.getPatch(1)->evaluate(input_model.getPatch(1),
                                            inputs,
                                            output_vector)) {
      std::cerr << "*** Error transfering nodal data" << std::endl;
      return 6;
    }
    output_model.setSolutions(output_vector);
  }

  Vector elm_pressures;
  IFEM::cout << "\tTransferring elemental data using nearest match" << std::endl;
  output_model.elementTransfer(input_model);

  IFEM::cout << "<<< Writing output data to file '" << output_model.getResultFile() << "'" << std::endl;
  output_model.writeResults();

  if (!input_model.getBoundaryFile().empty()) {
    IFEM::cout << "\n>>> Reading boundary data from file '"
               << input_model.getBoundaryFile() << "'" << std::endl;
    if (!input_model.readBoundaryData()) {
      std::cerr << "*** Error reading boundary data\n";
      return 7;
    }

    IFEM::cout << "\n\tTransferring surface roughness using nearest match" << std::endl;
    Vector z0 = output_model.nodalTransfer(input_model);

    output_model.addTerrainBC(z0);
    if (output_model.useBoundaryTransfer() && !input_model.hasNestingRegion()) {
      IFEM::cout << "\tTransferring boundary data";
      output_model.boundaryTransfer(input_model);
    } else {
      IFEM::cout <<"\tClassifying boundary nodes";
      output_model.classifyBoundaryNodes();
    }
    output_model.addAtmosphereBC();

    IFEM::cout << "\n<<< Writing boundary data to file '"
               << output_model.getBoundaryFile() << "'" << std::endl;
    output_model.writeBoundaryData();

    if (input_model.opt.format > -1) {
      int nBlock_in = 0;
      int geoBlk_in = 0;
      std::string vtf_in = argv[1];
      vtf_in.replace(vtf_in.size()-5, 5, "_input.xinp");
      if (!input_model.writeGlvG(geoBlk_in, vtf_in.c_str())) {
        std::cerr << "Error writing VTF file for input model." << std::endl;
        return 8;
      }

      int nBlock_out = 0;
      int geoBlk_out = 0;
      if (!output_model.writeGlvG(geoBlk_out, argv[1])) {
        std::cerr << "Error writing VTF file for output model." << std::endl;
        return 9;
      }

      input_model.writeSolutionVectors(nBlock_in, {});
      output_model.writeSolutionVectors(nBlock_out, {});
      input_model.writeGlvBC(nBlock_in);
      output_model.writeGlvBC(nBlock_out);
      input_model.writeGlvStep(input_model.currentStep(), 0.0, 0);
      output_model.writeGlvStep(output_model.currentStep(), 0.0, 0);
    }
  }

  return 0;
}
