// $Id$
//==============================================================================
//!
//! \file main_CompareSetup.C
//!
//! \date Mar 17 2021
//!
//! \author Arne Morten Kvarving
//!
//! \brief Main program for comparing SIMRA initial and boundary data.
//!
//==============================================================================

#include "IFEM.h"
#include "ASMs3DSimra.h"
#include "LogStream.h"
#include "SimraIO.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <strings.h>
#include <vector>


/*!
 \brief Compares two init.dat files.
 \param file1 First file
 \param file2 Second file
 \param tol Tolerance in floating point comparisons
 \param npoint Number of nodes
 \param nelem Number of elements
*/

template<class Float>
bool compareInit(const char* file1, const char* file2, double tol, int npoint, int nelem)
{
  SimraIO<Float> init1(npoint), init2(npoint);
  std::vector<Float> elmPressure1(nelem), elmPressure2(nelem);

  using std::ifstream;
  auto&& readInit = [](SimraIO<Float>& init, std::vector<Float>& elmPressure, const char* file)
  {
    std::ifstream ifs(file);
    if (!ifs.good()) {
      std::cerr << "*** Failed to open init file " << file << std::endl;
      return false;
    }
    init.readInit(ifs, elmPressure);
    return true;
  };

  if (!readInit(init1, elmPressure1, file1))
    return false;
  if (!readInit(init2, elmPressure2, file2))
    return false;


  auto&& compareFloatVec = [tol](const std::vector<Float>& vec1,
                                 const std::vector<Float>& vec2,
                                 const char* name)
  {
    if (vec1.size() != vec2.size()) {
      std::cerr << "Vector \"" << name << "\" have different sizes, " << vec1.size() << " vs " << vec2.size() << std::endl;
      return false;
    }
    for (size_t i = 0; i < vec1.size(); ++i)
      if (std::fabs(vec1[i] - vec2[i]) > tol) {
        std::cerr << "Vector \"" << name << "\" difference at index " << i << " exceeds tolerance " << tol << std::endl;
        return false;
      }

    return true;
  };

  if (!compareFloatVec(elmPressure1, elmPressure2, "Pressure"))
    return false;
  if (!compareFloatVec(init1.u1, init2.u1, "u1"))
    return false;
  if (!compareFloatVec(init1.u2, init2.u2, "u2"))
    return false;
  if (!compareFloatVec(init1.u3, init2.u3, "u3"))
    return false;
  if (!compareFloatVec(init1.ps, init2.ps, "ps"))
    return false;
  if (!compareFloatVec(init1.tk, init2.tk, "tk"))
    return false;
  if (!compareFloatVec(init1.td, init2.td, "td"))
    return false;
  if (!compareFloatVec(init1.vtef, init2.vtef, "vtef"))
    return false;
  if (!compareFloatVec(init1.pt, init2.pt, "pt"))
    return false;
  if (!compareFloatVec(init1.pts, init2.pts, "pts"))
    return false;
  if (!compareFloatVec(init1.rho, init2.rho, "rho"))
    return false;
  if (!compareFloatVec(init1.rhos, init2.rhos, "rhos"))
    return false;
  if (!compareFloatVec(init1.strat, init2.strat, "strat"))
    return false;

  return true;
}


/*!
 \brief Compares two boun.dat files.
 \param file1 First file
 \param file2 Second file
 \param tol Tolerance in floating point comparisons
*/
bool compareBC(const char* file1, const char* file2, double tol)
{
  SimraBoundaryIO bc1, bc2;

  auto&& readFile = [](SimraBoundaryIO& bc, const char* file)
  {
    std::ifstream ifs(file);
    if (!ifs.good()) {
      std::cerr << "Failed to open boundary file " << file << std::endl;
      return false;
    }
    return bc.readFile(ifs, {});
  };

  if (!readFile(bc1, file1))
    return false;
  if (!readFile(bc2, file2))
    return false;

  auto&& compareFloatVec = [tol](const std::vector<double>& vec1,
                                 const std::vector<double>& vec2,
                                 const char* name)
  {
    if (vec1.size() != vec2.size()) {
      std::cerr << "Vector \"" << name << "\" have different sizes, "
                << vec1.size() << " vs " << vec2.size() << std::endl;
      return false;
    }
    for (size_t i = 0; i < vec1.size(); ++i)
      if (std::fabs(vec1[i] - vec2[i]) > tol) {
        std::cerr << "Vector \"" << name << "\" difference at index "
                  << i << " exceeds tolerance " << tol << std::endl;
        return false;
      }

    return true;
  };

  auto&& compareIntVec = [](const std::vector<int>& vec1,
                            const std::vector<int>& vec2,
                            const char* name)
  {
    if (vec1.size() != vec2.size()) {
      std::cerr << "Vector \"" << name << "\" have different sizes, "
                << vec1.size() << " vs " << vec2.size() << std::endl;
      return false;
    }
    for (size_t i = 0; i < vec1.size(); ++i)
      if (vec1[i] != vec2[i]) {
        std::cerr << "Vector \"" << name << "\" difference at index " << i
                  << ": " << vec1[i] << " vs " << vec2[i] << std::endl;
        return false;
      }

    return true;
  };

  auto&& compareData =
  [compareIntVec,compareFloatVec](const SimraBoundaryIO::Data& data1,
                                  const SimraBoundaryIO::Data& data2,
                                  const char* name)
  {
    return compareIntVec(data1.idx, data2.idx, (std::string(name) + ".idx").c_str()) &&
           compareFloatVec(data1.val, data2.val, (std::string(name) + ".val").c_str());
  };

  if (!compareFloatVec(bc1.z0, bc2.z0, "z0"))
    return false;
  if (!compareIntVec(bc1.wallE, bc2.wallE, "wallE"))
    return false;
  if (!compareIntVec(bc1.wall, bc2.wall, "wall"))
    return false;
  if (!compareIntVec(bc1.log, bc2.log, "log"))
    return false;
  if (!compareData(bc1.u, bc2.u, "u1"))
    return false;
  if (!compareData(bc1.v, bc2.v, "u2"))
    return false;
  if (!compareData(bc1.w, bc2.w, "u3"))
    return false;
  if (!compareData(bc1.p, bc2.p, "p"))
    return false;
  if (!compareData(bc1.e, bc2.e, "e"))
    return false;
  if (!compareData(bc1.k, bc2.k, "k"))
    return false;
  if (!compareData(bc1.d, bc2.d, "d"))
    return false;

  return true;
}


/*!
  \brief Main program for the comparing SIMRA initial and boundary data.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -tol     : Tolerance for comparing
  \arg -double  : Use double precision data
  \arg -tolBC   : Tolerance for comparing boundary data
  \arg -tolInit : Tolerance for comparing initial values
*/

int main (int argc, char** argv)
{
  if (argc < 6)
  {
    IFEM::cout <<"usage: "<< argv[0] <<" <initfile1> <initfile2> <boundaryfile1> <boundaryfile2> <mesh file>\n"
               << "\t[-double] - Use double precision input data\n"
               << "\t[-tol <tolerance>] - Comparison tolerance\n"
               << "\t[-tolBC <tolerance>] - Boundary data tolerance\n"
               << "\t[-tolInit <tolerance>] - Initial data tolerance\n";
    return 1;
  }

  double tolBC = 1e-5;
  double tolInit = 1e-5;
  std::vector<char*> file;
  for (int i = 1; i < argc; i++)
    if (!strcasecmp(argv[i],"-double"))
      ASMs3DSimra::useDouble = true;
    else if (!strcasecmp(argv[i],"-tol"))
      tolBC = tolInit = atof(argv[++i]);
    else if (!strcasecmp(argv[i],"-tolBC"))
      tolBC = atof(argv[i]);
    else if (!strcasecmp(argv[i],"-tolInit"))
      tolInit = atof(argv[i]);
    else
      file.push_back(argv[i]);

  if (file.size() != 5) {
    std::cerr << "*** Too many input files given\n";
    return 1;
  }

  IFEM::cout << "Reading mesh information from " << file[4] << std::endl;
  std::ifstream ifs(file[4]);
  int header, npoint, nelem;
  ifs.read(reinterpret_cast<char*>(&header), 4);
  ifs.read(reinterpret_cast<char*>(&npoint), 4);
  ifs.read(reinterpret_cast<char*>(&nelem), 4);
  ifs.close();

  IFEM::cout << "Comparing init files " << file[0] << " and " << file[1] << std::endl;
  bool res;
  if (ASMs3DSimra::useDouble)
    res = compareInit<double>(file[0], file[1], tolInit, npoint, nelem);
  else
    res = compareInit<float>(file[0], file[1], tolInit, npoint, nelem);

  if (!res)
    return 2;

  IFEM::cout << "Comparing boundary files " << file[2] << " and " << file[3] << std::endl;
  res = compareBC(file[2], file[3], tolBC);
  if (!res)
    return 3;

  return 0;
}
