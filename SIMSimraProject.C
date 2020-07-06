//==============================================================================
//!
//! \file SIMSimraProject.C
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for projection of SIMRA results.
//!
//==============================================================================

#include "SIMSimraProject.h"
#include "IFEM.h"
#include "ASMs3DSimra.h"
#include "Utilities.h"

#include <tinyxml.h>

#include <fstream>


SIMSimraProject::SIMSimraProject() :
  SIM3D(12)
{
}


ASMbase* SIMSimraProject::readPatch(std::istream &isp, int pchInd,
                                    const CharVec &unf, const char *whiteSpace) const
{
  ASMs3DSimra* pch = new ASMs3DSimra;
  if (pch)
  {
    if (!pch->read(isp) || this->getLocalPatchIndex(pchInd+1) < 1)
    {
      delete pch;
      pch = nullptr;
    }
    else
    {
      if (whiteSpace)
        IFEM::cout << whiteSpace <<"Reading patch "<< pchInd+1 << std::endl;
      if (checkRHSys && dynamic_cast<ASM3D*>(pch)->checkRightHandSystem())
        IFEM::cout <<"\tSwapped."<< std::endl;
      pch->idx = myModel.size();
    }
  }

  return pch;
}


bool SIMSimraProject::parse(const TiXmlElement *elem)
{
  if (!strcasecmp(elem->Value(), "resultfile")) {
    resultFile = utl::getValue(elem, "resultfile");
    IFEM::cout << "\tReading results from file " << resultFile << std::endl;
    return true;
  }

  return this->SIM3D::parse(elem);
}


bool SIMSimraProject::readResults()
{
  if (resultFile.empty())
    return false;

  std::ifstream ifs(resultFile);

  size_t nNode = this->getPatch(1)->getNoNodes();
  solution.resize(nNode*12);

  int header;
  Result cr(this->getNoNodes());
  ifs.read(reinterpret_cast<char*>(&header), 4);
  ifs.read(reinterpret_cast<char*>(&cr.time), sizeof(float));
  for (size_t i = 0; i < this->getNoNodes(); ++i) {
    ifs.read(reinterpret_cast<char*>(&cr.u1[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.u2[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.u3[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.ps[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.tk[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.td[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.vtef[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.pt[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.pts[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.rho[i]), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&cr.rhos[i]), sizeof(float));
  }
  ifs.read(reinterpret_cast<char*>(&header), 4);

  ifs.read(reinterpret_cast<char*>(&header), 4);
  ifs.read(reinterpret_cast<char*>(cr.strat.data()), nNode*sizeof(float));
  ifs.read(reinterpret_cast<char*>(&header), 4);

  size_t k = 1;
  for (size_t i = 0; i < nNode; ++i) {
    solution(k++) = cr.u1[i];
    solution(k++) = cr.u2[i];
    solution(k++) = cr.u3[i];
    solution(k++) = cr.ps[i];
    solution(k++) = cr.tk[i];
    solution(k++) = cr.td[i];
    solution(k++) = cr.vtef[i];
    solution(k++) = cr.pt[i];
    solution(k++) = cr.pts[i];
    solution(k++) = cr.rho[i];
    solution(k++) = cr.rhos[i];
    solution(k++) = cr.strat[i];
  }
  solTime = cr.time;

  return true;
}


//bool SIMSimraProject::writeSolFields()
//{

//}
