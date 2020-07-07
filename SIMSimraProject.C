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
  SIM3D(4)
{
  myProblem = &itg;
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
  solution.resize(12);
  for (Vector& sol : solution)
    sol.resize(this->getNoNodes());

  auto&& copySolution = [](const auto& cr, Vectors& solution)
  {
    for (size_t i = 0; i < cr.u1.size(); ++i) {
      solution[0][i] = cr.u1[i];
      solution[1][i] = cr.u2[i];
      solution[2][i] = cr.u3[i];
      solution[3][i] = cr.ps[i];
      solution[4][i] = cr.tk[i];
      solution[5][i] = cr.td[i];
      solution[6][i]= cr.vtef[i];
      solution[7][i] = cr.pt[i];
      solution[8][i] = cr.pts[i];
      solution[9][i] = cr.rho[i];
      solution[10][i] = cr.rhos[i];
      solution[11][i] = cr.strat[i];
    }
    return cr.time;
  };

  if (ASMs3DSimra::useDouble) {
    Result<double> cr(this->getNoNodes());
    cr.read(ifs);
    solTime = copySolution(cr, solution);
  } else {
    Result<float> cr(this->getNoNodes());
    cr.read(ifs);
    solTime = copySolution(cr, solution);
  }

  return true;
}


bool SIMSimraProject::writeSolutionVectors(int& nBlock)
{
  static constexpr const char* names[] =
    {"u_x", "u_y", "u_z", "ps", "tk",
     "td", "vtef", "pt", "pts", "rho", "rhos", "strat"};

  bool result = true;
  for (size_t i = 0; i < 12; ++i)
    result &= this->writeGlvS(solution[i], names[i], 1, nBlock, i+2);

  return result;
}


Vector SIMSimraProject::getSolution() const
{
  Matrix stmp(4, this->getPatch(1)->getNoNodes());
  stmp.fillRow(1, solution[0].data());
  stmp.fillRow(2, solution[1].data());
  stmp.fillRow(3, solution[2].data());
  stmp.fillRow(4, solution[3].data());

  return stmp;
}


void SIMSimraProject::printSolutionNorms(const Vectors& gNorm) const
{
  if (gNorm.empty()) return;

  IFEM::cout << "\n>>> Norm summary <<<";
  size_t w = 36;
  IFEM::cout << "\n  L2 norm |p^h| = (p^h,p^h)^0.5"
             << utl::adjustRight(w-31,"") << gNorm[0](4);
  IFEM::cout << "\n  L2 norm |u^h| = (u^h,u^h)^0.5"
             << utl::adjustRight(w-31,"")<< gNorm[0](1);
  IFEM::cout << "\n  H1 norm |u^h| = a(u^h,u^h)^0.5"
             << utl::adjustRight(w-32,"") << gNorm[0](2);
  IFEM::cout << "\n  L2 norm |div u^h|"
             << utl::adjustRight(w-19,"") << gNorm[0](3);

  size_t j = 0;
  for (const auto& prj : this->opt.project)
    if (++j < gNorm.size())
      this->printNormGroup(gNorm[j],gNorm[0],prj.second);
}


void SIMSimraProject::printNormGroup(const Vector& rNorm,
                                     const Vector& fNorm,
                                     const std::string& name) const
{
  IFEM::cout << "\nError estimates based on >>> " << name << " <<<";
  size_t w = 36;
  IFEM::cout << "\n  H1 norm |u^*-u^h|"
             << utl::adjustRight(w-19,"") << rNorm(1)
             << "\n  L2 norm |p^*-p^h|"
             << utl::adjustRight(w-19,"") << rNorm(4)
             << "\n  Relative error: " << rNorm(1) / hypot(fNorm(2), rNorm(1))*100.0
             << "%" << std::endl;

}
