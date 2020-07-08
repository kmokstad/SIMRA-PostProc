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
  SIM3D(5)
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
  if (!strcasecmp(elem->Value(),"simra")) {
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(), "resultfile")) {
        resultFile = utl::getValue(child, "resultfile");
        IFEM::cout << "\tReading results from file " << resultFile << std::endl;
        return true;
      } else if (!strcasecmp(child->Value(),"scaling")) {
        utl::getAttribute(child, "u", uRef);
        utl::getAttribute(child, "L", lRef);
        IFEM::cout << "\tReference velocity U = " << uRef << std::endl;
        IFEM::cout << "\tReference length L = " << lRef << std::endl;
      }

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

  size_t nx, ny, nz;
  static_cast<const ASMs3DSimra*>(this->getPatch(1))->getNoStructNodes(nx,ny,nz);
  auto&& copySolution = [nx,ny,nz,this](const auto& cr, Vectors& solution)
  {
    size_t idx = 0;
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i, ++idx) {
          solution[0][idx] = cr.u1[k + i*nz + j*nz*nx] * uRef;
          solution[1][idx] = cr.u2[k + i*nz + j*nz*nx] * uRef;
          solution[2][idx] = cr.u3[k + i*nz + j*nz*nx] * uRef;
          solution[3][idx] = cr.ps[k + i*nz + j*nz*nx] * uRef * uRef;
          solution[4][idx] = cr.tk[k + i*nz + j*nz*nx] * uRef * uRef;
          solution[5][idx] = cr.td[k + i*nz + j*nz*nx] * pow(uRef,3) / lRef;
          solution[6][idx]= cr.vtef[k + i*nz + j*nz*nx] * uRef * lRef;
          solution[7][idx] = cr.pt[k + i*nz + j*nz*nx];
          solution[8][idx] = cr.pts[k + i*nz + j*nz*nx];
          solution[9][idx] = cr.rho[k + i*nz + j*nz*nx];
          solution[10][idx] = cr.rhos[k + i*nz + j*nz*nx];
          solution[11][idx] = cr.strat[k + i*nz + j*nz*nx];
        }
    return cr.time * lRef / uRef;
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
  if (solution.empty())
    return {};

  Matrix stmp(5, this->getPatch(1)->getNoNodes());
  stmp.fillRow(1, solution[0].data());
  stmp.fillRow(2, solution[1].data());
  stmp.fillRow(3, solution[2].data());
  stmp.fillRow(4, solution[3].data());
  stmp.fillRow(5, solution[6].data());

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
