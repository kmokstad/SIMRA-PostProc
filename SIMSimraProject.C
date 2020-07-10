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
#include "DataExporter.h"
#include "ASMs3DSimra.h"
#include "IFEM.h"
#include "Utilities.h"

#include <tinyxml.h>


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
  if (!strcasecmp(elem->Value(),"simra")) {
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(), "resultfile")) {
        resultFile = utl::getValue(child, "resultfile");
        std::string type;
        utl::getAttribute(child, "type", type);
        if (type == "history")
          rType = HISTORY_FILE;
        else
          rType = RESTART_FILE;
        IFEM::cout << "\tReading results from "
                   << (rType == HISTORY_FILE?"history":"restart")
                   << " file " << resultFile << std::endl;
      } else if (!strcasecmp(child->Value(),"scaling")) {
        utl::getAttribute(child, "u", uRef);
        utl::getAttribute(child, "L", lRef);
        IFEM::cout << "\tReference velocity U = " << uRef << std::endl;
        IFEM::cout << "\tReference length L = " << lRef << std::endl;
      } else if (!strcasecmp(child->Value(),"anasol")) {
          std::string type;
          utl::getAttribute(child,"type",type,true);
          if (type == "expression") {
            IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
            mySol = new AnaSol(child,false);
          } else
            std::cerr <<"  ** SIMSimraProject::parse: Invalid analytical solution "
                         << type <<" (ignored)"<< std::endl;
      }

    return true;
  }


  return this->SIM3D::parse(elem);
}


bool SIMSimraProject::readResults()
{
  if (resultFile.empty()) {
    std::cerr << "No result file specified." << std::endl;
    return false;
  }

  if (rType == RESTART_FILE && iStep > 1)
    return false;

  if (!ifs.is_open()) {
    ifs.open(resultFile);
    if (!ifs.good()) {
      std::cerr << "Error opening results file " << resultFile << std::endl;
      return false;
    }
    ifs.seekg(0, std::ios_base::end);
    fileSize = ifs.tellg();
    ifs.seekg(0, std::ios_base::beg);
  }

  if (!ifs.good() || ifs.eof() || ifs.tellg() == fileSize)
    return false;

  if (solution.empty()) {
    solution.resize(12);
    for (Vector& sol : solution)
      sol.resize(this->getNoNodes());
  }

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

  auto&& copyElmPressure = [nx,ny,nz](const auto& in,
                                      std::vector<double>& out)
  {
    size_t idx = 0;
    for (size_t k = 0; k < nz-1; ++k)
      for (size_t j = 0; j < ny-1; ++j)
        for (size_t i = 0; i < nx-1; ++i, ++idx) {
          out[idx] = in[k + i*(nz-1) + j*(nz-1)*(nx-1)];
        }
  };

  if (rType == HISTORY_FILE && iStep == 0) {
    int header, size;
    ifs.read(reinterpret_cast<char*>(&header), 4);
    ifs.read(reinterpret_cast<char*>(&size), 4);
    ifs.read(reinterpret_cast<char*>(&header), 4);
    if (size == 8)
      ASMs3DSimra::useDouble = true;
    itg.elmPressure.resize(this->getPatch(1)->getNoElms());
  }

  if (ASMs3DSimra::useDouble) {
    Result<double> cr(this->getNoNodes());
    if (rType == RESTART_FILE)
      cr.readRestart(ifs);
    else {
      std::vector<double> elmPressure(this->getPatch(1)->getNoElms());
      cr.readHistory(ifs, elmPressure);
      copyElmPressure(elmPressure, itg.elmPressure);
    }
    solTime = copySolution(cr, solution);
  } else {
    Result<float> cr(this->getNoNodes());
    if (rType == RESTART_FILE)
      cr.readRestart(ifs);
    else {
      std::vector<float> elmPressure(this->getPatch(1)->getNoElms());
      cr.readHistory(ifs, elmPressure);
      copyElmPressure(elmPressure, itg.elmPressure);
    }
    solTime = copySolution(cr, solution);
  }

  ++iStep;
  return true;
}


bool SIMSimraProject::writeSolutionVectors(int& nBlock, const Vector& sol)
{
  static constexpr const char* names[] =
    {"u_x", "u_y", "u_z", "ps", "tk",
     "td", "vtef", "pt", "pts", "rho", "rhos", "strat"};

  bool result = true;
  for (size_t i = 0; i < 12; ++i) {
    result &= this->writeGlvS(solution[i], names[i], iStep, nBlock, i+2);
    if ( i == 3 && !itg.elmPressure.empty())
      result &= this->writeGlvE(itg.elmPressure, iStep, nBlock, "p");
  }
  result &= this->writeGlvS2(sol, iStep, nBlock, this->getSolutionTime());

  return result;
}


Vector SIMSimraProject::getSolution() const
{
  if (solution.empty())
    return {};

  Matrix stmp(4, this->getPatch(1)->getNoNodes());
  stmp.fillRow(1, solution[0].data());
  stmp.fillRow(2, solution[1].data());
  stmp.fillRow(3, solution[2].data());
  stmp.fillRow(4, solution[6].data());

  return stmp;
}


 bool SIMSimraProject::postProcessNorms(Vectors& gNorm, Matrix* eNormp)
 {
   if (!eNormp)
     return true;

   NormBase* norm = itg.getNormIntegrand(mySol);
   if (!norm)
     return false;

   size_t ip = norm->getNoFields(1)+1;
   delete norm;

   Matrix& eNorm = *eNormp;
   for (size_t i = 1; i <= eNorm.cols(); ++ i) {
     eNorm(ip + SimraNorm::EFF_REC_VEL, i) = eNorm(i, ip+SimraNorm::H1_Ur_Uh) /
                                             eNorm(i, 1+SimraNorm::H1_E_U);
     eNorm(ip + SimraNorm::EFF_REC_PRESS, i) = eNorm(i, ip+SimraNorm::L2_Pr_Ph) /
                                               eNorm(i, 1+SimraNorm::L2_E_P);
     eNorm(ip + SimraNorm::EFF_REC_STRESS, i) = eNorm(i, ip+SimraNorm::L2_SIGMA_REC) /
                                                eNorm(i, 1+SimraNorm::L2_E_SIGMA);
     double tRel = hypot(eNorm(i, 1+SimraNorm::H1_E_U), eNorm(i, 1+SimraNorm::L2_E_P));
     if (tRel == 0)
       tRel = 1.0;
     eNorm(ip + SimraNorm::EFF_REC_TOTAL, i) = hypot(eNorm(i, ip+SimraNorm::H1_Ur_Uh), eNorm(i, ip+SimraNorm::L2_Pr_Ph)) / tRel;
   }

   return true;
 }


void SIMSimraProject::printSolutionNorms(const Vectors& gNorm) const
{
  if (gNorm.empty()) return;

  IFEM::cout << "\n>>> Norm summary <<<";
  size_t w = 36;
  IFEM::cout << "\n  L2 norm |p^h| = (p^h,p^h)^0.5"
             << utl::adjustRight(w-31,"") << gNorm[0][SimraNorm::L2_Ph]
             << "\n  L2 norm |u^h| = (u^h,u^h)^0.5"
             << utl::adjustRight(w-31,"")<< gNorm[0][SimraNorm::L2_Uh]
             << "\n  H1 norm |u^h| = a(u^h,u^h)^0.5"
             << utl::adjustRight(w-32,"") << gNorm[0][SimraNorm::H1_Uh]
             << "\n  L2 norm |div u^h|"
             << utl::adjustRight(w-19,"") << gNorm[0][SimraNorm::L2_DIV_Uh]
             << "\n  L2 norm |s^h| = (s^h,s^h)^0.5"
             << utl::adjustRight(w-31,"") << gNorm[0][SimraNorm::L2_SIGMAh];
  if (mySol)
    this->printExactNorms(gNorm[0], w);

  size_t j = 0;
  for (const auto& prj : this->opt.project)
    if (++j < gNorm.size())
      this->printNormGroup(gNorm[j],gNorm[0],prj.second);
}


void SIMSimraProject::printExactNorms(const Vector& gNorm, size_t w) const
{
  if (!mySol)
    return;

  if (mySol->getScalarSol()) {
    IFEM::cout << "\n  L2 norm |p| = (p,p)^0.5"
               << utl::adjustRight(w-25,"") << gNorm[SimraNorm::L2_P]
               << "\n  L2 norm |e| = (e,e)^0.5, e=p-p^h"
               << utl::adjustRight(w-34,"") << gNorm[SimraNorm::L2_E_P];
  }
  if (mySol->getVectorSol()) {
    IFEM::cout << "\n  L2 norm |u| = (u,u)^0.5"
               << utl::adjustRight(w-25,"") << gNorm[SimraNorm::L2_U]
               << "\n  L2 norm |e| = (e,e)^0.5, e=u-u^h"
               << utl::adjustRight(w-34,"") << gNorm[SimraNorm::L2_E_U];
  }
  if (mySol->getVectorSecSol()) {
    IFEM::cout << "\n  H1 norm |u| = a(u,u)^0.5"
               << utl::adjustRight(w-26,"") << gNorm[SimraNorm::H1_U]
               << "\n  H1 norm |e| = a(e,e)^0.5, e=u-u^h"
               << utl::adjustRight(w-35,"") << gNorm[SimraNorm::H1_E_U];
  }
  if (mySol->getVectorSecSol() && mySol->getScalarSol()) {
    IFEM::cout << "\n  L2 norm |sigma|"
               << utl::adjustRight(w-17,"") << gNorm[SimraNorm::L2_SIGMA]
               << "\n  L2 norm |e| = (e,e)^0.5, e=s-s^h"
               << utl::adjustRight(w-34,"") << gNorm[SimraNorm::L2_E_SIGMA];
  }
  if (mySol->getScalarSol()) {
    double pRel = gNorm[SimraNorm::L2_P];
    if (pRel == 0.0)
      pRel = 1.0;
    IFEM::cout << "\n  Exact pressure error (%)"
               << utl::adjustRight(w-26,"")
               << gNorm[SimraNorm::L2_E_P] / pRel * 100.0;
  }
  if (mySol->getVectorSol() && gNorm[SimraNorm::L2_U] != 0.0)
    IFEM::cout << "\n  Exact L2 velocity error (%)"
               << utl::adjustRight(w-29,"")
               << gNorm[SimraNorm::L2_E_U] / gNorm[SimraNorm::L2_U] * 100.0;
  if (mySol->getVectorSecSol() && gNorm[SimraNorm::H1_U] != 0.0)
    IFEM::cout << "\n  Exact H1 velocity error (%)"
               << utl::adjustRight(w-29,"")
               << gNorm[SimraNorm::H1_E_U] / gNorm[SimraNorm::H1_U] * 100.0;
  if (mySol->getVectorSecSol() && mySol->getScalarSol() && gNorm[SimraNorm::L2_SIGMA] != 0.0)
    IFEM::cout << "\n  Exact stress error (%)"
               << utl::adjustRight(w-24,"")
               << gNorm[SimraNorm::L2_E_SIGMA] / gNorm[SimraNorm::L2_SIGMA] * 100.0;
  if (mySol->getVectorSecSol() && mySol->getScalarSol() &&
      (gNorm[SimraNorm::L2_P] != 0.0 || gNorm[SimraNorm::H1_U] != 0.0))
    IFEM::cout << "\n  Exact total error (%)"
               << utl::adjustRight(w-23,"")
               << hypot(gNorm[SimraNorm::L2_E_P], gNorm[SimraNorm::H1_E_U]) /
                  hypot(gNorm[SimraNorm::L2_P], gNorm[SimraNorm::H1_U]) * 100.0;
}


void SIMSimraProject::printNormGroup(const Vector& rNorm,
                                     const Vector& fNorm,
                                     const std::string& name) const
{
  IFEM::cout << "\nError estimates based on >>> " << name << " <<<";
  size_t w = 36;
  double uRel;
  if (mySol && mySol->getVectorSecSol())
    uRel = fNorm[SimraNorm::H1_U];
  else
    uRel = hypot(fNorm[SimraNorm::H1_Uh], rNorm[SimraNorm::H1_Ur_Uh]);
  double pRel;
  if (mySol && mySol->getScalarSol())
    pRel = fNorm[SimraNorm::L2_P];
  else
    pRel = hypot(fNorm[SimraNorm::L2_Ph], rNorm[SimraNorm::L2_Pr_Ph]);
  if (pRel == 0.0)
    pRel = 1.0;
  double sRel;
  if (mySol && mySol->getVectorSecSol() && mySol->getScalarSol())
    sRel = fNorm[SimraNorm::L2_SIGMA];
  else
    sRel = hypot(fNorm[SimraNorm::L2_SIGMAh], rNorm[SimraNorm::L2_SIGMA_REC]);

  IFEM::cout << "\n  H1 norm |u^*-u^h|"
             << utl::adjustRight(w-19,"") << rNorm[SimraNorm::H1_Ur_Uh]
             << "\n  L2 norm |p^*-p^h|"
             << utl::adjustRight(w-19,"") << rNorm[SimraNorm::L2_Pr_Ph]
             << "\n  L2 norm |s^*-s^h|"
             << utl::adjustRight(w-19,"") << rNorm[SimraNorm::L2_SIGMA_REC];
  if (mySol) {
    IFEM::cout << "\n  H1 norm |u^*-u|"
               << utl::adjustRight(w-17,"") << rNorm[SimraNorm::H1_Ur_U]
               << "\n  L2 norm |p^*-p|"
               << utl::adjustRight(w-17,"") << rNorm[SimraNorm::L2_Pr_P]
               << "\n  L2 norm |s^*-s|"
               << utl::adjustRight(w-17,"") << rNorm[SimraNorm::L2_SIGMA_E_REC];
  }
  IFEM::cout << "\n  Relative velocity error (%)" << utl::adjustRight(w-29,"")
             << rNorm[SimraNorm::H1_Ur_Uh] / uRel * 100.0
             << "\n  Relative pressure error (%)" << utl::adjustRight(w-29,"")
             << rNorm[SimraNorm::L2_Pr_Ph] / pRel * 100.0
             << "\n  Relative stress error (%)" << utl::adjustRight(w-27,"")
             << rNorm[SimraNorm::L2_SIGMA_REC] / sRel * 100.0
             << "\n  Relative total error (%)" << utl::adjustRight(w-26,"")
             << hypot(rNorm[SimraNorm::H1_Ur_Uh], rNorm[SimraNorm::L2_Pr_Ph]) /
                hypot(uRel, pRel);

  if (mySol) {
    IFEM::cout << "\n  Effectivity index eta^u"
               << utl::adjustRight(w-25,"")
               << rNorm[SimraNorm::H1_Ur_Uh] / fNorm[SimraNorm::H1_E_U]
               << "\n  Effectivity index eta^p"
               << utl::adjustRight(w-25,"")
               << rNorm[SimraNorm::L2_Pr_Ph] / fNorm[SimraNorm::L2_E_P]
               << "\n  Effectivity index eta^s"
               << utl::adjustRight(w-25,"")
               << rNorm[SimraNorm::L2_SIGMA_REC] / fNorm[SimraNorm::L2_E_SIGMA]
               << "\n  Effectivity index eta^total"
               << utl::adjustRight(w-29,"")
               << hypot(rNorm[SimraNorm::H1_Ur_Uh], rNorm[SimraNorm::L2_Pr_Ph]) /
                  hypot(fNorm[SimraNorm::H1_E_U], fNorm[SimraNorm::L2_E_P]);
  }
  IFEM::cout << std::endl;
}


void SIMSimraProject::registerFields(DataExporter& exporter, const Vector& sol,
                                     const Vectors& projs, const Matrix& eNorm) const
{
  exporter.registerField("u and vtef", "u and vtef",
                         DataExporter::SIM,
                         DataExporter::PRIMARY |
                         DataExporter::SECONDARY |
                         DataExporter::NORMS);
  exporter.setFieldValue("u and vtef",this,&sol,&projs,&eNorm);
  static constexpr const char* names[] =
    {"u_x", "u_y", "u_z", "ps", "tk", "td",
     "vtef", "pt", "pts", "rho", "rhos", "strat"};
  for (int i : {3,4,5,7,8,9,10,11}) {
    exporter.registerField(names[i], names[i],
                           DataExporter::SIM, -DataExporter::PRIMARY, "", 1);
    exporter.setFieldValue(names[i], this, &solution[i]);
  }
  exporter.registerField("p", "p", DataExporter::KNOTSPAN);
  exporter.setFieldValue("p", this, &itg.elmPressure);
}
