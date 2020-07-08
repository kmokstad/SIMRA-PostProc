// $Id$
//==============================================================================
//!
//! \file SimraIntegrand.C
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for simra post-processing.
//!
//==============================================================================

#include "SimraIntegrand.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"


SimraIntegrand::SimraIntegrand () : IntegrandBase(3)
{
  primsol.resize(1);
  npv = 5;
}


void SimraIntegrand::velocityGradient(const FiniteElement& fe,
                                      Tensor& grad, const Vectors& vec) const
{
  for (unsigned short int j = 1; j <= 3; j++)
  {
    Vector dNdX = fe.grad(1).getColumn(j);
    for (unsigned short int i = 1; i <= nsd; i++)
      grad(i,j) = vec[0].dot(dNdX, i-1, 5);
  }
}


Vec3 SimraIntegrand::velocity(const FiniteElement& fe,
                              const Vectors& vec) const
{
  Vec3 result;
  for (size_t i = 0; i < 3; i++)
    result[i] = vec[0].dot(fe.N, i, 5);

  return result;
}


double SimraIntegrand::pressure(const FiniteElement& fe, const Vectors& vec) const
{
  return vec[0].dot(fe.N, 3, 5);
}


double SimraIntegrand::viscosity(const FiniteElement& fe, const Vectors& vec) const
{
  return vec[0].dot(fe.N, 4, 5);
}


void SimraIntegrand::stress (const FiniteElement& fe, Tensor& sigma,
                             const Vectors& vec) const
{
  this->velocityGradient(fe,sigma,vec);
}


void SimraIntegrand::strain (Tensor& eps) const
{
  if (0) // formulation & CFD::STRESS)
    eps.symmetrize();
  else
    eps *= 0.5;
}


bool SimraIntegrand::evalSol2(Vector& s, const Vectors& elmVec,
                              const FiniteElement &fe, const Vec3 &X) const
{
  Tensor grad(nsd);
  this->velocityGradient(fe, grad, elmVec);

  s.reserve(this->getNoFields(2));
  s = grad;

  s.push_back(this->pressure(fe, elmVec));

  return true;
}


size_t SimraIntegrand::getNoFields (int fld) const
{
  return fld < 2 ? 4 : nsd*nsd + 1;
}


std::string SimraIntegrand::getField1Name (size_t i, const char* prefix) const
{
  static const char* s[4] = { "u_x", "u_y", "u_z", "p"};

  return prefix ? prefix + std::string(" ") + s[i] : s[i];
}


std::string SimraIntegrand::getField2Name (size_t i, const char* prefix) const
{
  if (i >= this->getNoFields(2)) return "";

  static const char* s3[] = { "u_x,x", "u_y,x", "u_z,x",
                              "u_x,y", "u_y,y", "u_z,y",
                              "u_x,z", "u_y,z", "u_z,z", "p-rec" };
  std::string n(s3[i]);

  return prefix ? prefix + std::string(" ") + n : n;
}


NormBase* SimraIntegrand::getNormIntegrand (AnaSol* asol) const
{
  return new SimraNorm(*const_cast<SimraIntegrand*>(this));
}


SimraNorm::SimraNorm (SimraIntegrand& p) : NormBase(p)
{
  nrcmp = myProblem.getNoFields(2);
}


bool SimraNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& X) const
{
  SimraIntegrand& problem = static_cast<SimraIntegrand&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  size_t nsd   = fe.dNdX.cols();

  // Computed pressure
  double Ph = problem.pressure(fe, elmInt.vec);

  // Computed velocity
  Vec3 Uh = problem.velocity(fe, elmInt.vec);

  // Viscosity
  double mu = problem.viscosity(fe, elmInt.vec);

  // Numerical velocity gradient
  Tensor gradUh(nsd);
  problem.velocityGradient(fe,gradUh,elmInt.vec);
  Tensor epsUh(gradUh);
  problem.strain(epsUh);

  // L2-norm of the velocity
  pnorm[0] += Uh*Uh*fe.detJxW;

  // |u|_H1
  pnorm[1] += mu*epsUh.innerProd(epsUh)*fe.detJxW;

  // |div u|_L2
  pnorm[2] += gradUh.trace()*gradUh.trace()*fe.detJxW;

  // L2-norm of pressure
  pnorm[3] += Ph*Ph*fe.detJxW;

  size_t ip = this->getNoFields(1);
  if (!pnorm.psol.empty())
    for (size_t i = 0; i < pnorm.psol.size(); ++i) {
      Tensor gradUs(nsd);
      double Ps;
      size_t idx = 0;
      for (size_t k = 1; k <= nsd; ++k)
        for (size_t j = 1; j <= nsd; ++j)
          gradUs(j,k) = pnorm.psol[i].dot(fe.N,idx++,problem.getNoFields(2));

      Ps = pnorm.psol[i].dot(fe.N,idx++,problem.getNoFields(2));

      Tensor epssh(gradUs);
      epssh -= gradUh;
      double divus = gradUs.trace();
      problem.strain(epssh);

      pnorm[ip + 0] += mu*epssh.innerProd(epssh)*fe.detJxW;
      pnorm[ip + 1] += (mu*epssh.innerProd(epssh) +
                                   (Ps-Ph)*(Ps-Ph))*fe.detJxW;

      pnorm[ip + 2] += divus*divus*fe.detJxW;
      pnorm[ip + 3] += (Ps-Ph)*(Ps-Ph)*fe.detJxW;

      ip += this->getNoFields(i+2);
    }

  return true;



  return true;
}


size_t SimraNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return 4;
}


std::string SimraNorm::getName (size_t i, size_t j, const char* prefix) const
{
  static const char* u[4] = {
    "|u^h|_L2",
    "|u^h|_H1",
    "|div u^h|_L2",
    "|p^h|_L2",
  };

  static const char* rec[13] = {
    "|u^*-u^h|_H1",
    "|(u^*,p^*)-(u^h,p^h)|",
    "|div u^*|_L2",
    "|p^*-p^h|_L2",
  };

  std::string n;
  if (i < 2) {
    if (j > 4)
      return "";
    n = u[j-1];
  } else {
    if (j > 4)
      return "";
    n = rec[j-1];
  }

  return prefix ? prefix + std::string(" ") + n : n;
}
