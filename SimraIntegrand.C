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
#include "ElmMats.h"
#include "ElmNorm.h"
#include "FiniteElement.h"
#include "TensorFunction.h"
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
      grad(i,j) = vec[0].dot(dNdX, i-1, npv);
  }
}


Vec3 SimraIntegrand::velocity(const FiniteElement& fe,
                              const Vectors& vec) const
{
  Vec3 result;
  for (size_t i = 0; i < 3; i++)
    result[i] = vec[0].dot(fe.N, i, npv);

  return result;
}


double SimraIntegrand::temperature(const FiniteElement& fe,
                                   const Vectors& vec) const
{
  return vec[0].dot(fe.N, 4, npv);
}


Vec3 SimraIntegrand::temperatureGradient(const FiniteElement& fe,
                                         const Vectors& vec) const
{
  Vector dT;
  fe.grad(1).multiply(vec[0], dT, 1.0, 0.0, true, npv, 1, 4, 0);
  return Vec3(dT.data(), 3);
}


double SimraIntegrand::pressure(const FiniteElement& fe, const Vectors& vec) const
{
  return elmPressure.empty() ? 0.0 : elmPressure[fe.iel-1];
}


double SimraIntegrand::viscosity(const FiniteElement& fe, const Vectors& vec) const
{
  return vec[0].dot(fe.N, 3, npv);
}


void SimraIntegrand::stress (const FiniteElement& fe, Tensor& sigma,
                             const Vectors& vec) const
{
  this->velocityGradient(fe,sigma,vec);
  this->strain(sigma);
  sigma *= 2.0*this->viscosity(fe, vec);
  sigma -= this->pressure(fe,vec);
}


void SimraIntegrand::strain (Tensor& eps) const
{
  if (0) // formulation & CFD::STRESS)
    eps.symmetrize();
  else
    eps *= 0.5;
}


bool SimraIntegrand::evalSol2(Vector& s, const Vectors& elmVec,
                              const FiniteElement& fe, const Vec3& X) const
{
  Tensor grad(nsd);
  this->velocityGradient(fe, grad, elmVec);

  s.reserve(this->getNoFields(2));
  s = grad;

  s.push_back(this->pressure(fe, elmVec));

  Tensor sigma(nsd);
  this->stress(fe, sigma, elmVec);
  for (size_t i = 0; i < sigma.size(); ++i)
    s.push_back(sigma.ptr()[i]);

  Vec3 dT = this->temperatureGradient(fe, elmVec);
  s.push_back(dT[0]);
  s.push_back(dT[1]);
  s.push_back(dT[2]);

  return true;
}


size_t SimraIntegrand::getNoFields (int fld) const
{
  return fld < 2 ? 5 : nsd*nsd + 1 + nsd*nsd + nsd;
}


std::string SimraIntegrand::getField1Name (size_t i, const char* prefix) const
{
  static const char* s[5] = { "u_x", "u_y", "u_z", "vtef", "pT"};
  if (i == 11)
    return "u_x&&u_y&&u_z&&vtef&&pT";

  return prefix ? prefix + std::string(" ") + s[i] : s[i];
}


std::string SimraIntegrand::getField2Name (size_t i, const char* prefix) const
{
  if (i >= this->getNoFields(2)) return "";

  static const char* s3[] = { "u_x,x", "u_y,x", "u_z,x",
                              "u_x,y", "u_y,y", "u_z,y",
                              "u_x,z", "u_y,z", "u_z,z", "p-rec",
                              "sigma_xx", "sigma_yx", "sigma_zx",
                              "sigma_xy", "sigma_yy", "sigma_zy",
                              "sigma_xz", "sigma_yz", "sigma_zz",
                              "pT,x", "pT,y", "pT,z" };
  std::string n(s3[i]);

  return prefix ? prefix + std::string(" ") + n : n;
}


NormBase* SimraIntegrand::getNormIntegrand (AnaSol* asol) const
{
  return new SimraNorm(*const_cast<SimraIntegrand*>(this), asol);
}


SimraNorm::SimraNorm (SimraIntegrand& p, AnaSol* asol) :
  NormBase(p), aSol(asol)
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

  // Compute temperature
  double pTh = problem.temperature(fe, elmInt.vec);

  // Compute gradient temperature
  Vec3 dpTh = problem.temperatureGradient(fe, elmInt.vec);

  // Viscosity
  double mu = problem.viscosity(fe, elmInt.vec);

  // Stress
  Tensor sigma_h(nsd);
  problem.stress(fe, sigma_h, elmInt.vec);

  // Numerical velocity gradient
  Tensor gradUh(nsd);
  problem.velocityGradient(fe,gradUh,elmInt.vec);
  Tensor epsUh(gradUh);
  problem.strain(epsUh);

  pnorm[L2_Uh] += Uh*Uh*fe.detJxW;
  pnorm[H1_Uh] += mu*epsUh.innerProd(epsUh)*fe.detJxW;
  pnorm[L2_DIV_Uh] += gradUh.trace()*gradUh.trace()*fe.detJxW;
  pnorm[L2_Ph] += Ph*Ph*fe.detJxW;
  pnorm[L2_SIGMAh] += sigma_h.innerProd(sigma_h)*fe.detJxW;
  pnorm[L2_pTh] += pTh*pTh*fe.detJxW;
  pnorm[H1_pTh] += dpTh*dpTh*fe.detJxW;

  // Velocity norms
  if (aSol && aSol->getVectorSol()) {
    // Analytical velocity
    Vec3 U = (*aSol->getVectorSol())(X);
    pnorm[L2_U] += U*U*fe.detJxW;

    U -= Uh;
    pnorm[L2_E_U] += U*U*fe.detJxW;
  }

  // Analytical velocity gradient
  if (aSol && aSol->getVectorSecSol()) {
    Tensor gradU(nsd);
    gradU = (*aSol->getVectorSecSol())(X);
    problem.strain(gradU);
    pnorm[H1_U] += mu*gradU.innerProd(gradU)*fe.detJxW;
    epsUh -= gradU;
    pnorm[H1_E_U] += mu*epsUh.innerProd(epsUh)*fe.detJxW;
    pnorm[TOTAL_ERROR] += mu*epsUh.innerProd(epsUh)*fe.detJxW;
  }

  // Pressure norms
  double P = 0;
  if (aSol && aSol->getScalarSol(0)) {
    const RealFunc* psol = aSol->getScalarSol(0);
    // Analytical pressure
    P = (*psol)(X);

    // L2-norm of analytical pressure
    pnorm[L2_P] += P*P*fe.detJxW;

    // L2-norm of pressure error
    pnorm[L2_E_P] += (P-Ph)*(P-Ph)*fe.detJxW;
    pnorm[TOTAL_ERROR] += (P-Ph)*(P-Ph)*fe.detJxW;
  }

  // Analytical stress
  Tensor sigma(nsd);
  if (aSol && aSol->getVectorSecSol() && aSol->getScalarSol(0)) {
    sigma = (*aSol->getVectorSecSol())(X);
    problem.strain(sigma);
    sigma *= 2.0*mu;
    sigma -= P;
    pnorm[L2_SIGMA] += sigma.innerProd(sigma)*fe.detJxW;
    Tensor sigma_e = sigma;
    sigma_e -= sigma_h;
    pnorm[L2_E_SIGMA] += sigma_e.innerProd(sigma_e)*fe.detJxW;
  }

  // Analytical temperature
  double pT = 0;
  if (aSol && aSol->getScalarSol(1)) {
    pT = (*aSol->getScalarSol(1))(X);
    pnorm[L2_pT] += pT*pT*fe.detJxW;
    pTh -= pT;
    pnorm[L2_E_pT] += pTh*pTh*fe.detJxW;
  }

  // Analytical temperature gradient
  if (aSol && aSol->getScalarSecSol(0)) {
    Vec3 dpT = (*aSol->getScalarSecSol(0))(X);
    pnorm[H1_pT] += dpT*dpT*fe.detJxW;
    dpT -= dpTh;
    pnorm[H1_E_pT] += dpT*dpT*fe.detJxW;
    pnorm[TOTAL_ERROR] += dpT*dpT*fe.detJxW;
  }

  size_t ip = this->getNoFields(1);
  if (!pnorm.psol.empty())
    for (size_t i = 0; i < pnorm.psol.size(); ++i) {
      Tensor gradUs(nsd), sigma_r(nsd);
      double Pr;
      size_t idx = 0;
      for (size_t k = 1; k <= nsd; ++k)
        for (size_t j = 1; j <= nsd; ++j)
          gradUs(j,k) = pnorm.psol[i].dot(fe.N,idx++,problem.getNoFields(2));

      Pr = pnorm.psol[i].dot(fe.N,idx++,problem.getNoFields(2));
      for (size_t j = 1; j <= nsd; ++j)
        for (size_t k = 1; k <= nsd; ++k)
          sigma_r(k,j) = pnorm.psol[i].dot(fe.N,idx++,problem.getNoFields(2));

      Vec3 dpTr;
      for (size_t k = 0; k < 3; ++k)
        dpTr[k] = pnorm.psol[i].dot(fe.N,idx++,problem.getNoFields(2));

      Tensor epssh(gradUs);
      epssh -= gradUh;
      double divus = gradUs.trace();
      problem.strain(epssh);
      Tensor sigma_e = sigma_h;
      sigma_e -= sigma_r;

      pnorm[ip + H1_Ur_Uh] += mu*epssh.innerProd(epssh)*fe.detJxW;
      pnorm[ip + TOTAL_NORM_REC] += (mu*epssh.innerProd(epssh) +
                                    (Pr-Ph)*(Pr-Ph) +
                                    (dpTr-dpTh)*(dpTr-dpTh))*fe.detJxW;

      pnorm[ip + L2_DIV_Ur] += divus*divus*fe.detJxW;
      pnorm[ip + L2_Pr_Ph] += (Pr-Ph)*(Pr-Ph)*fe.detJxW;
      pnorm[ip + L2_SIGMA_REC] += sigma_e.innerProd(sigma_e)*fe.detJxW;
      pnorm[ip + H1_pTr_pTh] += (dpTr-dpTh)*(dpTr-dpTh)*fe.detJxW;

      if (aSol && aSol->getVectorSecSol() && aSol->getScalarSol(0)) {
        Tensor gradU = (*aSol->getVectorSecSol())(X);
        Tensor epsus(gradUs);
        epsus -= gradU;
        problem.strain(epsus);

        const RealFunc* psol = aSol->getScalarSol(0);

        double P = (*psol)(X);
        pnorm[ip+H1_Ur_U] += mu*epsus.innerProd(epsus)*fe.detJxW;
        pnorm[ip+TOTAL_E_REC] += (mu*epsus.innerProd(epsus)+(Pr-P)*(Pr-P))*fe.detJxW;
        pnorm[ip+L2_Pr_P] += (Pr-P)*(Pr-P)*fe.detJxW;

        sigma_r -= sigma;
        pnorm[ip+L2_SIGMA_E_REC] += sigma_r.innerProd(sigma_r)*fe.detJxW;
      }
      if (aSol && aSol->getScalarSecSol(0)) {
        Vec3 dpT = (*aSol->getScalarSecSol(0))(X);
        pnorm[ip+H1_pTr_pT] += (dpTr-dpT)*(dpTr-dpT)*fe.detJxW;
        pnorm[ip+TOTAL_E_REC] += (dpTr-dpT)*(dpTr-dpT)*fe.detJxW;
      }

      ip += this->getNoFields(i+2);
    }

  return true;
}


size_t SimraNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group == 1)
    return 7 + (aSol ? 13 : 0);
  else
    return 7 + (aSol ? 9 : 0);
}


std::string SimraNorm::getName (size_t i, size_t j, const char* prefix) const
{
  static const char* u[20] = {
    "|u^h|_L2",
    "|u^h|_H1",
    "|div u^h|_L2",
    "|p^h|_L2",
    "|s^h|_L2",
    "|pT^h|_L2",
    "|pT^h|_H1",
    "|u|_L2",
    "|e|_L2|, e=u-u^h",
    "|u|_H1",
    "|e|_H1|, e=u-u^h",
    "|p|_L2",
    "|e|_L2|, e=p-p^h",
    "|sigma|_L2",
    "|e|_L2, e=sigma-sigma^h",
    "|pT|_L2",
    "|e|_L2, e=pT-pT^h",
    "|pT|_H1",
    "|e|_H1, e=pT-pT^h",
    "|e_u|_H1 + |e_p|_L2 + |e_pT|_H1"
  };

  static const char* rec[16] = {
    "|u^*-u^h|_H1",
    "|(u^*,p^*,pT^*)-(u^h,p^h,pT^h)|",
    "|p^*-p^h|_L2",
    "|div u^*|_L2",
    "|s^*-s^h|_L2",
    "|pT^*-pT^h|_H1",
    "|u^*-u|_H1",
    "|p^*-p|_L2",
    "|sigma^*-sigma|_L2",
    "|pT^*-pT|_H1",
    "eta^u",
    "eta^p",
    "eta^s",
    "eta^pT",
    "eta^total"
  };

  std::string n;
  if (i < 2) {
    if (j > 20)
      return "";
    n = u[j-1];
  } else {
    if (j > 16)
      return "";
    n = rec[j-1];
  }

  return prefix ? prefix + std::string(" ") + n : n;
}
