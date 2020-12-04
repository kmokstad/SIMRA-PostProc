// $Id$
//==============================================================================
//!
//! \file FSWallDistance.C
//!
//! \date Jan 13 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for the Fares-Schröder
//!        wall distance function. See docs/FSWallDistance.pdf and
//!        http://onlinelibrary.wiley.com/doi/10.1002/fld.348/pdf.
//!
//==============================================================================

#include "FSWallDistance.h"
#include "FiniteElement.h"
#include "ElmMats.h"


FSWallDistance::FSWallDistance (unsigned short int n, double w)
  : IntegrandBase(n), G0(w), sigma(0.2), gamma(1.4)
{
  primsol.resize(1);
}


bool FSWallDistance::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                              const Vec3&) const
{
  Matrix& A = static_cast<ElmMats&>(elmInt).A.front();
  Vector& b = static_cast<ElmMats&>(elmInt).b.front();

  // Gradient in integration point
  Vector dUdX;
  if (!fe.dNdX.multiply(elmInt.vec.front(),dUdX,true)) // dUdX = dNdX^t * {u}
    return false;

  // Value in integration point
  double u = elmInt.vec.front().dot(fe.N);

  for (size_t i=1;i <= fe.N.size();++i) {
    double gradNi_gradu = 0.0, gradu_gradu = 0.0;
    for (size_t s=1;s<=nsd;++s) {
      gradNi_gradu += fe.dNdX(i,s)*dUdX(s);
      gradu_gradu  += dUdX(s)*dUdX(s);
    }

    for (size_t j=1;j <= fe.N.size();++j) {
      double gradNi_gradNj = 0.0, gradNj_gradu = 0.0;
      for (size_t s=1; s <= nsd; ++s) {
        gradNi_gradNj += fe.dNdX(i,s)*fe.dNdX(j,s);
        gradNj_gradu  += fe.dNdX(j,s)*dUdX(s);
      }

      // equation (10)
      double h  = u*gradNi_gradNj + fe.N(j)*gradNi_gradu;
      // equation (11)
      h += (sigma-1.0)*(fe.N(i)*gradNj_gradu+fe.N(j)*gradNi_gradu);
      h += (sigma-1.0)*(fe.N(i)*gradNj_gradu+u*gradNi_gradNj);
      // equation (12)
      h += 4.0*gamma*pow(u,3.0)*fe.N(i)*fe.N(j);

      A(i,j) += h*fe.detJxW;
    }

    double h  = u*gradNi_gradu;
    h += (sigma-1.0)*(u*gradNi_gradu+fe.N(i)*gradu_gradu);
    h += gamma*pow(u,4.0)*fe.N(i);

    b(i) += -h*fe.detJxW;
  }

  return true;
}


std::string FSWallDistance::getField1Name (size_t, const char* prefix) const
{
  return prefix ? prefix + std::string(" g") : std::string("g");
}


std::string FSWallDistance::getField2Name (size_t, const char* prefix) const
{
  return prefix ? prefix + std::string(" y") : std::string("y");
}


bool FSWallDistance::evalSol2 (Vector& s, const Vectors& eV,
                               const FiniteElement& fe, const Vec3&) const
{
  if (eV.empty())
    return false;

  s = std::vector<double>({ 1.0/fe.N.dot(eV.front()) - 1.0/G0 });

  return true;
}
