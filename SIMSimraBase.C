//==============================================================================
//!
//! \file SIMSimraBase.C
//!
//! \date Dec 4 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for simulation drivers using SIMRA patches.
//!
//==============================================================================

#include "SIMSimraBase.h"

#include "ASMs3DSimra.h"
#include "IFEM.h"
#include "LogStream.h"
#include "SIMdependency.h"

#include <ostream>
#include <vector>


class ASMbase;


SIMSimraBase::SIMSimraBase (size_t nf) :
  SIM3D(nf)
{
}


ASMbase* SIMSimraBase::readPatch (std::istream &isp,
                                  int pchInd,
                                  const CharVec &unf,
                                  const char *whiteSpace) const
{
  ASMs3DSimra* pch = new ASMs3DSimra(nf.front());
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
      pch->idx = myModel.size();
    }
  }

  return pch;
}
