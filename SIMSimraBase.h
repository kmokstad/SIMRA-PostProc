//==============================================================================
//!
//! \file SIMSimraBase.h
//!
//! \date Dec 4 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for simulation drivers using SIMRA patches.
//!
//==============================================================================

#ifndef SIM_SIMRA_BASE_H
#define SIM_SIMRA_BASE_H

#include "SIM3D.h"


/*!
  \brief Base class for simulation drivers using SIMRA patches.
*/

class SIMSimraBase : public SIM3D {
public:
  //! \brief Default constructor.
  //! \param nf Number of primary solution fields
  SIMSimraBase(size_t nf);

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] unf Number of unknowns per basis function for each field
  //! \param[in] whiteSpace For message formatting
  ASMbase* readPatch(std::istream& isp, int pchInd, const CharVec& unf,
                     const char* whiteSpace) const override;

};


#endif // SIM_SIMRA_BASE_H
