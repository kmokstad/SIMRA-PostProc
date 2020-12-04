// $Id$
//==============================================================================
//!
//! \file SIMFSWallDistance.h
//!
//! \date Jan 13 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for NURBS-based FEM calculation of the Fares-Schröder
//!        wall distance function on SIMRA patches.
//!
//==============================================================================

#ifndef _SIM_FS_WALL_DISTANCE_H
#define _SIM_FS_WALL_DISTANCE_H

#include "FSWallDistance.h"
#include "NonLinSIM.h"
#include "SIMSimraBase.h"


/*!
  \brief Driver class for NURBS-based FEM calculation of the
  Fares-Schröder wall distance approximation for use in turbulence models.
*/

class SIMFSWallDistance : public SIMSimraBase
{
public:
  //! \brief Default constructor.
  SIMFSWallDistance();

  //! \brief The destructor zero out the integrand pointer (deleted by parent).
  virtual ~SIMFSWallDistance();

  //! \brief Returns the wall value.
  double getWallValue() const { return fsw.getWallValue(); }

  //! \brief Solves the nonlinear wall distance problem.
  bool solve(Vector& distance);

  //! \brief Sets up and solves the nonlinear wall distance problem.
  static bool solveProblem(Vector& distance, const std::string& infile);

protected:
  using SIMSimraBase::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem);

private:
  FSWallDistance fsw;  //!< Fares-Schröder wall distance problem
  NonLinSIM      nSim; //!< The nonlinear solution driver
};

#endif
