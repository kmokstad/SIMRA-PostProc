// $Id$
//==============================================================================
//!
//! \file FSWallDistance.h
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

#ifndef _FS_WALL_DISTANCE_H
#define _FS_WALL_DISTANCE_H

#include "IntegrandBase.h"
#include "MatVec.h"

#include <cstddef>
#include <string>


class FiniteElement;
class LocalIntegral;
class Vec3;


/*!
  \brief Class representing the integrand of the Fares-Schröder
         wall distance problem.
*/

class FSWallDistance : public IntegrandBase
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] w Wall value
  FSWallDistance(unsigned short int n = 3, double w = 1.0);
  //! \brief Empty destructor.
  virtual ~FSWallDistance() {}

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] eV Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3&) const override;

  //! \brief Returns whether this integrand has explicit boundary contributions.
  bool hasBoundaryTerms() const override { return false; }

  //! \brief Returns the number of primary/secondary solution field components.
  size_t getNoFields(int) const override { return 1; }

  //! \brief Returns the name of a primary solution field component.
  //! \param[in] prefix Name prefix for all components
  std::string getField1Name(size_t, const char* prefix) const override;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t, const char* prefix) const override;

  //! \brief Defines the wall value.
  void setWallValue(double wall) { G0 = wall; }
  //! \brief Returns the wall value.
  double getWallValue() const { return G0; }

  //! \brief Set the sigma value.
  void setSigma(double s) { sigma = s; gamma = 1.0 + 2*sigma; }

protected:
  double G0;    //!< The wall value
  double sigma; //!< Integration parameter
  double gamma; //!< Integration parameter
};

#endif
