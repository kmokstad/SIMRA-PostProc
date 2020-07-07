// $Id$
//==============================================================================
//!
//! \file SimraIntegrand.h
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations simra result projection.
//!
//==============================================================================

#ifndef _SIMRA_INTEGRAND_H_
#define _SIMRA_INTEGRAND_H_

#include "IntegrandBase.h"
#include "GlobalIntegral.h"
#include "EqualOrderOperators.h"


/*!
  \brief Class representing the integrand of the Darcy problem.
*/

class SimraIntegrand : public IntegrandBase
{
public:
  //! \brief The constructor initializes all pointers to zero.
  explicit SimraIntegrand();

  //! \brief Empty destructor.
  virtual ~SimraIntegrand() {}

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] elmVec Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol2(Vector& s, const Vectors& elmVec,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group = 0) const override;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField1Name(size_t i, const char* prefix = 0) const override;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix = 0) const override;

  //! \brief Calculates the velocity at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  Vec3 velocity(const FiniteElement& fe, const Vectors& vec) const;

  //! \brief Calculates the velocity gradient at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[out] grad Velocity gradient tensor at current point
  //! \param[in] vec Primary solution vector for current element
  void velocityGradient(const FiniteElement& fe,
                        Tensor& grad, const Vectors& vec) const;

  //! \brief Symmetrize the supplied velocity gradient if necessary.
  //! \param[out] eps Strain tensor at current point
  void strain(Tensor& eps) const;

  //! \brief Calculates the (Cauchy) stress tensor at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[out] sigma Strain tensor at current point
  //! \param[in] vec Primary solution vector for current element
  void stress(const FiniteElement& fe,
              Tensor& sigma, const Vectors& vec) const;

  //! \brief Calculates the pressure at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] vec Primary solution vector for current element
  double pressure(const FiniteElement& fe, const Vectors& vec) const;
};


/*!
  \brief Class representing the integrand of Darcy energy norms.
*/

class SimraNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical heat flux (optional)
  SimraNorm(SimraIntegrand& p);
  //! \brief Empty destructor.
  virtual ~SimraNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;
};

#endif
