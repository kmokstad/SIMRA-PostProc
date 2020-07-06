//==============================================================================
//!
//! \file SIMSimraProject.h
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for projection of SIMRA results.
//!
//==============================================================================

#ifndef SIMSIMRAPROJECT_H
#define SIMSIMRAPROJECT_H

#include "SIM3D.h"


class SIMSimraProject : public SIM3D {
public:
  //! \brief Default constructor.
  SIMSimraProject();

  //! \brief Empty destructor.
  virtual ~SIMSimraProject() {}

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] unf Number of unknowns per basis function for each field
  //! \param[in] whiteSpace For message formatting
  ASMbase* readPatch(std::istream& isp, int pchInd, const CharVec& unf,
                     const char* whiteSpace) const override;

  //! \brief Read results from the result file.
  bool readResults();

  //! \brief Returns solution vector.
  const Vector& getSolution() const { return solution; }

  //! \brief Returns solution time.
  double getSolutionTime() const { return solTime; }

protected:
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement *elem) override;

  std::string resultFile; //!< File with result vectors
  Vector solution; //!< Solution vector
  double solTime; //!< Solution time

private:
  //! \brief Structure for temporary storage of results.
  struct Result {
    //! \brief Constructor resizes vectors.
    //! !param len Length of results
    Result(int len)
    {
      u1.resize(len);
      u2.resize(len);
      u3.resize(len);
      ps.resize(len);
      tk.resize(len);
      td.resize(len);
      vtef.resize(len);
      pt.resize(len);
      pts.resize(len);
      rho.resize(len);
      rhos.resize(len);
      strat.resize(len);
    }

    float time; //!< Time level for solution
    std::vector<float> u1; //!< X velocities
    std::vector<float> u2; //!< Y velocities
    std::vector<float> u3; //!< Z velocities
    std::vector<float> ps; //!< Pressures
    std::vector<float> tk; //!< Turbulent kinetic energy
    std::vector<float> td; //!< Turbulent dissipation
    std::vector<float> vtef; //!< Effective viscosity
    std::vector<float> pt; //!< Potential temperature
    std::vector<float> pts; //!< Potential temperature with stratification
    std::vector<float> rho; //!< Densities
    std::vector<float> rhos; //!< Stratified densities
    std::vector<float> strat; //!< Stratified temperature
  };
};


#endif // SIMSIMRAPROJECT_H
