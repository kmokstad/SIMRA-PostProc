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
#include "SimraIntegrand.h"
#include "SIMSimraBase.h"

#include <fstream>


class DataExporter;

/*!
  \brief Simulation driver for projection of SIMRA results.
*/

class SIMSimraProject : public SIMSimraBase {
public:
  //! \brief Default constructor.
  SIMSimraProject();

  //! \brief Empty destructor.
  virtual ~SIMSimraProject() { myInts.clear(); myProblem = nullptr; }

  //! \brief Read results from the result file.
  bool readResults();

  //! \brief Returns solution time.
  double getSolutionTime() const { return solTime; }

  //! \brief Returns primary solution vector.
  Vector getSolution() const;

  //! \brief Write all solution vectors to VTF.
  //! \param nBlock Running VTF block counter
  //! \param sol Primary solution vectors
  bool writeSolutionVectors(int& nBlock, const Vector& sol);

  //! \brief Print solution norms to terminal.
  //! \param gNorm Global norms
  void printSolutionNorms(const Vectors& gNorm) const;

  //! \brief Apply post-processing tasks on norms.
  bool postProcessNorms(Vectors& gNorm, Matrix* eNormp) override;

  //! \brief Register fields for exporter output.
  void registerFields(DataExporter& exporter, const Vector& sol,
                      const Vectors& projs, const Matrix& eNorm) const;

  //! \brief Name for data exporter.
  std::string getName() const override { return "Simra"; }

protected:
  //! \brief Prints a norm group to the log stream.
  void printNormGroup(const Vector& rNorm,
                      const Vector& fNorm,
                      const std::string& name) const;

  //! \brief Print exact solution and error norms.
  void printExactNorms(const Vector& gNorm, size_t w) const;

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement *elem) override;

  //! \brief Enumeration of result file types.
  enum ResultsType {
    RESTART_FILE, //!< Restart file, holds one time step
    HISTORY_FILE  //!< History file, holds multiple time steps
  };

  ResultsType rType = RESTART_FILE; //!< Type for result file
  std::string resultFile; //!< File with result vectors
  Vectors solution; //!< Solution vector
  bool stratified = true; //!< Include temperature
  double solTime; //!< Solution time
  double uRef = 1.0; //!< Reference velocity value
  double lRef = 1.0; //!< Reference length scale

  SimraIntegrand itg; //!< Integrand to use

private:
  int iStep = 0; //!< Current time step to read
  std::ifstream ifs; //!< File stream for reading
  std::ifstream::pos_type fileSize = 0; //!< Size of file

  //! \brief Structure for temporary storage of results.
  template<class T>
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

    //! \brief Reads results from a restart file.
    //! \param ifs File stream to read from
    void readRestart(std::istream& ifs)
    {
      int header;
      ifs.read(reinterpret_cast<char*>(&header), 4);
      ifs.read(reinterpret_cast<char*>(&time), sizeof(T));
      for (size_t i = 0; i < u1.size(); ++i) {
        ifs.read(reinterpret_cast<char*>(&u1[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&u2[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&u3[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&ps[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&tk[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&td[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&pt[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&pts[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&rho[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
      }
      ifs.read(reinterpret_cast<char*>(&header), 4);

      ifs.read(reinterpret_cast<char*>(&header), 4);
      ifs.read(reinterpret_cast<char*>(strat.data()), strat.size()*sizeof(T));
      ifs.read(reinterpret_cast<char*>(&header), 4);
    }

    //! \brief Reads results for a single time step from a history file.
    //! \param ifs File stream to read from
    void readHistory(std::istream& ifs, std::vector<T>& elmPressures)
    {
      int header;
      ifs.read(reinterpret_cast<char*>(&header), 4);
      ifs.read(reinterpret_cast<char*>(&time), sizeof(T));
      for (size_t i = 0; i < u1.size(); ++i) {
        ifs.read(reinterpret_cast<char*>(&u1[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&u2[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&u3[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&ps[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&tk[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&td[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&pt[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&pts[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&rho[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
        ifs.read(reinterpret_cast<char*>(&strat[i]), sizeof(T));
      }
      ifs.read(reinterpret_cast<char*>(&header), 4);

      ifs.read(reinterpret_cast<char*>(&header), 4);
      ifs.read(reinterpret_cast<char*>(elmPressures.data()),elmPressures.size()*sizeof(T));
      ifs.read(reinterpret_cast<char*>(&header), 4);
    }

    T time; //!< Time level for solution
    std::vector<T> u1; //!< X velocities
    std::vector<T> u2; //!< Y velocities
    std::vector<T> u3; //!< Z velocities
    std::vector<T> ps; //!< Pressures
    std::vector<T> tk; //!< Turbulent kinetic energy
    std::vector<T> td; //!< Turbulent dissipation
    std::vector<T> vtef; //!< Effective viscosity
    std::vector<T> pt; //!< Potential temperature
    std::vector<T> pts; //!< Potential temperature with stratification
    std::vector<T> rho; //!< Densities
    std::vector<T> rhos; //!< Stratified densities
    std::vector<T> strat; //!< Stratified temperature
  };
};


#endif // SIMSIMRAPROJECT_H
