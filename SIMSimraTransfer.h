//==============================================================================
//!
//! \file SIMSimraTransfer.h
//!
//! \date Apr 12 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for transfer of SIMRA results.
//!
//==============================================================================

#ifndef SIM_SIMRA_TRANSFER_H
#define SIM_SIMRA_TRANSFER_H

#include "SIMSimraProject.h"


/*!
  \brief Simulation driver for transfer of SIMRA results.
*/

class SIMSimraTransfer : public SIMSimraProject {
public:
  //! \brief Default constructor.
  //! \param context The base xml tag to parse
  SIMSimraTransfer(const std::string& context = "simra");

  //! \brief Write results to the result file.
  bool writeResults();

  //! \brief Returns full reordered solution vector.
  Matrix getSolutions() const;

  //! \brief Transfer elemental pressures from another mesh.
  //! \param input The simulator with the data to transfer
  //! \param[out] elmPressures the new elemental pressures
  void elementTransfer(const SIMSimraTransfer& input);

  //! \brief Transfer nodal surface roughness from another mesh.
  //! \param input Input simulator
  Vector nodalTransfer(const SIMSimraTransfer& input);

  //! \brief Transfer boundary conditions from another mesh.
  //! \param in Input simulator
  void boundaryTransfer(const SIMSimraTransfer& in);

  //! \brief Classify boundary nodes by flow direction.
  void classifyBoundaryNodes();

  //! \brief Whether or not to use transfer for boundary data.
  bool useBoundaryTransfer() const { return useBndTransfer; }

  //! \brief Returns the file name of the boundary file.
  const std::string& getBoundaryFile() const { return boundaryFile; }

  //! \brief Reads the boundary conditions.
  bool readBoundaryData();

  //! \brief Write boundary condition data.
  //! \param z0 Surface roughness
  void writeBoundaryData(const std::vector<double>& z0);

  //! \brief Write boundary conditions to VTF file.
  //! \param nBlock Running VTF block counter
  bool writeGlvBC(int& nBlock);

protected:
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement *elem) override;

  //! \brief Enumeration of boundary file types.
  enum BoundaryType {
    SERIAL, //!< Serial format
    PARALLEL //!< Parallel format
  };

  BoundaryType bType = SERIAL; //!< Type for boundary file
  std::string boundaryFile; //!< File with boundary conditions
  SimraBoundaryIO bnd; //!< Boundary data
  bool useBndTransfer = false; //!< Whether or not to use transfer for boundary data

private:
  //! \brief Obtain a reference to the IFEM-to-simra node mapping.
  const std::vector<int>& getNodeMapping() const;

  mutable std::vector<int> nodeMap; //!< Maps nodes numbers from IFEM to simra
};


#endif // SIM_SIMRA_TRANSFER_H
