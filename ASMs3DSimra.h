// $Id$
//==============================================================================
//!
//! \file ASMs3DSimra.h
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of structured 3D SIMRA models.
//!
//==============================================================================

#ifndef _ASM_S3D_SIMRA_H
#define _ASM_S3D_SIMRA_H

#include "ASMs3DLag.h"


/*!
  \brief Driver for assembly of structured 3D Lagrange FE models.
  \details This class contains methods for reading from a SIMRA mesh.dat
*/

class ASMs3DSimra : public ASMs3DLag
{
public:
  //! \brief Default constructor.
  ASMs3DSimra();
  //! \brief Empty destructor.
  virtual ~ASMs3DSimra() {}

  //! \brief Creates an instance by reading the given input stream.
  bool read(std::istream& is) override;

  //! \brief Checks if this patch is empty.
  bool empty() const override;

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  bool tesselate(ElementBlock& grid, const int* npe) const override;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the nodal coordinate array, as well as global node and element numbers.
  bool generateFEMTopology() override;

private:
  //! \brief Temporary structure for holding mesh structure.
  struct Mesh {
    int npoint; //!< Number of points in mesh
    int nelem;  //!< Number of elements in mesh
    int imax;   //!< Number of nodes in x
    int jmax;   //!< Number of nodes in y
    int kmax;   //!< Number of nodes in z
    int nelM;   //!< Number of macro elements
    std::vector<std::array<float,3>> coords; //!< Coordinates of nodes
    std::vector<std::array<int,8>> elms;  //!< Elemental node connectivities
  };
};

#endif
