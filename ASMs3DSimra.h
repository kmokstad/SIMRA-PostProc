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
  static bool useDouble; //!< True to use double precision

  //! \brief Default constructor.
  //! \param nf Number of primary solution fields
  ASMs3DSimra(size_t nf);
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

  //! \brief Returns number of nodes in each direction.
  //! \param n1 Number of nodes in x
  //! \param n2 Number of nodes in y
  //! \param n3 Number of nodes in z
  void getNoStructNodes(size_t& n1, size_t& n2, size_t& n3) const;

  //! \brief Perform transfer of elemental values from another mesh.
  //! \param[in] from The mesh to transfer from
  //! \param[in] oldValues Old elemental values
  //! \param[out] newValues New elemental values
  void elementTransfer(const ASMs3DSimra* from,
                       const Vector& oldValues, Vector& newValues);

  //! \brief Get parameters for a node.
  std::array<double,3> getNodeParams(int node) const;

  //! \brief Find element for parameter, and optionally calculate local coordinates.
  int findElement(double u, double v, double w, double* xi = nullptr,
                  double* eta = nullptr, double* zeta = nullptr) const override;

  //! \brief Evaluates and interpolates a field over a given geometry.
  //! \param[in] basis The basis of the field to evaluate
  //! \param[in] locVec The coefficients of the field to evaluate
  //! \param[out] vec The obtained coefficients after interpolation
  //! \param[in] basisNum The basis to evaluate for (mixed)
  bool evaluate (const ASMbase* basis, const Vector& locVec,
                 RealArray& vec, int basisNum) const override;

  //! \brief Obtain parameters.
  const std::array<std::vector<double>,3>& getParams() const { return param; }

  //! \brief Obtain parameters.
  std::array<std::vector<double>,3>& getParams() { return param; }

  //! \brief Returns node matching coordinates.
  //! \param coord Coordinates of node to find
  //! \param tol Tolerance for accepted node
  std::array<int,3> findNode(const Vec3& coord, double tol) const;


  using ASMs3DLag::getBoundaryNodes;
  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face
  //! \param nodes Array of node numbers
  //! \param[in] basis Which basis to grab nodes for (for mixed methods)
  //! \param[in] thick Thickness of connection
  //! \param[in] local If \e true, return patch-local node numbers
  void getBoundaryNodes(int lIndex, IntVec& nodes,
                        int basis, int thick,
                        int, bool local) const override;

private:
  //! \brief Temporary structure for holding mesh structure.
  template<class T>
  struct Mesh {
    int npoint; //!< Number of points in mesh
    int nelem;  //!< Number of elements in mesh
    int imax;   //!< Number of nodes in x
    int jmax;   //!< Number of nodes in y
    int kmax;   //!< Number of nodes in z
    int nelM;   //!< Number of macro elements
    std::vector<std::array<T,3>> coords; //!< Coordinates of nodes
    std::vector<std::array<int,8>> elms;  //!< Elemental node connectivities

    void read(std::istream &is)
    {
      int header;
      is.read(reinterpret_cast<char*>(&header), 4);
      is.read(reinterpret_cast<char*>(&npoint), 4);
      is.read(reinterpret_cast<char*>(&nelem), 4);
      is.read(reinterpret_cast<char*>(&imax), 4);
      is.read(reinterpret_cast<char*>(&jmax), 4);
      is.read(reinterpret_cast<char*>(&kmax), 4);
      is.read(reinterpret_cast<char*>(&nelM), 4);
      is.read(reinterpret_cast<char*>(&header), 4);

      coords.resize(npoint);
      is.read(reinterpret_cast<char*>(&header), 4);
      for (int i = 0; i < npoint; ++i)
        is.read(reinterpret_cast<char*>(coords[i].data()), 3*sizeof(T));
      is.read(reinterpret_cast<char*>(&header), 4);

      elms.resize(nelem);
      is.read(reinterpret_cast<char*>(&header), 4);
      for (int i = 0; i < nelem; ++i)
        is.read(reinterpret_cast<char*>(elms[i].data()), 8*4);
    }
  };

  std::array<std::vector<double>,3> param;
};

#endif
