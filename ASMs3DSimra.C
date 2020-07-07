// $Id$
//==============================================================================
//!
//! \file ASMs3DSimra.C
//!
//! \date Jul 6 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for assembly of SIMRA models.
//!
//==============================================================================

#include "ASMs3DSimra.h"
#include "ElementBlock.h"

#include <numeric>

bool ASMs3DSimra::useDouble = false; //!< Use double precision


ASMs3DSimra::ASMs3DSimra ()
  : ASMs3DLag(4)
{
}


bool ASMs3DSimra::read(std::istream &is)
{
  if (is.tellg() != 0)
    return false;

  auto&& copyMesh = [&](const auto& mesh)
  {
    nel = mesh.nelem;
    nx = mesh.imax;
    ny = mesh.jmax;
    nz = mesh.kmax;
    myMNPC.resize(nel);

    p1 = p2 = p3 = 2;

    size_t idx = 0;
    for (const auto& node : mesh.coords) {
      Vec3 c;
      c.x = node[0]; c.y = node[1]; c.z = node[2];
      this->setCoord(++idx, c);
    }

    idx = 0;
    for (const std::array<int,8>& mnpc : mesh.elms) {
      myMNPC[idx++] = {mnpc.begin(), mnpc.end()};
      for (int& i : myMNPC[idx-1])
        --i;
      std::swap(myMNPC[idx-1][2], myMNPC[idx-1][3]);
      std::swap(myMNPC[idx-1][6], myMNPC[idx-1][7]);
    }
  };

  if (useDouble) {
    Mesh<double> mesh;
    mesh.read(is);
    copyMesh(mesh);
  } else {
    Mesh<float> mesh;
    mesh.read(is);
    copyMesh(mesh);
  }

  return true;
}


bool ASMs3DSimra::empty() const
{
  return myMNPC.empty();
}


bool ASMs3DSimra::tesselate(ElementBlock& grid, const int* npe) const
{
  grid.unStructResize(nel, nnod);

  for (size_t i = 0; i < nnod; i++)
   grid.setCoor(i,this->getCoord(1+i));

  size_t k = 0;
  for (size_t el = 0; el < nel; ++el) {
    std::vector<int> mnpc = myMNPC[el];
    std::swap(mnpc[2], mnpc[3]);
    std::swap(mnpc[6], mnpc[7]);
    for (size_t i = 0; i < 8; ++i)
      grid.setNode(k++,mnpc[i]);
  }

  return true;
}


bool ASMs3DSimra::generateFEMTopology ()
{
  // Number of elements in each direction
  const int nelx = (nx-1)/(p1-1);
  const int nely = (ny-1)/(p2-1);
  const int nelz = (nz-1)/(p3-1);

  myMLGN.resize(nnod);
  std::iota(myMLGN.begin(), myMLGN.end(), 1);

  // Number of elements in patch
  nel = nelx*nely*nelz;
  // Number of nodes per element
  const int nen = p1*p2*p3;
  // Number of nodes in a xy-surface of an element
  const int ct  = p1*p2;

  // Connectivity array: local --> global node relation
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  int i, j, k, a, b, c, iel = 0;
  for (k = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
        myMLGE[iel] = ++gEl;
        myMNPC[iel].resize(nen);
        // First node in current element
        int corner = (p3-1)*(nx*ny)*k + (p2-1)*nx*j + (p1-1)*i;

        for (c = 0; c < p3; c++)
        {
          int cornod = ct*c;
          myMNPC[iel][cornod] = corner + c*nx*ny;
          for (b = 1; b < p2; b++)
          {
            int facenod = cornod + b*p1;
            myMNPC[iel][facenod] = myMNPC[iel][cornod] + b*nx;
            for (a = 1; a < p1; a++)
            {
              myMNPC[iel][facenod+a] = myMNPC[iel][facenod] + a;
              myMNPC[iel][cornod+a]  = myMNPC[iel][cornod] + a;
            }
          }
        }
      }

  return true;
}


/*
bool ASMs3DSimra::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                      const IntegrandBase& integrand,
                                      bool continuous) const
{
  const size_t nnod = this->getNoProjectionNodes();

  int pmax = p1 > p2 ? p1 : p2;
  if (pmax < p3) pmax = p3;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? this->getNoGaussPt(pmax,true) : p1 - 1;
  const int ng2 = continuous ? ng1 : p2 - 1;
  const int ng3 = continuous ? ng2 : p3 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* zg = GaussQuadrature::getCoord(ng3);
  const double* wg = continuous ? GaussQuadrature::getWeight(ng1) : 0;
  if (!xg || !yg || !zg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  std::array<RealArray,3> gpar;
  gpar[0] = this->getGaussPointParameters(gp,0,ng1,xg);
  gpar[1] = this->getGaussPointParameters(gp,1,ng2,yg);
  gpar[2] = this->getGaussPointParameters(gp,2,ng3,zg);

  // Evaluate the secondary solution at all integration points
  Matrix sField;
  if (!this->evalSolution(sField,integrand,gpar.data()))
  {
    std::cerr <<" *** ASMs3DLag::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar[0].size()*gpar[1].size()*gpar[2].size()
              << std::endl;
    return false;
  }

  double dV = 1.0;
  Vector phi(p1*p2*p3), phi2(g1*g2*g3);
  Matrix dNdu, Xnod, J;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0;
  for (int i3 = 0; i3 < nel3; i3++)
    for (int i2 = 0; i2 < nel2; i2++)
      for (int i1 = 0; i1 < nel1; i1++, iel++)
      {
        if (MLGE[iel] < 1) continue; // zero-volume element

        if (continuous)
        {
          // Set up control point (nodal) coordinates for current element
          if (!this->getElementCoordinates(Xnod,1+iel))
            return false;
          else if ((dV = 0.125*this->getParametricVolume(1+iel)) < 0.0)
            return false; // topology error (probably logic error)
        }

        const IntVec& mnpc = MNPC[iel];

        // --- Integration loop over all Gauss points in each direction --------

        Matrix eA(p1*p2*p3, p1*p2*p3);
        Vectors eB(sField.rows(), Vector(p1*p2*p3));
        for (int k = 0; k < ng3; k++, ip += ng2*(nel2-1)*ng1*nel1)
          for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
            for (int i = 0; i < ng1; i++, ip++)
            {
              if (continuous)
              {
                SplineUtils::extractBasis(spl1[ip],phi,dNdu);
                SplineUtils::extractBasis(spl2[ip],phi2,dNdu);
              }
              else
                phi = spl0[ip].basisValues;

              // Compute the Jacobian inverse and derivatives
              double dJw = dV;
              if (continuous)
              {
                dJw *= wg[i]*wg[j]*wg[k]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
                if (dJw == 0.0) continue; // skip singular points
              }

              // Integrate the mass matrix
              eA.outer_product(phi, phi, true, dJw);

              // Integrate the rhs vector B
              for (size_t r = 1; r <= sField.rows(); r++)
                eB[r-1].add(phi,sField(r,ip+1)*dJw);
            }

        for (int i = 0; i < p1*p2*p3; ++i) {
          for (int j = 0; j < p1*p2*p3; ++j)
            A(mnpc[i]+1, mnpc[j]+1) += eA(i+1, j+1);

          int jp = mnpc[i]+1;
          for (size_t r = 0; r < sField.rows(); r++, jp += nnod)
            B(jp) += eB[r](1+i);
        }
      }

  return true;
}
*/
