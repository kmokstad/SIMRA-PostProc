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
#include "Vec3Oper.h"

#include <numeric>

bool ASMs3DSimra::useDouble = false;


ASMs3DSimra::ASMs3DSimra (size_t nf)
  : ASMs3DLag(nf)
{
}


bool ASMs3DSimra::read (std::istream& is)
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
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i) {
          const auto& node = mesh.coords[k + i*nz + j*nz*nx];
          Vec3 c;
          c.x = node[0]; c.y = node[1]; c.z = node[2];
          this->setCoord(++idx, c);
        }

    idx = 0;
    for (size_t k = 0; k < nz-1; ++k)
      for (size_t j = 0; j < ny-1; ++j)
        for (size_t i = 0; i < nx-1; ++i) {
          const std::array<int,8>& mnpc = mesh.elms[k + i*(nz-1) + j*(nz-1)*(nx-1)];
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

  Vec3 delta = this->getCoord(nx) - this->getCoord(1);
  param[0].resize(nx);
  for (size_t i = 0; i < nx; ++i)
    param[0][i] = (this->getCoord(i+1).x - this->getCoord(1).x) / delta.x;

  delta = this->getCoord((ny-1)*nx+1) - this->getCoord(1);
  param[1].resize(ny);
  for (size_t i = 0; i < ny; ++i)
    param[1][i] = (this->getCoord(i*nx+1).y - this->getCoord(1).y) / delta.y;

  delta = this->getCoord((nz-1)*nx*ny + 1) - this->getCoord(1);
  param[2].resize(nz);
  for (size_t i = 0; i < nz; ++i)
    param[2][i] = (this->getCoord(i*nx*ny+1).z - this->getCoord(1).z) / delta.z;

  return true;
}


bool ASMs3DSimra::empty () const
{
  return myMNPC.empty();
}


bool ASMs3DSimra::tesselate (ElementBlock& grid, const int* npe) const
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


void ASMs3DSimra::getNoStructNodes (size_t& n1, size_t& n2, size_t& n3) const
{
  n1 = nx;
  n2 = ny;
  n3 = nz;
}


void ASMs3DSimra::elementTransfer (const ASMs3DSimra* from,
                                   const Vector& oldValues,
                                   Vector& newValues)
{
  size_t idx = 1;
  for (size_t k = 0; k < nz-1; ++k) {
    double w = (param[2][k+1] + param[2][k]) / 2.0;
    for (size_t j = 0; j < ny-1; ++j) {
      double v = (param[1][j+1] + param[1][j]) / 2.0;
      for (size_t i = 0; i < nx-1; ++i, ++idx) {
        double u = (param[0][i+1] + param[0][i]) / 2.0;
        newValues(idx) = oldValues(from->findElement(u,v,w));
      }
    }
  }
}


std::array<double,3> ASMs3DSimra::getNodeParams(int node) const
{
  int x = (node-1) % nx;
  int y = (node-1) % (nx*ny) / nx;
  int z = (node-1) / (nx*ny);

  return {param[0][x], param[1][y], param[2][z]};
}


int ASMs3DSimra::findElement(double u, double v, double w,
                             double* xi, double* eta, double* zeta) const
{
  int elmx = std::upper_bound(param[0].begin(), param[0].end()-1, u) - param[0].begin() - 1;
  int elmy = std::upper_bound(param[1].begin(), param[1].end()-1, v) - param[1].begin() - 1;
  int elmz = std::upper_bound(param[2].begin(), param[2].end()-1, w) - param[2].begin() - 1;

  if (xi) {
    double du = param[0][elmx+1] - param[0][elmx];
    *xi   = -1.0 + (u - param[0][elmx])*2.0 / du;
  }
  if (eta) {
    double dv = param[1][elmy+1] - param[1][elmy];
    *eta  = -1.0 + (v - param[1][elmy])*2.0 / dv;
  }
  if (zeta) {
    double dw = param[2][elmz+1] - param[2][elmz];
    *zeta = -1.0 + (w - param[2][elmz])*2.0 / dw;
  }

  return 1 + elmx + elmy*(nx-1) + elmz*(ny-1)*(nx-1);
}


bool ASMs3DSimra::evaluate (const ASMbase* basis, const Vector& locVec,
                            RealArray& vec, int basisNum) const
{
  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!basis->evalSolution(sValues,locVec,param.data()))
    return false;
  vec = sValues;
  return true;
}


std::array<int,3> ASMs3DSimra::findNode(const Vec3& nodeCoord, double tol) const
{
  double m = 1e10;
  auto it = std::find_if(coord.begin(), coord.end(),
                         [nodeCoord,&m,tol](const Vec3& X)
                         {
                           m = std::min(m, (nodeCoord-X).length());
                           return (nodeCoord-X).length() < tol;
                         });
  if (it == coord.end()) {
    std::cerr << "\n** Closest node was " << m;
    return {-1,-1,-1};
  }

  int node = it-coord.begin();

  int x = node % nx;
  int y = node % (nx*ny) / nx;
  int z = node / (nx*ny);

  return {x+1,y+1,z+1};
}


void ASMs3DSimra::getBoundaryNodes(int lIndex, IntVec& nodes, int basis, int thick, int, bool local) const
{
  IntVec tmpNodes;
  this->ASMs3D::getBoundaryNodes(lIndex, tmpNodes, basis);
  nodes.resize(tmpNodes.size());

  size_t nx, ny, nz;
  this->getNoStructNodes(nx,ny,nz);

  size_t idx = 0;
  if (lIndex == 3 || lIndex == 4) {
    for (size_t i = 0; i < nx; ++i)
      for (size_t k = 0; k < nz; ++k)
        nodes[idx++] = tmpNodes[k*nx+i];
  } else if (lIndex == 1 || lIndex == 2) {
    for (size_t j = 0; j < ny; ++j)
      for (size_t k = 0; k < nz; ++k)
        nodes[idx++] = tmpNodes[k*ny+j];
  } else if (lIndex == 6) {
    for (size_t i = 0; i < nx; ++i)
      for (size_t j = 0; j < ny; ++j)
        nodes[idx++] = tmpNodes[j*nx+i];
  } else
    nodes = tmpNodes;
}
