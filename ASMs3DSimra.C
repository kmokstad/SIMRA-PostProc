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


ASMs3DSimra::ASMs3DSimra ()
  : ASMs3DLag(12)
{
}


bool ASMs3DSimra::read(std::istream &is)
{
  if (is.tellg() != 0)
    return false;

  Mesh mesh;
  int header;
  is.read(reinterpret_cast<char*>(&header), 4);
  is.read(reinterpret_cast<char*>(&mesh.npoint), 4);
  is.read(reinterpret_cast<char*>(&mesh.nelem), 4);
  is.read(reinterpret_cast<char*>(&mesh.imax), 4);
  is.read(reinterpret_cast<char*>(&mesh.jmax), 4);
  is.read(reinterpret_cast<char*>(&mesh.kmax), 4);
  is.read(reinterpret_cast<char*>(&mesh.nelM), 4);
  is.read(reinterpret_cast<char*>(&header), 4);

  mesh.coords.resize(mesh.npoint);
  is.read(reinterpret_cast<char*>(&header), 4);
  for (int i = 0; i < mesh.npoint; ++i)
    is.read(reinterpret_cast<char*>(mesh.coords[i].data()), 3*sizeof(float));
  is.read(reinterpret_cast<char*>(&header), 4);

  mesh.elms.resize(mesh.nelem);
  is.read(reinterpret_cast<char*>(&header), 4);
  for (int i = 0; i < mesh.nelem; ++i)
    is.read(reinterpret_cast<char*>(mesh.elms[i].data()), 8*4);

  nel = mesh.nelem;
  nx = mesh.imax;
  ny = mesh.jmax;
  nz = mesh.kmax;
  myMNPC.resize(nel);

  p1 = p2 = p3 = 2;

  size_t idx = 0;
  for (const std::array<float,3>& node : mesh.coords) {
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


