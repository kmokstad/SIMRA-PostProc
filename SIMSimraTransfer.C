//==============================================================================
//!
//! \file SIMSimraTransfer.C
//!
//! \date Apr 12 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for transfer of SIMRA results.
//!
//==============================================================================

#include "SIMSimraTransfer.h"
#include "ASMs3DSimra.h"
#include "SimraIntegrand.h"
#include "SimraIO.h"

#include "ASMbase.h"
#include "IFEM.h"
#include "LogStream.h"
#include "Utilities.h"
#include "Vec3.h"
#include "Vec3Oper.h"

#include <algorithm>
#include <cstddef>
#include <ext/alloc_traits.h>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <strings.h>
#include <tinyxml2.h>


SIMSimraTransfer::SIMSimraTransfer (const std::string& context) :
  SIMSimraProject(context)
{
  for (std::array<int,3>& it : configuredRegion)
    std::fill(it.begin(), it.end(), -1);
}


bool SIMSimraTransfer::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),inputContext.c_str())) {
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(), "boundaryfile")) {
        boundaryFile = utl::getValue(child, "boundaryfile");
        std::string type;
        utl::getAttribute(child, "type", type);
        if (type == "parallel")
          bType = PARALLEL;
        else
          bType = SERIAL;
        IFEM::cout << "\tBoundary condition file (type="
                    << (bType == SERIAL?"serial":"parallel")
                    << "): " << boundaryFile << std::endl;
      } else if (!strcasecmp(child->Value(),"use_boundary_transfer")) {
        useBndTransfer = true;
        IFEM::cout << "\tUsing transfer for boundary values" << std::endl;
      } else if (!strcasecmp(child->Value(), "nesting")) {
        NestedRegion& cr = configuredRegion;
        const char* start = child->Attribute("start");
        if (start) {
          std::vector<int> tmp;
          utl::parseIntegers(tmp, start);
          std::copy(tmp.begin(), tmp.begin()+3, cr[0].begin());
        }
        const char* end = child->Attribute("end");
        if (end) {
          std::vector<int> tmp;
          utl::parseIntegers(tmp, end);
          std::copy(tmp.begin(), tmp.begin()+3, cr[1].begin());
        }
        utl::getAttribute(child,"tol", nestingTolerance);
        if (std::find(cr[0].begin(), cr[0].end(), -1) == cr[0].end() &&
            std::find(cr[1].begin(), cr[1].end(), -1) == cr[1].end())
          IFEM::cout << "\tConfigured nesting region =  ("
                     << cr[0][0] << ',' << cr[0][1] << ',' << cr[0][2] << " -> ("
                     << cr[1][0] << ',' << cr[1][1] << ',' << cr[1][2] << ')' << std::endl;
        if (nestingTolerance != -1.0)
          IFEM::cout << "\tTolerance for detecting nesting region = "
                     << nestingTolerance << std::endl;
      } else if (!strcasecmp(child->Value(), "inflow")) {
          utl::parseIntegers(inflow_faces, utl::getValue(child,"inflow"));
          IFEM::cout << "\tConfigured inflow faces:";
          for (int face : inflow_faces)
            IFEM::cout << " " << face;
          IFEM::cout << std::endl;
      }
  }

  return this->SIMSimraProject::parse(elem);
}


bool SIMSimraTransfer::readBoundaryData ()
{
  std::ifstream ifs(boundaryFile);
  return bnd.readFile(ifs, bType == PARALLEL ? this->getNoNodes() : 0);
}


Matrix SIMSimraTransfer::getSolutions () const
{
  Matrix result(solution.size(), this->getNoNodes());
  for (size_t i = 0; i < solution.size(); ++i)
    result.fillRow(i+1, solution[i].data());

  return result;
}


void SIMSimraTransfer::elementTransfer (const SIMSimraTransfer& input)
{
  itg.elmPressure.resize(this->getNoElms());

  if (input.getElmPressures().empty())
    return;

  ASMs3DSimra* opch = static_cast<ASMs3DSimra*>(this->getPatch(1));
  ASMs3DSimra* ipch = static_cast<ASMs3DSimra*>(input.getPatch(1));
  opch->elementTransfer(ipch, input.getElmPressures(), itg.elmPressure);
}


Vector SIMSimraTransfer::nodalTransfer (const SIMSimraTransfer& input)
{
  ASMs3DSimra* opch = static_cast<ASMs3DSimra*>(this->getPatch(1));
  const ASMs3DSimra* ipch = static_cast<const ASMs3DSimra*>(input.getPatch(1));
  const SimraBoundaryIO& ci = input.bnd;

  size_t from_nx, from_ny, from_nz;
  ipch->getNoStructNodes(from_nx,from_ny,from_nz);
  size_t to_nx, to_ny, to_nz;
  opch->getNoStructNodes(to_nx,to_ny,to_nz);

  const std::vector<double>& from_u = ipch->getParams()[0];
  const std::vector<double>& from_v = ipch->getParams()[1];
  const std::vector<double>& to_u = opch->getParams()[0];
  const std::vector<double>& to_v = opch->getParams()[1];

  const std::vector<int>& from_map = input.getNodeMapping();
  Vector result(to_nx*to_ny);
  size_t idx = 1;
  for (size_t j = 0; j < to_ny; ++j) {
    size_t f_j = utl::find_closest(from_v,to_v[j]);
    for (size_t i = 0; i < to_nx; ++i, ++idx) {
      size_t f_i = utl::find_closest(from_u,to_u[i]);
      int n = f_i + f_j*from_nx;
      int simra_node = from_map[n];
      int ze = utl::findIndex(ci.wall, simra_node);
      result(idx) = ci.z0[ze];
    }
  }

  return result;
}


void SIMSimraTransfer::boundaryTransfer (const SIMSimraTransfer& in)
{
  const std::vector<int>& nodeMapIn = in.getNodeMapping();
  const std::vector<int>& nodeMapOut = this->getNodeMapping();

  size_t nx, ny, nz;
  const ASMs3DSimra* inpch = static_cast<const ASMs3DSimra*>(in.getPatch(1));
  inpch->getNoStructNodes(nx,ny,nz);

  size_t nx_out, ny_out, nz_out;
  const ASMs3DSimra* pch = static_cast<const ASMs3DSimra*>(this->getPatch(1));
  pch->getNoStructNodes(nx_out,ny_out,nz_out);

  for (int face : {1,2,3,4}) {
    IntVec outNodes;
    this->getPatch(1)->getBoundaryNodes(face, outNodes);
    IFEM::cout << "\n\t\tFace " << face << ": " << outNodes.size() << " nodes, " << std::flush;

    std::vector<int> numBC(outNodes.size());

#pragma omp parallel for schedule(static)
    for (size_t j = 0; j < nz_out; ++j) {
      for (size_t i = 0; i < (face < 3 ? ny_out : nx_out); ++i) {
        size_t idx = j * (face < 3 ? ny_out : nx_out) + i;
        auto c = pch->getNodeParams(outNodes[idx]);
        int inElm = inpch->findElement(c[0], c[1], c[2]);
        IntVec inNodes;
        in.getElmNodes(inNodes, inElm);
        static const std::vector<std::vector<int>> checkNodes = {{1,3,5,7},
                                                                 {2,4,6,8},
                                                                 {1,2,5,6},
                                                                 {3,4,7,8}};
        for (int nNode : checkNodes[face-1]) {
          int simraNode = nodeMapIn[inNodes[nNode-1]-1];
          if (std::find(in.bnd.u.idx.begin(),
                        in.bnd.u.idx.end(), simraNode) != in.bnd.u.idx.end())
            ++numBC[idx];
        }
      }
    }

    size_t nInFlow = 0;
    for (size_t idx = 0; idx < outNodes.size(); ++idx) {
      if (numBC[idx] > 1) {
        int n = outNodes[idx];
        int fNode = nodeMapOut[n-1];
        bnd.u.idx.push_back(fNode);
        bnd.v.idx.push_back(fNode);
        bnd.w.idx.push_back(fNode);
        bnd.k.idx.push_back(fNode);
        bnd.d.idx.push_back(fNode);
        bnd.e.idx.push_back(fNode);
        bnd.u.val.push_back(solution[U_X](n));
        bnd.v.val.push_back(solution[U_Y](n));
        bnd.w.val.push_back(solution[U_Z](n));
        bnd.k.val.push_back(solution[TK](n));
        bnd.d.val.push_back(solution[TD](n));
        bnd.e.val.push_back(solution[PT](n));
        ++nInFlow;
      }
    }

    IFEM::cout << nInFlow << " marked as inflow nodes.";
  }
}


void SIMSimraTransfer::classifyBoundaryNodes ()
{
  const std::vector<int>& nodeMap = this->getNodeMapping();
  static const std::vector<int> all = {3,4,2,1};
  const std::vector<int>& boundaries = inflow_faces.empty() ? all : inflow_faces;
  for (int face : boundaries) {
    IntVec outNodes;
    this->getPatch(1)->getBoundaryNodes(face, outNodes);
    IFEM::cout << "\n\t\tFace " << face << ": " << outNodes.size() << " nodes, " << std::flush;
    Vec3 normal;
    if (face < 3)
      normal.x = -1.0 + (face-1)*2.0;
    else
      normal.y = -1.0 + (face-3)*2.0;

    size_t nInFlow = 0;
    for (int n : outNodes) {
      Vec3 u;
      u.x = solution[U_X](n);
      u.y = solution[U_Y](n);
      u.z = solution[U_Z](n);
      if (u*normal < 0) {
        int fNode = nodeMap[n-1];
        bnd.u.idx.push_back(fNode);
        bnd.v.idx.push_back(fNode);
        bnd.w.idx.push_back(fNode);
        bnd.k.idx.push_back(fNode);
        bnd.d.idx.push_back(fNode);
        bnd.e.idx.push_back(fNode);
        bnd.u.val.push_back(u.x);
        bnd.v.val.push_back(u.y);
        bnd.w.val.push_back(u.z);
        bnd.k.val.push_back(solution[TK](n));
        bnd.d.val.push_back(solution[TD](n));
        bnd.e.val.push_back(solution[PT](n));
        ++nInFlow;
      }
    }
    IFEM::cout << nInFlow << " marked as inflow nodes.";
  }
}


bool SIMSimraTransfer::writeResults ()
{
  if (resultFile.empty()) {
    std::cerr << "Cannot write results to file " << resultFile << std::endl;
    return false;
  }

  if (rType == HISTORY_FILE) {
    std::cerr << "Can not write results to .hist files" << std::endl;
    return false;
  }

  std::ofstream ofs(resultFile);

  size_t nx, ny, nz;
  static_cast<const ASMs3DSimra*>(this->getPatch(1))->getNoStructNodes(nx,ny,nz);
  auto&& copySolution = [nx,ny,nz,this](auto& cr)
  {
    // init.dat uses physical units - no scaling
    size_t idx = 1;
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i, ++idx) {
          cr.u1[k + i*nz + j*nz*nx] = solution[U_X](idx);
          cr.u2[k + i*nz + j*nz*nx] = solution[U_Y](idx);
          cr.u3[k + i*nz + j*nz*nx] = solution[U_Z](idx);
          cr.ps[k + i*nz + j*nz*nx] = solution[PS](idx);
          cr.tk[k + i*nz + j*nz*nx] = solution[TK](idx);
          cr.td[k + i*nz + j*nz*nx] = solution[TD](idx);
          cr.vtef[k + i*nz + j*nz*nx] = solution[VTEF](idx);
          cr.pt[k + i*nz + j*nz*nx] = solution[PT](idx);
          cr.pts[k + i*nz + j*nz*nx] = solution[PTS](idx);
          cr.rho[k + i*nz + j*nz*nx] = solution[RHO](idx);
          cr.rhos[k + i*nz + j*nz*nx] = solution[RHOS](idx);
          cr.strat[k + i*nz + j*nz*nx] = solution[STRAT](idx);
        }
  };

  auto&& copyElmPressure = [nx,ny,nz](const auto& in,
                                      auto& out)
  {
    size_t idx = 0;
    for (size_t k = 0; k < nz-1; ++k)
      for (size_t j = 0; j < ny-1; ++j)
        for (size_t i = 0; i < nx-1; ++i, ++idx) {
          out[k + i*(nz-1) + j*(nz-1)*(nx-1)] = in[idx];
        }
  };

  if (ASMs3DSimra::useDouble) {
    SimraIO<double> cr(this->getNoNodes());
    std::vector<double> ePress(this->getNoElms());
    copySolution(cr);
    copyElmPressure(itg.elmPressure, ePress);
    if (rType == INIT_FILE)
      cr.writeInit(ofs, ePress);
    else
      cr.writeRestart(ofs);
  } else {
    SimraIO<float> cr(this->getNoNodes());
    std::vector<float> ePress(this->getNoElms());
    copySolution(cr);
    copyElmPressure(itg.elmPressure, ePress);
    if (rType == INIT_FILE)
      cr.writeInit(ofs, ePress);
    else
      cr.writeRestart(ofs);
  }
  ++iStep;

  return true;
}


void SIMSimraTransfer::addAtmosphereBC()
{
  IntVec nodes;
  this->getPatch(1)->getBoundaryNodes(6, nodes);
  const std::vector<int>& nodeMap = this->getNodeMapping();
  for (int n : nodes) {
    int fNode = nodeMap[n-1];
    bnd.w.idx.push_back(fNode);
    bnd.k.idx.push_back(fNode);
    bnd.d.idx.push_back(fNode);
    bnd.w.val.push_back(0.0);
    bnd.k.val.push_back(solution[TK](n));
    bnd.d.val.push_back(solution[TD](n));
  }
}


void SIMSimraTransfer::addTerrainBC(const std::vector<double>& z0)
{
  // no slip on bottom
  this->getPatch(1)->getBoundaryNodes(5, bnd.wall);
  const std::vector<int>& nodeMap = this->getNodeMapping();
  for (int n : bnd.wall) {
    int fNode = nodeMap[n-1];
    bnd.u.idx.push_back(fNode);
    bnd.v.idx.push_back(fNode);
    bnd.w.idx.push_back(fNode);
    bnd.e.idx.push_back(fNode);
    bnd.u.val.push_back(0.0);
    bnd.v.val.push_back(0.0);
    bnd.w.val.push_back(0.0);
    bnd.e.val.push_back(solution[PT](n));
  }

  bnd.z0 = z0;

  size_t nx, ny, nz;
  static_cast<const ASMs3DSimra*>(this->getPatch(1))->getNoStructNodes(nx,ny,nz);

  bnd.log = bnd.wall;
  for (int& n : bnd.log)
    n += nx*ny;

  // log layer
  for (int n : bnd.log) {
    int fNode = nodeMap[n-1];
    bnd.k.idx.push_back(fNode);
    bnd.d.idx.push_back(fNode);
    bnd.k.val.push_back(0.0);
    bnd.d.val.push_back(0.0);
  }

  // translate log nodes
  bnd.log.resize(bnd.wall.size());
  for (size_t i = 0; i < bnd.wall.size(); ++i) {
    bnd.wall[i] = nodeMap[bnd.wall[i]-1];
    bnd.log[i] = nodeMap[bnd.log[i]-1];
  }
}


void SIMSimraTransfer::writeBoundaryData ()
{
  std::ofstream ofs(boundaryFile);
  bnd.writeFile(ofs);
}


const std::vector<int>& SIMSimraTransfer::getNodeMapping () const
{
  if (nodeMap.empty()) {
    size_t nx, ny, nz;
    static_cast<const ASMs3DSimra*>(this->getPatch(1))->getNoStructNodes(nx,ny,nz);

    nodeMap.resize(nx*ny*nz);
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i)
          nodeMap[i + j*nx + k*nx*ny] = 1 + k + i*nz + j*nz*nx;
  }

  return nodeMap;
}


bool SIMSimraTransfer::writeGlvBC (int& nBlock)
{
  bool result = true;

  const auto& nmap = this->getNodeMapping();

  std::vector<int> simra_to_ifem(nmap.size());
  for (size_t i = 0; i < nmap.size(); ++i)
    simra_to_ifem[nmap[i]-1] = i+1;

  int field = 14;
  Vector mask(this->getNoNodes()), val(this->getNoNodes());

  auto&& writeVal =
  [&simra_to_ifem,&mask,&val,&field,&nBlock,this](const SimraBoundaryIO::Data& bnd,
                                                  const std::string& prefix)
  {
    mask.fill(0.0);
    val.fill(0.0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < bnd.idx.size(); ++i) {
      int idx = simra_to_ifem[bnd.idx[i]-1];
      mask(idx) = 1.0;
      val(idx) = bnd.val[i];
    }
    bool result = this->writeGlvS(mask, (prefix+"-mask").c_str(), iStep, nBlock, field++);
    result &= this->writeGlvS(val, (prefix+"-vals").c_str(), iStep, nBlock, field++);

    return result;
  };

  result &= writeVal(bnd.u, "u");
  result &= writeVal(bnd.v, "v");
  result &= writeVal(bnd.w, "w");
  result &= writeVal(bnd.k, "k");
  result &= writeVal(bnd.d, "d");
  result &= writeVal(bnd.e, "e");

  return result;
}


SIMSimraTransfer::NestedRegion
SIMSimraTransfer::determineNestingRegion (const SIMSimraTransfer& out) const
{
    const NestedRegion& cr = configuredRegion;
    if (std::find(cr[0].begin(), cr[0].end(), -1) == cr[0].end() &&
        std::find(cr[1].begin(), cr[1].end(), -1) == cr[1].end())
      return configuredRegion;

  if (nestingTolerance == -1.0) {
    std::cerr << "*** Have to detect nesting region but no tolerance is specified." << std::endl;
    return NestedRegion{{{-1}}};
  }

  const ASMs3DSimra* pch = static_cast<const ASMs3DSimra*>(this->getPatch(1));

  return {pch->findNode(out.getNodeCoord(1), nestingTolerance),
          pch->findNode(out.getNodeCoord(out.getNoNodes()), nestingTolerance)};
}


void SIMSimraTransfer::setParams (const SIMSimraTransfer::NestedRegion& region,
                                  const SIMSimraTransfer& in)
{
  ASMs3DSimra* pch = static_cast<ASMs3DSimra*>(this->getPatch(1));
  const ASMs3DSimra* ipch = static_cast<const ASMs3DSimra*>(in.getPatch(1));
  size_t nx, ny, nz;
  pch->getNoStructNodes(nx,ny,nz);
  size_t nx_in, ny_in, nz_in;
  ipch->getNoStructNodes(nx_in,ny_in,nz_in);

  size_t startInNode = region[0][0] + (region[0][1]-1) * ny_in + (region[0][2]-1)*nx_in*ny_in;
  size_t endInNode = region[1][0] + (region[1][1]-1) * ny_in + (region[1][2]-1)*nx_in*ny_in;

  auto start_param = ipch->getNodeParams(startInNode);
  auto end_param = ipch->getNodeParams(endInNode);

  auto& param = pch->getParams();
  for (size_t i = 0; i < nx; ++i)
    param[0][i] = start_param[0] + param[0][i] * (end_param[0]-start_param[0]);

  for (size_t i = 0; i < ny; ++i)
    param[1][i] = start_param[1] + param[1][i] * (end_param[1]-start_param[1]);

  // TODO: Refinement in z?
  param[2] = ipch->getParams()[2];
}


bool SIMSimraTransfer::hasNestingRegion () const
{
  const NestedRegion& cr = configuredRegion;
  return nestingTolerance != -1.0 ||
         (std::find(cr[0].begin(), cr[0].end(), -1) == cr[0].end() &&
          std::find(cr[1].begin(), cr[1].end(), -1) == cr[1].end());
}
