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
#include "SimraIO.h"
#include "ASMs3DSimra.h"
#include "IFEM.h"
#include "Utilities.h"
#include "Vec3Oper.h"

#include <tinyxml.h>


SIMSimraTransfer::SIMSimraTransfer (const std::string& context) :
  SIMSimraProject(context)
{
}


bool SIMSimraTransfer::parse (const TiXmlElement *elem)
{
  if (!strcasecmp(elem->Value(),inputContext.c_str())) {
    const TiXmlElement* child = elem->FirstChildElement();
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
        bnd.u.val.push_back(solution[0](n));
        bnd.v.val.push_back(solution[1](n));
        bnd.w.val.push_back(solution[2](n));
        bnd.k.val.push_back(solution[4](n));
        bnd.d.val.push_back(solution[5](n));
        bnd.e.val.push_back(solution[7](n));
        ++nInFlow;
      }
    }

    IFEM::cout << nInFlow << " marked as inflow nodes.";
  }
}


void SIMSimraTransfer::classifyBoundaryNodes ()
{
  const std::vector<int>& nodeMap = this->getNodeMapping();
  for (int face : {1,2,3,4}) {
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
      u.x = solution[0](n);
      u.y = solution[1](n);
      u.z = solution[2](n);
      if (u*normal < 0) {
        int fNode = nodeMap[n-1];
        bnd.u.idx.push_back(fNode);
        bnd.v.idx.push_back(fNode);
        bnd.w.idx.push_back(fNode);
        bnd.k.idx.push_back(fNode);
        bnd.d.idx.push_back(fNode);
        bnd.e.idx.push_back(fNode);
        bnd.u.val.push_back(solution[0](n));
        bnd.v.val.push_back(solution[1](n));
        bnd.w.val.push_back(solution[2](n));
        bnd.k.val.push_back(solution[4](n));
        bnd.d.val.push_back(solution[5](n));
        bnd.e.val.push_back(solution[7](n));
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
    size_t idx = 1;
    for (size_t k = 0; k < nz; ++k)
      for (size_t j = 0; j < ny; ++j)
        for (size_t i = 0; i < nx; ++i, ++idx) {
          cr.u1[k + i*nz + j*nz*nx] = solution[0](idx) / uRef;
          cr.u2[k + i*nz + j*nz*nx] = solution[1](idx) / uRef;
          cr.u3[k + i*nz + j*nz*nx] = solution[2](idx) / uRef;
          cr.ps[k + i*nz + j*nz*nx] = solution[3](idx) / (uRef*uRef);
          cr.tk[k + i*nz + j*nz*nx] = solution[4](idx) / (uRef*uRef);
          cr.td[k + i*nz + j*nz*nx] = solution[5](idx) * lRef / pow(uRef,3.0);
          cr.vtef[k + i*nz + j*nz*nx] = solution[6](idx) / (uRef*lRef);
          cr.pt[k + i*nz + j*nz*nx] = solution[7](idx);
          cr.pts[k + i*nz + j*nz*nx] = solution[8](idx);
          cr.rho[k + i*nz + j*nz*nx] = solution[9](idx);
          cr.rhos[k + i*nz + j*nz*nx] = solution[10](idx);
          cr.strat[k + i*nz + j*nz*nx] = solution[11](idx);
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

  return true;
}


void SIMSimraTransfer::writeBoundaryData (const std::vector<double>& z0)
{
  size_t nx, ny, nz;
  static_cast<const ASMs3DSimra*>(this->getPatch(1))->getNoStructNodes(nx,ny,nz);

  const std::vector<int>& nodeMap = this->getNodeMapping();

  // no slip on bottom
  this->getPatch(1)->getBoundaryNodes(5, bnd.wall);
  for (int n : bnd.wall) {
    int fNode = nodeMap[n-1];
    bnd.u.idx.push_back(fNode);
    bnd.v.idx.push_back(fNode);
    bnd.w.idx.push_back(fNode);
    bnd.e.idx.push_back(fNode);
    bnd.u.val.push_back(0.0);
    bnd.v.val.push_back(0.0);
    bnd.w.val.push_back(0.0);
    bnd.e.val.push_back(solution[7](n));
  }

  bnd.z0 = z0;

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

  // log nodes
  bnd.log.resize(bnd.wall.size());
  for (size_t i = 0; i < bnd.wall.size(); ++i) {
    bnd.wall[i] = nodeMap[bnd.wall[i]-1];
    bnd.log[i] = nodeMap[bnd.log[i]-1];
  }

  // slip on top
  IntVec nodes;
  this->getPatch(1)->getBoundaryNodes(6, nodes);
  for (int n : nodes) {
    int fNode = nodeMap[n-1];
    bnd.w.idx.push_back(fNode);
    bnd.k.idx.push_back(fNode);
    bnd.d.idx.push_back(fNode);
    bnd.w.val.push_back(0.0);
    bnd.k.val.push_back(solution[4](n));
    bnd.d.val.push_back(solution[4](n));
  }

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
    bool result = this->writeGlvS(mask, (prefix+"-mask").c_str(), 1, nBlock, field++);
    result &= this->writeGlvS(val, (prefix+"-vals").c_str(), 1, nBlock, field++);

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
