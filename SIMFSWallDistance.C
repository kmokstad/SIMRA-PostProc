// $Id$
//==============================================================================
//!
//! \file SIMFSWallDistance.C
//!
//! \date Dec 4 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for NURBS-based FEM calculation of the Fares-Schröder
//!        wall distance function on SIMRA patches.
//!
//==============================================================================

#include "SIMFSWallDistance.h"

#include "IFEM.h"
#include "LogStream.h"
#include "SIMadmin.h"
#include "SIMbase.h"
#include "SIMenums.h"
#include "SIMoptions.h"
#include "Utilities.h"

#include <cstdlib>
#include <ostream>
#include <strings.h>
#include <tinyxml2.h>
#include <vector>


SIMFSWallDistance::SIMFSWallDistance() :
  SIMSimraBase(1),
  fsw(3),
  nSim(*this,NonLinSIM::L2)
{
  myProblem = &fsw;
}


SIMFSWallDistance::~SIMFSWallDistance()
{
  myProblem = nullptr;
  myInts.clear();
 }


bool SIMFSWallDistance::solve(Vector& distance)
{
  Matrix tmp;
  this->initSystem(opt.solver,1,1,false);
  nSim.init(1,RealArray(this->getNoDOFs(),0.001*fsw.getWallValue()));
  if (nSim.solve() != SIM::CONVERGED)
    return false;
  if (!this->project(tmp,nSim.getSolution()))
    return false;

  distance = tmp;
  return true;
}


bool SIMFSWallDistance::solveProblem(Vector& distance, const std::string& infile)
{
  int msgbak = SIMadmin::msgLevel;
  SIMFSWallDistance fs;
    IFEM::cout <<"\nSolving the terrain-distance problem."
               <<"\n======================================"<< std::endl;

  bool ok = fs.read(infile.c_str()) && fs.preprocess() && fs.solve(distance);
  SIMadmin::msgLevel = msgbak;
  return ok;
}


bool SIMFSWallDistance::parse(const tinyxml2::XMLElement* elem)
{
  bool result = true;
  if (!strcasecmp(elem->Value(),"fswalldistance")) {
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; child && result; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"wall_value")) {
        const char* w = utl::getValue(child,"wall_value");
        if (w) {
          fsw.setWallValue(atof(w));
          IFEM::cout <<"\tWall value = " << atof(w) << std::endl;
        }
      } else if (!strcasecmp(child->Value(),"sigma")) {
        const char* s = utl::getValue(child,"sigma");
        if (s) {
          fsw.setSigma(atof(s));
          IFEM::cout << "\tSigma = " << atof(s) << std::endl;
        }
      } else if (!strcasecmp(child->Value(),"nonlinearsolver"))
        nSim.parse(child);
      else
        result = this->SIMSimraBase::parse(child);
  } else
    result = this->SIMSimraBase::parse(elem);

  return result;
}
