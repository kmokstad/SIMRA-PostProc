//==============================================================================
//!
//! \file SimraIO.C
//!
//! \date Mar 19 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for reading from and writing to SIMRA files.
//!
//==============================================================================

#include "SimraIO.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>


std::string SimraBoundaryIO::make_real_string(double value)
{
  char buffer [21 + 1];
  std::snprintf (buffer, sizeof buffer, "%10.4E", value);

  if (std::fabs(value) < 1e-7)
    value = 0.0;

  if (value == 0.0) {
    return "  0.0000E+00";
  } else {
    if (std::isnan(value))
      return "NAN";

    if (std::isinf(value)) {
      if (value > 0)
        return "INF";
      else
        return "-INF";
    }

    std::string tmpstr(buffer);
    int exp = value < 0.0 ? std::stoi(tmpstr.substr(8, 4)) : std::stoi(tmpstr.substr(7, 4));
    const bool use_exp_char = (exp >= -100) && (exp < 99);

    if (value < 0.0) {
      if (use_exp_char) {
        tmpstr = " -0." + tmpstr.substr(1, 1) + tmpstr.substr(3, 3) + "E";
      } else {
        tmpstr = " -0." + tmpstr.substr(1, 1) + tmpstr.substr(3, 3);
      }
    } else {
      if (use_exp_char) {
        tmpstr = "  0." + tmpstr.substr(0, 1) + tmpstr.substr(2, 3) + "E";
      } else {
        tmpstr = "  0." + tmpstr.substr(0, 1) + tmpstr.substr(2, 3);
      }
    }

    std::snprintf(buffer, sizeof buffer, "%+03i", exp + 1);
    tmpstr = tmpstr + buffer;
    return tmpstr;
  }
}



bool SimraBoundaryIO::readFile (std::istream & ifs, int nNodes)
{
  char header[32];
  ifs.getline(&header[0],32,'\n');
  std::string h = header;
  if (h != "Boundary conditions") {
    std::cerr << "Not a boundary condition file\n";
    return false;
  }

  int nfixu, nfixv, nfixw, nfixp, nfixe, nfixk, nwelm=0, nlog;
  float z0_const;
  ifs >> nfixu >> nfixv >> nfixw >> nfixp >> nfixe >> nfixk;
  if (nNodes)
    ifs >> nwelm;
  ifs >> nlog >> z0_const >> std::ws;
  z0.resize(nlog);
  for (int i = 0; i < nlog; ++i)
    ifs >> z0[i];

  auto&& readComp = [&ifs,&header](const std::string& h, int num, Data& data)
  {
    ifs >> std::ws;
    ifs.getline(&header[0],32,'\n');
    if (h + " :" != header) {
      std::cerr << "Failed to read " << h << " from boundary condition file\n";
      return false;
    }

    ifs >> std::ws;
    data.idx.resize(num);
    data.val.resize(num);
    for (int i = 0; i < num; ++i)
      ifs >> data.idx[i] >> std::ws >> data.val[i];

    return true;
  };

  if (!readComp("U-component",nfixu,u))
    return false;
  if (!readComp("V-component",nfixv,v))
    return false;
  if (!readComp("W-component",nfixw,w))
    return false;
  if (!readComp("P-points",nfixp,p))
    return false;

  ifs >> std::ws;
  ifs.getline(&header[0],32,'\n');
  h = "Wall/Log-points:";
  if (h != header) {
    std::cerr << "Failed to read " << h << " from boundary condition file\n";
    return false;
  }

  if (nNodes && nwelm) {
    wallE.resize(nwelm);
    for (int i = 0; i < nwelm; ++i)
      ifs >> wallE[i];
    ifs >> std::ws;
  }

  wall.resize(nlog);
  log.resize(nlog);
  for (int i = 0; i < nlog; ++i)
    ifs >> wall[i] >> log[i];
  ifs >> std::ws;

  ifs.getline(&header[0],32,'\n');
  h = header;
  if (h != "(K,e)-values:") {
    std::cerr << "Failed to read (K-E)-values from boundary condition file\n";
    return false;
  }

  k.idx.resize(nfixk);
  for (int i = 0; i < nfixk; ++i)
    ifs >> k.idx[i];
  ifs >> std::ws;
  d.idx = k.idx;

  k.val.resize(nfixk);
  d.val.resize(nfixk);
  for (int i = 0; i < nfixk; ++i)
    ifs >> k.val[i] >> d.val[i];

  if (nNodes) {
    ifs >> std::ws;
    distance.resize(nNodes);
    for (int i = 0; i < nNodes; ++i)
      ifs >> distance[i];
  }

  ifs >> std::ws;
  ifs.getline(&header[0],32,'\n');
  h = header;
  if (h != "Pot.temp :") {
    std::cerr << "Failed to read temperature values from boundary condition file\n";
    return false;
  }

  e.idx.resize(nfixe);
  for (int i = 0; i < nfixe; ++i)
    ifs >> e.idx[i];
  ifs >> std::ws;

  e.val.resize(nfixe);
  for (int i = 0; i < nfixe; ++i)
    ifs >> e.val[i];

  return true;
}


void SimraBoundaryIO::writeFile (std::ostream& ofs)
{
  ofs << "Boundary conditions\n";
  ofs << std::setw(8) << u.idx.size()
      << std::setw(8) << v.idx.size()
      << std::setw(8) << w.idx.size()
      << std::setw(8) << p.idx.size()
      << std::setw(8) << e.idx.size()
      << std::setw(8) << k.idx.size()
      << std::setw(8) << wall.size()
      << make_real_string(z0.front()) << "\n";

  for (size_t i = 0; i < z0.size(); ++i)
    ofs << make_real_string(z0[i])
        << ((i % 6 == 5 && i != z0.size()-1) ? "\n":"");

  ofs << "\nU-component :\n";
  for (size_t i = 0; i < u.idx.size(); ++i)
    ofs << std::setw(8) << u.idx[i] << make_real_string(u.val[i])
        << ((i % 3 == 2 && i != u.idx.size()-1) ? "\n":"");
  ofs << "\nV-component :\n";
  for (size_t i = 0; i < v.idx.size(); ++i)
    ofs << std::setw(8) << v.idx[i] << make_real_string(v.val[i])
        << ((i % 3 == 2 && i != v.idx.size()-1) ? "\n":"");
  ofs << "\nW-component :\n";
  for (size_t i = 0; i < w.idx.size(); ++i)
    ofs << std::setw(8) << w.idx[i] << make_real_string(w.val[i])
        << ((i % 3 == 2 && i != w.idx.size()-1) ? "\n":"");
  ofs << "\nP-points :\n";
  ofs << "\nWall/Log-points:\n";
  size_t num = 0;
  for (size_t i = 0; i < wall.size(); ++i) {
    ofs << std::setw(8) << wall[i] << (num % 9 == 8 ? "\n":"");
    ++num;
    ofs << std::setw(8) << log[i] << (num % 9 == 8 && i != wall.size()-1 ? "\n":"");
    ++num;
  }
  ofs << "\n(K,e)-values:\n";
  for (size_t i = 0; i < k.idx.size(); ++i)
    ofs << std::setw(8) << k.idx[i] << (i % 9 == 8 && i != k.idx.size()-1 ? "\n":"");
  ofs << '\n';
  for (size_t i = 0; i < k.val.size(); ++i)
    ofs << make_real_string(k.val[i]) << make_real_string(d.val[i]) << ((i % 3 == 2 && i != k.val.size()-1) ? "\n":"");
  ofs << "\nPot.temp :\n";
  for (size_t i = 0; i < e.idx.size(); ++i)
    ofs << std::setw(8) << e.idx[i] << ((i % 9 == 8 && i != e.idx.size()-1) ? "\n":"");
  ofs << '\n';
  for (size_t i = 0; i < e.val.size(); ++i)
    ofs << make_real_string(e.val[i]) << ((i % 6 == 5 && i != e.val.size()-1) ? "\n":"");
}
