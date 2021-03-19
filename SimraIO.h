//==============================================================================
//!
//! \file SimraIO.h
//!
//! \date Mar 19 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for reading from and writing to SIMRA files.
//!
//==============================================================================
#ifndef SIMRA_IO_H_
#define SIMRA_IO_H_

#include <istream>
#include <ostream>
#include <vector>

//! \brief Structure for temporary storage of results.
template<class T> struct SimraIO {
  //! \brief Constructor resizes vectors.
  //! !param len Length of results
  SimraIO(int len)
  {
    u1.resize(len);
    u2.resize(len);
    u3.resize(len);
    ps.resize(len);
    tk.resize(len);
    td.resize(len);
    vtef.resize(len);
    pt.resize(len);
    pts.resize(len);
    rho.resize(len);
    rhos.resize(len);
    strat.resize(len);
  }

  //! \brief Reads results from a restart file.
  //! \param ifs File stream to read from
  void readRestart(std::istream& ifs)
  {
    int header;
    ifs.read(reinterpret_cast<char*>(&header), 4);
    ifs.read(reinterpret_cast<char*>(&time), sizeof(T));
    for (size_t i = 0; i < u1.size(); ++i) {
      ifs.read(reinterpret_cast<char*>(&u1[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&u2[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&u3[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&ps[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&tk[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&td[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&pt[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&pts[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&rho[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
    }
    ifs.read(reinterpret_cast<char*>(&header), 4);

    ifs.read(reinterpret_cast<char*>(&header), 4);
    ifs.read(reinterpret_cast<char*>(strat.data()), strat.size()*sizeof(T));
    ifs.read(reinterpret_cast<char*>(&header), 4);
  }

  //! \brief Write results to a restart file.
  //! \param ifs File stream to write to
  void writeRestart(std::ostream& ofs)
  {
    int header = u1.size()*sizeof(T)*11+sizeof(T);
    ofs.write(reinterpret_cast<char*>(&header), 4);
    time = 0;
    ofs.write(reinterpret_cast<char*>(&time), sizeof(T));
    for (size_t i = 0; i < u1.size(); ++i) {
      ofs.write(reinterpret_cast<char*>(&u1[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&u2[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&u3[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&ps[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&tk[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&td[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&pt[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&pts[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&rho[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
    }
    ofs.write(reinterpret_cast<char*>(&header), 4);

    header = strat.size()*sizeof(T);
    ofs.write(reinterpret_cast<char*>(&header), 4);
    ofs.write(reinterpret_cast<char*>(strat.data()), strat.size()*sizeof(T));
    ofs.write(reinterpret_cast<char*>(&header), 4);
  }

  //! \brief Reads results for a single time step from a history file.
  //! \param ifs File stream to read from
  void readHistory(std::istream& ifs, std::vector<T>& elmPressures)
  {
    int header;
    ifs.read(reinterpret_cast<char*>(&header), 4);
    ifs.read(reinterpret_cast<char*>(&time), sizeof(T));
    for (size_t i = 0; i < u1.size(); ++i) {
      ifs.read(reinterpret_cast<char*>(&u1[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&u2[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&u3[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&ps[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&tk[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&td[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&pt[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&pts[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&rho[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&strat[i]), sizeof(T));
    }
    ifs.read(reinterpret_cast<char*>(&header), 4);

    ifs.read(reinterpret_cast<char*>(&header), 4);
    ifs.read(reinterpret_cast<char*>(elmPressures.data()),elmPressures.size()*sizeof(T));
    ifs.read(reinterpret_cast<char*>(&header), 4);
  }

  //! \brief Skip a history step.
  void skipHistory(std::istream& ifs, std::vector<T>& elmPressures)
  {
    ifs.seekg(4, std::ios::cur);
    ifs.seekg(sizeof(T), std::ios::cur);
    ifs.seekg(u1.size()*sizeof(T)*12, std::ios::cur);
    ifs.seekg(4, std::ios::cur);
    ifs.seekg(4, std::ios::cur);
    ifs.seekg(elmPressures.size()*sizeof(T), std::ios::cur);
    ifs.seekg(4, std::ios::cur);
  }

  //! \brief Reads results for a single time step from a init file.
  //! \param ifs File stream to read from
  void readInit(std::istream& ifs, std::vector<T>& elmPressures)
  {
    int header;
    ifs.read(reinterpret_cast<char*>(&header), 4);
    for (size_t i = 0; i < u1.size(); ++i) {
      ifs.read(reinterpret_cast<char*>(&u1[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&u2[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&u3[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&ps[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&tk[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&td[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&pt[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&pts[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&rho[i]), sizeof(T));
      ifs.read(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
    }
    ifs.read(reinterpret_cast<char*>(&header), 4);

    ifs.read(reinterpret_cast<char*>(&header), 4);
    ifs.read(reinterpret_cast<char*>(elmPressures.data()),elmPressures.size()*sizeof(T));
    ifs.read(reinterpret_cast<char*>(&header), 4);
  }

  //! \brief Write results for a initial condition file.
  //! \param ofs File stream to write to
  void writeInit(std::ostream& ofs, const std::vector<T>& elmPressures)
  {
    int header = u1.size()*sizeof(T)*11;
    ofs.write(reinterpret_cast<char*>(&header), 4);
    for (size_t i = 0; i < u1.size(); ++i) {
      ofs.write(reinterpret_cast<char*>(&u1[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&u2[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&u3[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&ps[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&tk[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&td[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&vtef[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&pt[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&pts[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&rho[i]), sizeof(T));
      ofs.write(reinterpret_cast<char*>(&rhos[i]), sizeof(T));
    }
    ofs.write(reinterpret_cast<char*>(&header), 4);
    header = elmPressures.size()*sizeof(T);
    ofs.write(reinterpret_cast<char*>(&header), 4);
    ofs.write(reinterpret_cast<const char*>(elmPressures.data()), elmPressures.size()*sizeof(T));
    ofs.write(reinterpret_cast<char*>(&header), 4);
  }

  T time; //!< Time level for solution
  std::vector<T> u1; //!< X velocities
  std::vector<T> u2; //!< Y velocities
  std::vector<T> u3; //!< Z velocities
  std::vector<T> ps; //!< Pressures
  std::vector<T> tk; //!< Turbulent kinetic energy
  std::vector<T> td; //!< Turbulent dissipation
  std::vector<T> vtef; //!< Effective viscosity
  std::vector<T> pt; //!< Potential temperature
  std::vector<T> pts; //!< Potential temperature with stratification
  std::vector<T> rho; //!< Densities
  std::vector<T> rhos; //!< Stratified densities
  std::vector<T> strat; //!< Stratified temperature
};


//! \brief Structure for storage of boundary data.
struct SimraBoundaryIO {

  //! \brief Print a properly formatted float value.
  std::string make_real_string(double value) const;

  //! \brief Reads a 'boun.dat' file.
  //! \param ifs Stream to read from
  //! \param nNodes If non-null read a parallel-formatted boun.dat
  bool readFile(std::istream& ifs, int nNodes);

  //! \brief Writes a 'boun.dat' file.
  //! \details Only serial format is currently supported
  void writeFile(std::ostream& ofs);

  //! \brief Struct holding info for a field.
  struct Data {
    std::vector<int> idx; //!< Indices for entries
    std::vector<double> val; //!< Values for entries
  };

  Data u; //!< X velocities
  Data v; //!< Y velocities
  Data w; //!< Z velocities
  Data p; //!< Pressures
  Data e; //!< Potential temperature
  Data k; //!< Turbulent kinetic energy
  Data d; //!< Turbulent dissipation
  std::vector<double> z0; //!< z0 values
  std::vector<int> wallE; //!< Wall elements
  std::vector<int> wall; //!< Wall nodes
  std::vector<int> log; //!< Log nodes
  std::vector<double> distance; //!< Wall distance
};

#endif
