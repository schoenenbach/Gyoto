/**
 * \file GyotoPatternDiskBB.h
 * \brief A PatternDisk object with black body spectrum and
 *  a power law extension up to some rmax_
 *
 *  A sort of composite of Astrobj::PatternDisk and
 *  Astrobj::ThinDiskPL. Behaves like the latter for r >
 *  PatternDiskBB::rPL_) and mostly like the former for r <
 *  PatternDiskBB::rPL_ except PatternDiskBB::emission_ does not yield
 *  directly I<SUB>&nu;</SUB> but temperature.
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __GyotoPatternDiskBB_H_ 
#define __GyotoPatternDiskBB_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class PatternDiskBB; }
}

//#include <GyotoMetric.h>
#include <GyotoPatternDisk.h>
#include <GyotoBlackBodySpectrum.h>

/**
 * \class Gyoto::Astrobj::PatternDiskBB
 * \brief Geometrically thin disk read from FITS file with black body 
 * spectrum and a power law extension up to some rmax_
 * 
 *   This class describes a disk contained in the z=0 (equatorial)
 *   plane, extending from r=r_ISCO to r=rmax_.  The flux emitted
 *   at radius r and longitude phi at frequency nu is given in a FITS
 *   file.
 *
 *  A sort of composite of Astrobj::PatternDisk and
 *  Astrobj::ThinDiskPL. Behaves like the latter for r >
 *  PatternDiskBB::rPL_) and mostly like the former for r <
 *  PatternDiskBB::rPL_ except PatternDiskBB::emission_ does not yield
 *  directly I<SUB>&nu;</SUB> but temperature.
 */
class Gyoto::Astrobj::PatternDiskBB : public Astrobj::PatternDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::PatternDiskBB>;
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< disk black body
  ///< emission law
 private:
  /**
   * \brief 1 if spectral emission.
   *
   * XML: SpectralEmission
   *
   */
  int SpectralEmission_;
  int PLDisk_; ///< Whether the disk has a power law extension
  double PLSlope_; ///< Mass density power-law slope
  double PLRho_; ///< Mass density reference value
  double rPL_; ///< Mass density reference radius
  //  double rmax_; ///< Maximal extension of the disk !!! see Astrobj::rmax_
  // Constructors - Destructor
  // -------------------------
 public:
  PatternDiskBB(); ///< Standard constructor
  
  PatternDiskBB(const PatternDiskBB& ) ;///< Copy constructor
  virtual PatternDiskBB* clone () const;///< Cloner
  
  virtual ~PatternDiskBB() ;            ///< Destructor
  
  // Accessors
  // ---------
 public:

  int setParameter(std::string name, std::string content, std::string unit);

 public:
  using PatternDisk::emission;
  double emission(double nu_em, double dsem,
			  double c_ph[8], double c_obj[8]) const;

  double const * getVelocity() const ;
  void getVelocity(double const pos[4], double vel[4])  ;

  void setMetric(SmartPointer<Metric::Generic> gg); ///< Insures metric is KerrBL

 public:
#ifdef GYOTO_USE_XERCES
  void fillElement(FactoryMessenger *fmp) const ;
#endif

};

#endif
