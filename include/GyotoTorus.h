/**
 * \file GyotoTorus.h
 * \brief A simple torus
 *
 */

/*
    Copyright 2011 Thibaut Paumard

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

#ifndef __GyotoTorus_H_ 
#define __GyotoTorus_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

namespace Gyoto{
  namespace Astrobj { class Torus; }
}

#include <GyotoStandardAstrobj.h>
#include <GyotoSpectrum.h>
#include <GyotoUtils.h>

/**
 * \class Gyoto::Astrobj::Torus
 * \brief Optically thin or thick torus in circular rotation
 *
 * Any Metric::Generic is acceptable as long as it implements
 * Metric::Generic::circularVelocity().
 */
class Gyoto::Astrobj::Torus : public Gyoto::Astrobj::Standard {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Torus>;
  
  
  // Data : 
  // -----
 protected:
  /**
   * Distance from the center of the coordinate system to the center
   * of the torus tube. The (square of the) radius of a vertical
   * cross-section is stored in critical_value_.
   */
  double c_; ///< Large Radius

  SmartPointer<Spectrum::Generic> spectrum_; ///< Emission law
  SmartPointer<Spectrum::Generic> opacity_; ///< Absorption law

  // Constructors - Destructor
  // -------------------------
 public:
  /**
   *  kind_ =  "Torus", c_ = 3.5, a_=0.5
   */
  Torus(); ///< Default constructor.

  Torus(const Torus& ) ; ///< Copy constructor.
  virtual Torus* clone() const; ///< "Virtual" copy constructor
  
  virtual ~Torus() ; ///< Destructor: does nothing.

  // Accessors
  // ---------
 public:
  /**
   * Get large radius Torus::c_ in geometrical units
   */
  double getLargeRadius() const;

  /**
   * Get large radius Torus::c_ in specified unit
   */
  double getLargeRadius(std::string unit) const;

  /**
   * Get small radius in geometrical units
   */
  double getSmallRadius() const;

  /**
   * Get small radius in specified unit
   */
  double getSmallRadius(std::string unit) const;

  /**
   * \brief Set large radius Torus::c_
   */
  void setLargeRadius(double c);

  /**
   * \brief Set small radius
   */
  void setSmallRadius(double a);

  /**
   * \brief Set large radius Torus::c_ in specified unit
   */
  void setLargeRadius(double c, std::string unit);

  /**
   * \brief Set small radius in specified unit
   */
  void setSmallRadius(double a, std::string unit);

  /**
   * \brief Set Torus::spectrum_
   */
  virtual void setSpectrum(SmartPointer<Spectrum::Generic>);

  /**
   * \brief Get Torus::spectrum_
   */
  virtual SmartPointer<Spectrum::Generic> getSpectrum() const;

  /**
   * \brief Set Torus::opacity_
   */
  virtual void setOpacity(SmartPointer<Spectrum::Generic>);

  /**
   * \brief Get Torus::opacity_
   */
  virtual SmartPointer<Spectrum::Generic> getOpacity() const;

  using Standard::getRmax;
  virtual double getRmax();

  //XML I/O
 public:
  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit) ;

#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ;
  virtual void setParameters(FactoryMessenger *fmp) ;
#endif
  
  // Outputs
  // -------
 public:
  virtual double operator()(double const coord[4]) ;
  
 protected:
  virtual void getVelocity(double const pos[4], double vel[4]) ;

  using Standard::emission;
  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL) const ;
  using Standard::integrateEmission;
  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const;

  virtual double transmission(double nuem, double dsem, double coord[8]) const ;
  
};

#endif
