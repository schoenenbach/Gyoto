/**
 * \file GyotoPowerLawSpectrum.h
 * \brief A power law spectrum : I_nu=constant_*nu^exponent_
 *
 *  Light emitted by an astronomical object
 */

/*
    Copyright 2011, 2013 Thibaut Paumard

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

#ifndef __GyotoPowerLawSpectrum_H_ 
#define __GyotoPowerLawSpectrum_H_ 
#include <GyotoSpectrum.h>

namespace Gyoto {
  namespace Spectrum {
    class PowerLaw;
  }
}


/**
 * \class Gyoto::Spectrum::PowerLaw
 * \brief I_nu=constant_*nu^exponent_
 *
 *  Light emitted by e.g. a Star.
 *
 *  XML stanza:
 *  \code
 *    <Spectrum kind="PowerLaw">
 *      <Exponent> 0. </Exponent>
 *      <Constant> 1. </Constant>
 *    </Spectrum>
 *  \endcode
 */
class Gyoto::Spectrum::PowerLaw : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::PowerLaw>;
 protected:
  double constant_; ///< I_nu=constant_*nu^exponent_
  double exponent_; ///< I_nu=constant_*nu^exponent_

 public:
  PowerLaw();

  /**
   * \brief Constructor setting exponent_ and optionally constant_
   */
  PowerLaw(double exponent, double constant=1.);
  //  PowerLaw(const Spectrum &);
  virtual PowerLaw * clone() const; ///< Cloner

  double getConstant() const; ///< Get constant_
  void setConstant(double); ///< Set constant_
  double getExponent() const; ///< Get exponent_
  void setExponent(double); ///< Set exponent_

  using Gyoto::Spectrum::Generic::operator();
  virtual double operator()(double nu) const;

#ifdef GYOTO_USE_XERCES
  virtual void setParameter(std::string name,
			    std::string content,
			    std::string unit);

  virtual void fillElement(FactoryMessenger *fmp) const ;
#endif
};

#endif
