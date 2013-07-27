/**
 * \file GyotoBlackBodySpectrum.h
 * \brief I_nu(nu, T) = cst_*2*h*nu^3/c^2/(exp(h*nu/k*T)-1.);
 *
 * h = 6.62606896e-34 J.s; J.s = kg.m^2/s
 *
 * k = 1.3806504e-23 J.K-1;   h/k : s*K = K/Hz.
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

#ifndef __GyotoBlackBodySpectrum_H_ 
#define __GyotoBlackBodySpectrum_H_ 
#include "GyotoSpectrum.h"

namespace Gyoto {
  namespace Spectrum {
    class BlackBody;
  }
}

/**
 * \class Gyoto::Spectrum::BlackBody
 * \brief Black Body
 *
 *  Light emitted by e.g. a Star.
 *
 *  Example XML entity:
 *  \code
 *   <Spectrum kind="BlackBody">
 *     <Temperature> 6000 </Temperature>
 *     <Scaling> 1. </Scaling>
 *   </Spectrum>
 *  \endcode
 *
 */
class Gyoto::Spectrum::BlackBody : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::BlackBody>;
 protected:
  double T_; ///< Temperature (K)
  double cst_; ///< Scaling constant
  double Tm1_; ///< 1./T_;

 public:
  BlackBody();

  /**
   * \brief Constructor setting T_ and cst_
   */
  BlackBody(double T, double scaling=1.);
  virtual BlackBody * clone() const; ///< Cloner

  double getTemperature() const; ///< Get constant
  void setTemperature(double); ///< Set constant
  double getScaling() const; ///< Get exponent
  void setScaling(double); ///< Set exponent

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
