/**
 * \file GyotoRotStar3_1.h
 * \brief Numerical metric around a rotating star in 3+1 formalism
 * 
 *
 */

/*
    Copyright 2011 Frederic Vincent

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


#ifndef __GyotoRotStar3_1_H_
#define __GyotoRotStar3_1_H_ 

#include <iostream>
#include <fstream>

class Star_rot;

namespace Gyoto {
  namespace Metric { class RotStar3_1; }
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>
#include <GyotoSmartPointer.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

/**
 * \class Gyoto::Metric::RotStar3_1
 * \brief Numerical metric around a rotating star in 3+1 formalism
 */
class Gyoto::Metric::RotStar3_1 : public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::RotStar3_1>;

 private:
  char* filename_; ///< Lorene output file name
  Star_rot * star_; ///< Pointer to underlying Lorene Star_rot instance 
  int integ_kind_;///< 1 if RotStar3_1::myrk4(), 0 if Metric::myrk4()
 
 public:

  RotStar3_1(); ///< Constructor
  RotStar3_1(const RotStar3_1& ) ;                ///< Copy constructor
  virtual ~RotStar3_1() ;        ///< Destructor
  virtual RotStar3_1* clone() const ;
           ///< Cloner (uses RotStar3_1(file, integ_kind))

  void setFileName(char const *); ///< Set filename_
  char const * getFileName() const; ///< Get filename_

  void setIntegKind(int); ///< Set integ_kind_
  int getIntegKind() const ; ///< Get integ_kind_

  using Metric::Generic::myrk4;


  /**
   * \brief RK4 integrator
   *
   * NB: we use the 6-coordinates, here, unlike the inherited function.
   */

  int myrk4(const double coord[6], double h, double res[6]) const;

  
  /**
   * \brief Adaptive RK4 integrator
   *
   * Dispatches between Generic::myrk4_adaptive() and  myrk4_adaptive(const double coor[6], double lastnorm, double normref, double coornew[6], double cst[2], double& tdot_used, double h0, double& h1, double& hused) const depending on RotStar3_1::integ_kind_
   */
  int myrk4_adaptive(Gyoto::Worldline* line, const double coord[8], double lastnorm, double normref, double coordnew[8], double h0, double& h1) const;

  /**
   * \brief RK4 integrator (3+1)
   *
   * NB: we use the 6-coordinates, here, unlike the inherited function.
   */
  int myrk4_adaptive(const double coor[6], double lastnorm, double normref, double coornew[6], double cst[2], double& tdot_used, double h0, double& h1, double& hused) const;

  /**
   * \brief F function such as dy/dtau=F(y,cst)
   */
  int diff(const double coord[8], double res[8]) const ;

  /**
   * \brief Alternate version of diff(const double coord[8], double res[8]) const
   *
   * Using only 6 parameters. Last int is not used: it is only here to
   * distinguish the signature of the two methods. Could have been
   * done choosing another name, too, but hey...
   */
  int diff(const double y[6], double res[6], int) const ;


  /**
   * \brief Tweak coordinates to insure conservation of cst
   */
  void Normalize4v(const double coordin[6], double coordout[6], const double cst[2], double& tdot_used) const;

  double gmunu(const double * x, int mu, int nu) const ;

  double christoffel(const double coord[8], const int alpha, const int mu, 
		     const int nu) const ;

  double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const ;

  virtual void setParameter(std::string, std::string, std::string);
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp); ///< called from Factory
  virtual void setParameters(Gyoto::FactoryMessenger *fmp) ;
#endif

};

#endif
