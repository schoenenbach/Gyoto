/**
 * \file GyotoMetric.h
 * \brief Base class for metric description
 * 
 * Classes which represent a metric (e.g. Gyoto::Kerr) should inherit
 * from Gyoto::Metric::Generic and implement all of the virtual methods. 
 *
 */

/*
    Copyright 2011, 2013 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoMetric_H_
#define __GyotoMetric_H_ 

#include <iostream>
#include <fstream>
#include <string>

#include <GyotoSmartPointer.h>
#include <GyotoAstrobj.h>
#include <GyotoRegister.h>
#include <GyotoHooks.h>

namespace Gyoto {
  namespace Metric {
    class Generic;

    /// A function to build instances of a specific Metric::Generic sub-class
    /**
     * This is a more specific version of the
     * SmartPointee::Subcontractor_t type. A Metric::Subcontrator_t is
     * called by the Gyoto::Factory to build an instance of the kind
     * of metric specified in an XML file (see Register()). The
     * Factory and Subcontractor_t function communicate through a
     * Gyoto::FactoryMessenger.
     */
    typedef SmartPointer<Metric::Generic> Subcontractor_t(FactoryMessenger*);


    /** 
     * \brief Subcontractor template
     *
     * Instead of reimplementing the wheel, your subcontractor can simply be
     * Gyoto::Metric::Subcontractor<MyKind>
     *
     * \tparam T Sub-class of Metric::Generic 
     */
    template<typename T> SmartPointer<Metric::Generic> Subcontractor
      (FactoryMessenger* fmp) {
      SmartPointer<T> gg = new T();
      gg -> setParameters(fmp);
      return gg;
    }

    /// Query the Metric register
    /**
     * Query the Metric register to get the Metric::Subcontractor_t
     * correspondig to a given kind name. This function is normally
     * called only from the Factory.
     *
     * \param name e.g. "KerrBL"
     * \param errmode int=0. If errmode==0, failure to find a
     *        registered Metric by that name is an error. Else, simply
     *        return NULL pointer in that case.
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Metric::Subcontractor_t* getSubcontractor(std::string name,
						     int errmode=0);

    /// The Metric register
    /**
     * Use the Metric::initRegister() once in your program to
     * initiliaze it, the Metric::Register() function to fill it, and
     * the Metric::getSubcontractor() function to query it.
     */
    extern Register::Entry * Register_;

    /// Make a Metric kind known to the Factory
    /**
     * Register a new Metric::Generic sub-class so that the
     * Gyoto::Factory knows it.
     *
     * \param kind The kind name which identifies this object type in
     * an XML file, as in &lt;Metric kind="name"&gt;
     *
     * \param scp A pointer to the subcontractor, which will
     * communicate whith the Gyoto::Factory to build an instance of
     * the class from its XML description
     */
     void Register(std::string kind, Gyoto::Metric::Subcontractor_t* scp);

     /// Empty the Metric register.
     /**
      *  This must be called once. It is called by
      *  Gyoto::Register::init().
      */
     void initRegister();

  }

  /* Documented elswhere */
  class Worldline;
}

/**
 * \namespace Gyoto::Metric
 * \brief Access to metrics
 * 
 * Objects which describe space-time geometry must inherit from the
 * Gyoto::Metric::Generic class.
 *
 * To be usable, a Metric::Generic sub-class should register a
 * Metric::Subcontractor_t function using the Metric::Register()
 * function. See also \ref writing_plugins_page .
 */
/**
 * \class Gyoto::Metric::Generic
 * \brief Base class for metrics
 *
 * Example: class Gyoto::Metric::KerrBL
 *
 * See Gyoto::Metric for an introduction.
 *
 */
class Gyoto::Metric::Generic
: protected Gyoto::SmartPointee,
  public Gyoto::Hook::Teller
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::Generic>;

 private:
  std::string kind_; ///< Metric kind name (e.g. "KerrBL")
  double mass_;     ///< Mass yielding geometrical unit (in kg).
  int coordkind_; ///< Kind of coordinates (cartesian-like, spherical-like, unspecified)

 public:
  const std::string getKind() const; ///< Get kind_
  void setKind(const std::string); ///< Set kind_
  int getRefCount();
  
  // Constructors - Destructor
  // -------------------------
  //Metric(const Metric& ) ;                ///< Copy constructor
  Generic();
  Generic(const int coordkind); ///< Constructor setting Generic::coordkind_
  Generic(const double mass, const int coordkind);
      ///<  Constructor setting Generic::mass_ and Generic::coordkind_

  virtual ~Generic() ;                        ///< Destructor
  
  // Mutators / assignment
  // ---------------------
  virtual Generic * clone() const ; ///< Virtual copy constructor

  void setMass(const double);        ///< Set mass used in unitLength()
  void setMass(const double, const std::string &unit);        ///< Set mass used in unitLength()

  // Accessors

  int getCoordKind() const; ///< Get coordinate kind
  void setCoordKind(int coordkind); ///< Set coordinate kind

  double getMass() const;        ///< Get mass used in unitLength()
  double getMass(const std::string &unit) const; ///< Get mass used in unitLength()

  /**
   * Metrics implementations are free to express lengths and distances
   * in whatever unit they see fit (presumably most often geometrical
   * units). This function returns this unit in SI (meters).
   */
  double unitLength() const ; ///< M * G / c^2, M is in kg, unitLength in meters
  double unitLength(const std::string &unit) const ; ///< unitLength expressed in specified unit


  virtual void cartesianVelocity(double const coord[8], double vel[3]);
  ///< Compute xprime, yprime and zprime from 8-coordinates

  /**
   * \param coord 4-position (geometrical units);
   * \param v     3-velocity dx1/dx0, dx2/dx0, dx3/dx0;
   * \return tdot = dx0/dtau.
   */
  virtual double SysPrimeToTdot(const double coord[4], const double v[3]) const;
  ///<Compute tdot as a function of dr/dt, dtheta/dt and dphi/dt. Everything is in geometrical units.

  /**
   * \brief Yield circular valocity at a given position (projected on
   * equatorial plane).
   *
   * \param pos input: position,
   * \param vel output: velocity,
   * \param dir 1 for corotating, -1 for counterrotating.
   */
  virtual void circularVelocity(double const pos[4], double vel[4],
				double dir=1.) const ;

  /**
   * Set coord[4] so that the 4-velocity coord[4:7] is lightlike,
   * i.e. of norm 0. There may be up to two solutions. coord[4] is set
   * to the hightest. The lowest can be retrieved using
   * nullifyCoord(double coord[8], double& tdot2) const. Everything is
   * expressed in geometrical units.
   *
   * \param[in,out] coord 8-position, coord[4] will be set according
   * to the other elements;
   */
  virtual void nullifyCoord(double coord[8]) const;
  ///< Set tdot (coord[4]) such that coord is light-like. Everything is in geometrical units.

  /**
   * Set coord[4] so that the 4-velocity coord[4:7] is lightlike,
   * i.e. of norm 0. There may be up to two solutions. coord[4] is set
   * to the hightest. The lowest can be retrieved in tdot2. Everything
   * is expressed in geometrical units.
   *
   * \param[in,out] coord 8-position, coord[4] will be set according
   * to the other elements;
   * \param[out] tdot2    will be set to the smallest solution
   */
  virtual void nullifyCoord(double coord[8], double& tdot2) const;
  ///< Set tdot (coord[4]) such that coord is light-like and return other possible tdot


  /**
   * Compute the scalarproduct of the two quadrivectors u1 and u2 in
   * this Metric, at point pos expressed in coordinate system sys.
   * \param pos 4-position;
   * \param u1 1st quadrivector;
   * \param u2 2nd quadrivector;
   * \return u1*u2
   */
  virtual double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const; ///< Scalar product

  virtual double Norm3D(double* pos) const; ///< not clear
 

  /**
   * \brief Set parameter by name
   *
   * Assume MyKind is a subclass of Metric::Generic which has two
   * members (a string StringMember and a double DoubleMember):
\code
int MyKind::setParameter(std::string name, std::string content, std::string unit) {
 if      (name=="StringMember") setStringMember(content);
 else if (name=="DoubleMember") setDoubleMemeber(atof(content.c_str()), unit);
 else return Generic::setParameter(name, content, unit);
 return 0;
}
\endcode
   * If MyKind is not a direct subclass of Generic, it should call the
   * corresponding setParameter() implementation instead of
   * Generic::setParameter().
   *
   * \param name XML name of the parameter
   * \param content string representation of the value
   * \param unit string representation of the unit
   * \return 0 if this parameter is known, 1 if it is not.
   */
  virtual void setParameter(std::string name,
			    std::string content,
			    std::string unit);

  // Outputs
#ifdef GYOTO_USE_XERCES

  /**
   * \brief Main loop in Subcontractor_t function
   *
   * The Subcontractor_t function for each Metric kind should look
   * somewhat like this (templated as
   * Gyoto::Metric::Subcontractor<MyKind>):
\code
SmartPointer<Metric::Generic>
Gyoto::Metric::MyKind::Subcontractor(FactoryMessenger* fmp) {
  SmartPointer<MyKind> gg = new MyKind();
  gg -> setParameters(fmp);
  return gg;
}
\endcode
   *
   * Each metric kind should implement setParameter(string name,
   * string content, string unit) to interpret the individual XML
   * elements. setParameters() can be overloaded in case the specific
   * Metric class needs low level access to the FactoryMessenger. See
   * Gyoto::Astrobj::UniformSphere::setParameters().
   */
  virtual void setParameters(Gyoto::FactoryMessenger *fmp) ;

  /**
   * Metrics implementations should impement fillElement to save their
   * parameters to XML and call the Metric::fillElement(fmp) for the
   * shared properties
   */

  virtual void fillElement(FactoryMessenger *fmp) ; ///< called from Factory

  /**
   * \brief Process generic XML parameters 
   */
  void processGenericParameters(Gyoto::FactoryMessenger *fmp) ;
#endif

  /**
   * \brief Metric coefficients
   *
   * \param x  4-position at which to compute the coefficient;
   * \param mu 1st index of coefficient, 0&le;&mu;&le;3;
   * \param nu 2nd index of coefficient, 0&le;&nu;&le;3;
   * \return Metric coefficient g<SUB>&mu;,&nu;</SUB> at point x 
   */
  virtual double gmunu(const double * x,
		       int mu, int nu) const
    = 0 ;

  /**
   * \brief Chistoffel symbol
   *
   * Value of Christoffel symbol
   * &Gamma;<SUP>&alpha;</SUP><SUB>&mu;&nu;</SUB> at point
   * (x<SUB>1</SUB>, x<SUB>2</SUB>, x<SUB>3</SUB>).
   */  
  virtual double christoffel(const double coord[8],
			     const int alpha, const int mu, const int nu) const = 0;

  /**
   * \brief RK4 integrator
   */
  virtual int myrk4(Worldline * line, const double coord[8], double h, double res[8]) const;
  
  /**
   * \brief RK4 integrator with adaptive step
   */
  virtual int myrk4_adaptive(Gyoto::Worldline* line, const double coord[8],
			     double lastnorm, double normref,
			     double coordnew[8], double h0, double& h1) const;

  /**
   * \brief Check whether integration should stop
   *
   * The integrating loop will ask this the Metric through this method
   * whether or not it is happy to conitnue the integration.
   * Typically, the Metric should answer 0 when everything is fine, 1
   * when too close to the event horizon, inside the BH...
   *
   * \param coord 8-coordinate vector to check.
   */
  virtual int isStopCondition(double const * const coord) const;

  /**
   * \brief F function such as dy/dtau=F(y,cst)
   */
  virtual int diff(const double y[8], double res[8]) const ;

  /**
   * \brief Set Metric-specific constants of motion. Used e.g. in KerrBL.
   */
  virtual void setParticleProperties(Gyoto::Worldline* line,
				     const double * coord) const;
  

};

#endif
