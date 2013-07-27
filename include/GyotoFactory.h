#ifdef GYOTO_USE_XERCES

/**
 * \file GyotoFactory.h
 * \brief XML I/O
 *
 * The Factory is a place where objects are built.
 *
 */
/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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

#ifndef __GyotoFactory_H_
#define __GyotoFactory_H_

#include "GyotoConfig.h"

/// Xerces internal.
/**
 * For some reason it sometimes need to be set to 0 instead of being
 * undefined.
 */
#ifndef XERCES_INCLUDE_WCHAR_H
#define XERCES_INCLUDE_WCHAR_H 0
#endif

#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <GyotoScenery.h>
#include <GyotoPhoton.h>
#include <GyotoSpectrum.h>
#include <sstream>
#include <string>

namespace Gyoto {
  class Factory;
  class FactoryMessenger;
  namespace Spectrometer {
    class Generic;
    class Uniform;
  }
}

/**
 * \class Gyoto::Factory
 * \brief XML input/output
 *
 * The Factory is responsible from building objects from their XML
 * description, and from saving an XML description of existing
 * objects. Since the Factory doesn't know how to build the variety of
 * objects available in Gyoto and in external plug-ins, the Factory
 * orders Metric, Astrobj and Spectrum objects from registered
 * subcontractors (see SmartPointee::Subcontractor_t). The factory an the
 * various subcontractors communicate through a FactoryMessenger.
 *
 * To read an XML file, you simply create an instance of the Factory
 * with a filename, and get whichever object type you are interested
 * in:
 * \code
 *  Gyoto::Factory * factory = new Gyoto::Factory("some/input/file.xml");
 *  const std::string kind = factory->getKind();
 *  if (kind.compare("Scenery")) Gyoto::throwError("I wan't a Scenery");
 *  Gyoto::SmartPointer<Gyoto::Scenery> scenery = factory -> getScenery();
 *  Gyoto::SmartPointer<Gyoto::Screen>  screen = scenery->getScreen();
 *  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> object = scenery->getAstrobj();
 *  Gyoto::SmartPointer<Gyoto::Metric::Generic> metric = scenery->getMetric();
 *  delete factory; factory=NULL;
 * \endcode or, for a single object and without checking the kind
 * (getKind()) first:
 * \code
 *  Gyoto::SmartPointer<Gyoto::Scenery> scenery = Factory("some/input/file.xml").getScenery();
 * \endcode
 *
 *
 * Writing an object to a file is even easier. Assuming "object" below
 * is a Gyoto::SmartPointer<class> where "class" is one of Scenery,
 * Metric::Generic, Astrobj::Generic, Spectrum::Generic, Screen,
 * Photon or Spectrometer:
 * \code
 *  Gyoto::Factory * factory = new Gyoto::Factory(object);
 *  factory -> write("some/output/file.xml");
 *  delete factory; factory=NULL;
 * \endcode
 *
 * or, for short:
 * \code
 *  Gyoto::Factory(object).write("some/output/file.xml");
 * \endcode
 *
 * You can also directly display the object to stdout:
 * \code
 *  std::cout << Gyoto::Factory(object).format() << std::endl;
 * \endcode
 */
class Gyoto::Factory
{
  friend class Gyoto::FactoryMessenger;

 protected:
  // XERCES MACHINERY
  /// Xerces error handler 
  xercesc::ErrorHandler *reporter_;
  /// The document being read or written
  xercesc::DOMDocument *doc_;
  /// Root element in Factory::doc_
  xercesc::DOMElement *root_;
  /// Xerces parser
  xercesc::XercesDOMParser *parser_;
  /// Xerces resolver
  xercesc::DOMXPathNSResolver* resolver_;
  /// Xerces implementation
  xercesc::DOMImplementation* impl_;

  // Elements which must happen only once in a file
  // but may happen about anywhere
  /// XML element representing the Metric
  xercesc::DOMElement *gg_el_;
  /// XML element representing the Astrobj
  xercesc::DOMElement *obj_el_;
  /// XML element representing the Photon
  xercesc::DOMElement *ph_el_;

  // GYOTO elements
  /// The Scenery read from or written to Factory::doc_
  SmartPointer<Scenery> scenery_;
  /// The Metric read from or written to Factory::doc_
  SmartPointer<Metric::Generic> gg_;
  /// The Screen read from or written to Factory::doc_
  SmartPointer<Screen> screen_; 
  /// The Astrobj read from or written to Factory::doc_
  SmartPointer<Astrobj::Generic> obj_;
  /// The Photon read from or written to Factory::doc_
  SmartPointer<Photon> photon_;
  /// The Spectrometer read from or written to Factory::doc_
  SmartPointer<Spectrometer::Generic> spectro_;

  // Factory stuff
  /// XML file name, if actually reading from or writting to file.
  std::string filename_;
  /// Kind of root element (Scenery, Metric etc.)
  std::string kind_;

 public:
  /// Constructor for reading a file
  Factory(char * filename);

  /// Constructor for saving (or printing) a Scenery
  Factory(SmartPointer<Scenery> sc);
  /// Constructor for saving (or printing) a Metric
  Factory(SmartPointer<Metric::Generic> gg);
  /// Constructor for saving (or printing) an Astrobj
  Factory(SmartPointer<Astrobj::Generic> ao);
  /// Constructor for saving (or printing) a Spectrum
  Factory(SmartPointer<Spectrum::Generic> sp);
  /// Constructor for saving (or printing) a Screen
  Factory(SmartPointer<Screen> screen);
  /// Constructor for saving (or printing) a Photon
  Factory(SmartPointer<Photon> photon);
  /// Constructor for saving (or printing) a Spectrometer
  Factory(SmartPointer<Spectrometer::Generic> Spectrometer);

  /// Destructor
  ~Factory();

 private:
  /// Set Xerces reporter
  void setReporter(xercesc::ErrorHandler*);
  /// Get Factory::root_
  xercesc::DOMElement * getRoot();
  /// Get Factory::doc_
  xercesc::DOMDocument* getDoc();

 public:
  /// Get Factory::kind_
  const std::string getKind();

  /// Find Scenery element, instanciate it and get it.
  /**
   * Scenery must be the root element. getScenery() will call
   * getMetric(), getAstrobj() and getScreen().
   */
  Gyoto::SmartPointer<Gyoto::Scenery> getScenery();

  /// Find Metric element, instanciate it and get it.
  /**
   * Metric may be either the root element or directly within the root
   * element.
   */
  Gyoto::SmartPointer<Gyoto::Metric::Generic>  getMetric();

  /// Find Screen element, instanciate it and get it.
  /**
   * Screen may be either the root element or directly within the root
   * element.
   */
  Gyoto::SmartPointer<Gyoto::Screen>  getScreen();

  /// Find Astrobj element, instanciate it and get it.
  /**
   * Astrobj may be either the root element or directly within the root
   * element.
   */
  Gyoto::SmartPointer<Gyoto::Astrobj::Generic> getAstrobj();

  /// Find Photon element, instanciate it and get it.
  /**
   * Photon may be either the root element or directly within the root
   * element.
   */
  Gyoto::SmartPointer<Gyoto::Photon>  getPhoton();

  /// Find Photon element, instanciate it and get it.
  /**
   * Photon may be either the root element or directly within the root
   * element.
   */
  Gyoto::SmartPointer<Gyoto::Spectrum::Generic>  getSpectrum();

  /// Find Spectrometer element, instanciate it and get it.
  /**
   * Spectrometer may be either the root element or directly within the root
   * element.
   */
  Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>  getSpectrometer();

  // XML OUTPUT
  /// Write constructed XML representation to file
  void write(const char* const fname=0);

  /// Get constructed XML representation as std::string
  std::string format();

  // Setting elements
  /// Set Metric for this document.
  /**
   * If called several times for the same document, the metric
   * SmartPointers must point to the same instance or an error will be
   * thrown using Gyoto::throwError().
   */
  void setMetric(SmartPointer<Metric::Generic> gg, xercesc::DOMElement *el);

  /// Set Astrobj for this document.
  /**
   * If called several times for the same document, the astrobj
   * SmartPointers must point to the same instance or an error will be
   * thrown using Gyoto::throwError().
   */
  void setAstrobj(SmartPointer<Astrobj::Generic> ao, xercesc::DOMElement *el);

  /// Set Screen for this document.
  /**
   * If called several times for the same document, the screen
   * SmartPointers must point to the same instance or an error will be
   * thrown using Gyoto::throwError().
   */
  void setScreen(SmartPointer<Screen> scr, xercesc::DOMElement *el);

  /// Set text content of XML element
  void setContent(std::string content, xercesc::DOMElement *el);

  /// Create new XML element without content
  /**
   * E.g.
   * \code
   * <OpticallyThin/>
   * \endcode
   * \param name XML entity name.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, xercesc::DOMElement *pel);

  /// Create new XML element with double value
  /**
   * E.g.
   * \code
   * <Radius> 2. </Radius>
   * \endcode
   * \param name XML entity name.
   * \param value Entity content.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, double value,
		    xercesc::DOMElement *pel);

  /// Create new XML element with integer value
  /**
   * E.g.
   * \code
   * <IntParameter> 7 </IntParameter>
   * \endcode
   * \param name XML entity name.
   * \param value Entity content.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, int value,
		    xercesc::DOMElement *pel);

  /// Create new XML element with integer value
  /**
   * E.g.
   * \code
   * <IntParameter> 7 </IntParameter>
   * \endcode
   * \param name XML entity name.
   * \param value Entity content.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, unsigned int value,
		    xercesc::DOMElement *pel);

  /// Create new XML element with integer value
  /**
   * E.g.
   * \code
   * <IntParameter> 7 </IntParameter>
   * \endcode
   * \param name XML entity name.
   * \param value Entity content.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, long value,
		    xercesc::DOMElement *pel);

  /// Create new XML element with integer value
  /**
   * E.g.
   * \code
   * <IntParameter> 7 </IntParameter>
   * \endcode
   * \param name XML entity name.
   * \param value Entity content.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, unsigned long value,
		    xercesc::DOMElement *pel);

  /// Create new XML element with string content
  /**
   * E.g.
   * \code
   * <StringParameter> Text </StringParameter>
   * \endcode
   *
   * Any parameter can acually be set this way if total control over
   * Text formatting is wished for.
   * \param name XML entity name.
   * \param value Entity content.
   * \param pel Parent XML element.
   */
  void setParameter(std::string name, std::string value,
		    xercesc::DOMElement*pel);

  /// Create new XML element with array content
  /**
   * E.g.
   * \code
   * <Position> 0. 10. 3.14. 0. </Position>
   * \endcode
   * \param[in] name XML entity name.
   * \param[in] val Entity content.
   * \param[in] nelem Number of elements in val to output.
   * \param[in] pel Parent XML element.
   * \param[out] child If not NULL, set to a new Gyoto::FactoryMessenger
   * pointing to the new element. It must be deleted later.
   */
  void setParameter(std::string name, double val[], size_t nelem,
		    xercesc::DOMElement* pel,
		    FactoryMessenger **child = NULL);

  /// Transform relative path into absolute path.
  /**
   * \param relpath Path relative to XML file.
   * \return Absolute path to same file.
   */
  std::string fullPath(std::string relpath);
};

#endif
#endif
