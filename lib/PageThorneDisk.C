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

#include "GyotoPhoton.h"
#include "GyotoPageThorneDisk.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>
#include <cstring>
#include <time.h> 

#include <gsl/gsl_integration.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

PageThorneDisk::PageThorneDisk() :
  ThinDisk("PageThorneDisk"), aa_(0.), aa2_(0.),
  x0_(0.), x1_(0.), x2_(0.), x3_(0.), blackbody_(0), mdot_(0),
  uniflux_(0), spectrumBB_(NULL)
{
  if (debug()) cerr << "DEBUG: PageThorneDisk Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

PageThorneDisk::PageThorneDisk(const PageThorneDisk& o) :
  ThinDisk(o), aa_(o.aa_), aa2_(o.aa2_),
  x0_(o.x0_), x1_(o.x1_), x2_(o.x2_), x3_(o.x3_),
  blackbody_(o.blackbody_), uniflux_(o.uniflux_), 
  mdot_(0), spectrumBB_(NULL)
{
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;
  gg_->hook(this);
}
PageThorneDisk* PageThorneDisk::clone() const
{ return new PageThorneDisk(*this); }

PageThorneDisk::~PageThorneDisk() {
  GYOTO_DEBUG<<endl;
  if (gg_) gg_->unhook(this);
}

void PageThorneDisk::updateSpin() {
  if (!gg_) return;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    aa_ = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getSpin();
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    aa_ = static_cast<SmartPointer<Metric::KerrKS> >(gg_) -> getSpin();
    break;
  default:
    throwError("PageThorneDisk::getSpin(): unknown COORDKIND");
  }
  aa2_=aa_*aa_;
  double z1 =1.+pow((1.-aa2_),1./3.)*(pow((1.+ aa_),1./3.)+pow((1.-aa_),1./3.));
  double z2 = pow(3.*aa2_ + z1*z1,0.5);
  double acosaao3= acos(aa_)/3.;

  x0_ = sqrt((3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),0.5)));
  x1_ = 2.*cos(acosaao3 - M_PI/3.);
  x2_ = 2.*cos(acosaao3 + M_PI/3.); 
  x3_ = -2.*cos(acosaao3);
  if (rin_==0.) rin_=(3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
  pseudoB_ = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getPseudoB();
  if (pseudoB_ != 0.0)
  {
   if (aa_ == 0.0) 
   {
     rin_= 5.24392;
     x0_ = sqrt(rin_);
   }
   else if (aa_ == 0.1) 
   {
     rin_= 4.82365;
     x0_ = sqrt(rin_);
   }
   else if (aa_ == 0.2) 
   {
     rin_= 4.35976;
     x0_ = sqrt(rin_);
   }
   else if (aa_ == 0.3) 
   {
     rin_= 3.81529;
     x0_ = sqrt(rin_);
   }
   else if (aa_ == 0.4) 
   {
     rin_= 2.99911;
     x0_ = sqrt(rin_);
   }
   else
   {
    rin_= 1.334;
     x0_ = sqrt(rin_);
   }
//    else if (aa_ == 0.5) 
//    {
//      rin_= 1.43487;
//      x0_ = sqrt(rin_);
//    }
//    else if (aa_ == 0.6) 
//    {
//      rin_= 1.46850;
//      x0_ = sqrt(rin_);
//    }
//    else if (aa_ == 0.7) 
//    {
//      rin_= 1.502229;
//      x0_ = sqrt(rin_);
//    }
//    else if (aa_ == 0.8) 
//    {
//      rin_= 1.53476;
//      x0_ = sqrt(rin_);
//    }
//    else if (aa_ == 0.9) 
//    {
//      rin_= 1.56514;
//      x0_ = sqrt(rin_);
//    }
//    else if (aa_ == 1.0) 
//    {
//      rin_= 1.59268;
//      x0_ = sqrt(rin_);
//    }
  }
}

void PageThorneDisk::updateB() {
  pseudoB_ = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getPseudoB();
//  if (pseudoB_ != 0)  x0_ = 1.34;
}

void PageThorneDisk::setMetric(SmartPointer<Metric::Generic> gg) {
  if (gg_) gg_->unhook(this);
  string kind = gg->getKind();
  if (kind != "KerrBL" && kind != "KerrKS" && kind != "ChernSimons")
    throwError
      ("PageThorneDisk::setMetric(): metric must be KerrBL or KerrKS");
  ThinDisk::setMetric(gg);
  updateSpin();
  updateB();
  gg->hook(this);
}

double PageThorneDisk::emission(double nu_em, double dsem,
				    double *,
				    double coord_obj[8]) const{
  if (!blackbody_) {
    throwError("In PageThorneDisk::emission: "
	       "blackbody is necessary to compute emission, "
	       "else, use bolometricEmission");
  }

  double Ibolo=bolometricEmission(nu_em,dsem,coord_obj);
  /*
    From Ibolo, find T, then Bnu(T)
   */
  double mass=gg_->getMass()*1e3; // in cgs
  double c6=GYOTO_C_CGS*GYOTO_C_CGS*GYOTO_C_CGS
    *GYOTO_C_CGS*GYOTO_C_CGS*GYOTO_C_CGS;
  double g2m2=GYOTO_G_CGS*GYOTO_G_CGS*mass*mass;
  Ibolo*=mdot_*c6/g2m2; // Ibolo in cgs
  //F = sigma * T^4 (and F=pi*I)
  double TT=pow(Ibolo*M_PI/GYOTO_STEFANBOLTZMANN_CGS,0.25);
  spectrumBB_->setTemperature(TT);
  double Iem=(*spectrumBB_)(nu_em);
  //cout << "r T nu Iem = " << coord_obj[1] << " " << TT << " " << nu_em << " " << Iem << endl;
  if (Iem < 0.) throwError("In PageThorneDisk::emission"
			   " blackbody emission is negative!");
  return Iem;
}

Quantity_t PageThorneDisk::getDefaultQuantities() {
  return GYOTO_QUANTITY_USER4;
}

// Create the integrand for the Flux integral. This is not the complete function f of 
// Page and Thorne but only the integrand in the integral in f: (E-omega*L)d/dr(L).
double PageThorneDisk::integrand (double r, void * params) {
  // receive two parameters, spin a and pseudo-complex B
  double a = *((double *) params);
  double B = *((double *) params +1);
//   cout << "Spin a = "<< a << " --- B =" << B << endl;
  // Now the integrand  CForm of Mathematica output
  double f = (6*pow(a,3)*pow(r,1.5)*sqrt(1/(-3*B + 4*pow(r,2)))*(15*pow(B,2) - 32*B*pow(r,2) + 16*pow(r,4)) + 
     2*pow(r,3)*(15*pow(B,2) - 8*(-6 + r)*pow(r,4) - 2*B*pow(r,2)*(20 + 3*r)) - 
     4*pow(a,2)*(18*pow(B,2)*r + 4*(8 - 3*r)*pow(r,5) + 3*B*pow(r,3)*(-16 + 5*r)) + 
     a*(18*pow(B,2)*(10 - 7*r)*pow(r,2.5)*sqrt(1/(-3*B + 4*pow(r,2))) + 96*(2 - 3*r)*pow(r,6.5)*sqrt(1/(-3*B + 4*pow(r,2))) + 
        16*B*pow(r,4.5)*(-19 + 24*r)*sqrt(1/(-3*B + 4*pow(r,2))) - 45*pow(B,3)*sqrt(r/(-3*B + 4*pow(r,2)))))/
   (2.*(3*B - 4*pow(r,2))*(a + 2*sqrt(pow(r,5)/(-3*B + 4*pow(r,2))))*
     (4*pow(r,4)*(-3 + r + 4*a*sqrt(r/(-3*B + 4*pow(r,2)))) + B*(5*pow(r,2) - 12*a*sqrt(pow(r,5)/(-3*B + 4*pow(r,2))))));
  return f;
}



double PageThorneDisk::bolometricEmission(double nuem, double dsem,
				    double coord_obj[8]) const{
  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  // Important remark: this emision function gives I(r),
  // not I_nu(r). And I(r)/nu^4 is conserved.
  //  cout << "r hit= " << coord_obj[1] << endl;
  if (uniflux_) return 1;
  double xx;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    xx=sqrt(coord_obj[1]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xx=pow(coord_obj[1]*coord_obj[1]+coord_obj[2]*coord_obj[2]-aa2_, 0.25);
    break;
  default:
    throwError("Unknown coordinate system kind");
    xx=0;
  }
/*  
//*+++++++++++++++++++++++++++++++++++++++++++++++++++  
  // Ein pointer namens "w" des Typs gls_integration_workspace
// wird erzeugt und auf den Wert gsl_integration_workspace_alloc (1000) gesetzt
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, error;  
  double alpha = 1.0, expected = 4.0;

// Erzeugung eines Objekts(?) "F" des Types gsl_function
  gsl_function F;
// F.function wird auf die Adresse von f gesetzt.  
  F.function = &PageThorneDisk::f;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error); 
  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals =  %d\n", w->size);

  gsl_integration_workspace_free (w);
//*+++++++++++++++++++++++++++++++++++++++++++++++++++    
*/  
  
  // the above formula assume M=1 (x=sqrt(r/M)=sqrt(r))
  double x2=xx*xx;
  
/*  double ff=
    3./(2.)*1./(xx*xx*(xx*xx*xx-3.*xx+2.*aa_))
    *( 
      xx-x0_-3./2.*aa_*log(xx/x0_)
      -3.*(x1_-aa_)*(x1_-aa_)/(x1_*(x1_-x2_)*(x1_-x3_))*log((xx-x1_)
							    /(x0_-x1_)) 
      -3.*(x2_-aa_)*(x2_-aa_)/(x2_*(x2_-x1_)*(x2_-x3_))*log((xx-x2_)
							    /(x0_-x2_)) 
      -3.*(x3_-aa_)*(x3_-aa_)/(x3_*(x3_-x1_)*(x3_-x2_))*log((xx-x3_)
							    /(x0_-x3_))
       );
  // f of Page&Thorne, in units M=1
  
  double Iem=ff/(4.*M_PI*M_PI*x2);
  */
  
  /*
    with Mdot=1 (NB: Mdot is a constant)
    Assuming isotropic emission: 
    flux at r = (I at r)* \int cos\theta dOmega = pi*I
    thus intensity is only: 1/pi * flux
    NB: this is frequency integrated (bolometric) intensity, not I_nu
    NBB: the cgs value of I is c^6/G^2*Mdot/M^2 * Iem, it can be
    recovered from the dimensionless Iem a posteriori if needed
  */
  /*
   * TS: This stays the same. The determinant of the metric
   * is unaffected by the pc-parameter B and thus the connection
   * between intensity and flux stays unchanged.
   */

  // TS: Modification of ff. Need numerical integration here.  
  // Initializing a gsl integrator
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, result2, error;
  // Parameters to store a and B  
  double myPars[] = {aa_, pseudoB_};
  double *cur_par = myPars;
  
  // prefactor of the Integral (E-omega*L)^(-2)*d/dr(omega)
  double B= myPars[1];
  double r=x2;
  double prefactor =
  (-3*pow(r,1.5)*(-5*B + 4*pow(r,2))*sqrt(1/(-3*B + 4*pow(r,2))))/
   (4*pow(r,4)*(-3 + r + 4*aa_*sqrt(r/(-3*B + 4*pow(r,2)))) + B*(5*pow(r,2) - 12*aa_*sqrt(pow(r,5)/(-3*B + 4*pow(r,2)))));
  
//   if (prefactor > 0) cout << "Prefactor " << prefactor << " and the corresponding r: " << r << "\n";
  // create the gsl_function
  gsl_function F; 
  F.function = &PageThorneDisk::integrand;
  F.params = cur_par;  
  
  // Adapt the ISCO, if B!=0
  // This is ugly, as there is no closed form for the ISCO
  //Sqrt of ISCO = rms
  // Terms are rounded up to avoid singularities 
  // isco_classic = [6.0, 5.6693, 5.32944, 4.97862, 4.61434, 4.233, 3.82907, 3.39313, 
  //              2.90664, 2.32088, 1.0]
  double rms = x0_*x0_;
  /*   if (aa_ == 0.0) 
   {
     rms= 6.0;
   }
   else if (aa_ == 0.1) 
   {
     rms= 5.6693;
   }
   else if (aa_ == 0.2) 
   {
     rms= 5.32944;
   }
   else if (aa_ == 0.3) 
   {
     rms= 4.97862;
   }
   else if (aa_ == 0.4) 
   {
     rms= 4.61434;
   }
   else if (aa_ == 0.5) 
   {
     rms= 4.233;
   }
   else if (aa_ == 0.6) 
   {
     rms= 3.82907;
   }
   else if (aa_ == 0.7) 
   {
     rms= 3.39313;
   }
   else if (aa_ == 0.8) 
   {
     rms= 2.90664;
   }
   else if (aa_ == 0.9) 
   {
     rms= 2.32088;
   }
   else
   {
     rms= 1.0;     
   }
   */

   //TS: Adapt the inner radius of the integration
   // the radius is the radius, where the angular frequency has its maximum
   // really depends on the value of B, so be careful here. the value here
//    is for B= 64/27 m^3
   if (aa_ >= 0.5 && pseudoB_ != 0.0) 
   {
     rms= 1.72132593165;
   }

  // Do the integration
  // rms=2.05743  --- r= 32.8107
//   cout << "rms= " << rms << " --- r= "<< x2 << endl;
  gsl_integration_qags (&F, rms, x2, 0, 1e-7, 1000,
                        w, &result, &error); 
  
  gsl_integration_workspace_free (w);
  // TS 20. Januar. Fluss wird negativ....
//   result2 = - prefactor * result;
  result2 = -prefactor * result;
  double Iem= result2/(4.*M_PI*M_PI*x2);
// TS: end  
  
  if (flag_radtransf_) Iem *= dsem;
  GYOTO_DEBUG_EXPR(Iem);
//       cout << "Prefactor=" << prefactor << "\n";
//     cout << "Integrand=" << F.function << "\n";
  
//   if (Iem < 0){
//     cout << "Iem=" << Iem << "\n";
//   }
  return Iem;
  
}

void PageThorneDisk::processHitQuantities(Photon* ph, double* coord_ph_hit,
				     double* coord_obj_hit, double dt,
				     Properties* data) const {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  /*
      NB: freqObs is the observer's frequency chosen in
      Screen::getRayCoord for the actual computation of the geodesic ;
      the physical value of nuobs will be used in spectrum
      computations by resorting to the xml specifications of the user
      (see below) ; this freqObs is used to transform the null
      worldline parameter dlambda (see below)
  */
  double freqObs=ph->getFreqObs(); // this is a useless quantity, always 1
  SmartPointer<Spectrometer::Generic> spr = ph -> getSpectrometer();
  size_t nbnuobs = spr() ? spr -> getNSamples() : 0 ;
  double const * const nuobs = nbnuobs ? spr -> getMidpoints() : NULL;
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(coord_ph_hit,coord_obj_hit+4,
				    coord_ph_hit+4);// / 1.; 
                                       //this is nu_em/nu_obs
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  if (uniflux_) ggred=1.;
  double dsem = dlambda*ggredm1; // *1.
  double inc =0.;
//   cout << "data requested. " 
// 	      << ", ggredm1=" << ggredm1
// 	      << ", ggred=" << ggred
// 	      << endl;
// 	      
//   cout << "Koordinaten, 4er Geschwindigkeit, 4er Impuls der Photonen" << endl;
//   cout << "Koordinaten: " << "t: "<< coord_ph_hit[0] << " r: "<< coord_ph_hit[1] 
//        << " theta: "<< coord_ph_hit[2] << " phi: "  << coord_ph_hit[3] << endl;
//   cout << "4er Geschwindigkeit: " << "ut: "<< coord_obj_hit[4] << " ur: "<< coord_obj_hit[5] 
//        << " utheta: "<< coord_obj_hit[6] << " uphi: "  << coord_obj_hit[7] << endl;
//   cout << "4er Impuls: " << "pt: "<< coord_ph_hit[4] << " pr: "<< coord_ph_hit[5] 
//        << " ptheta: "<< coord_ph_hit[6] << " pphi: "  << coord_ph_hit[7] << endl;
       

  if (data) {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "data requested. " 
	      << ", ggredm1=" << ggredm1
	      << ", ggred=" << ggred
	      << endl;
#endif

    if (data->redshift) {
      *data->redshift=ggred;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->redshift);
#endif
    }
    if (data->time) {
      *data->time=coord_ph_hit[0];
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->time);
#endif
    }
    if (data->impactcoords) {
      memcpy(data->impactcoords, coord_obj_hit, 8 * sizeof(double));
      memcpy(data->impactcoords+8, coord_ph_hit, 8 * sizeof(double));
    }
#if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "dlambda = (dt="<< dt << ")/(tdot="<< coord_ph_hit[4]
		<< ") = " << dlambda << ", dsem=" << dsem << endl;
#endif
    if (data->intensity) throwError("unimplemented");
    if (data->user4) {
      inc = (bolometricEmission(freqObs*ggredm1, dsem, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* ggred*ggred*ggred*ggred; // I/nu^4 invariant
      *data->user4 += inc;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->user4);
#endif

    }
    if (data->binspectrum) throwError("unimplemented");
    if (data->spectrum)  {
      if (!blackbody_) {
	throwError("In PageThorneDisk::process: "
		   "blackbody is necessary to compute spectrum");
      }
      for (size_t ii=0; ii<nbnuobs; ++ii) {
	double nuem=nuobs[ii]*ggredm1;
	double * tmp;
	inc = (emission(nuem, dsem, tmp, coord_obj_hit))
	  * (ph -> getTransmission(size_t(-1)))
	  * ggred*ggred*ggred; // Inu/nu^3 invariant
	data->spectrum[ii*data->offset] += inc;
	//cout << "in spec stored= " << ggred << " " << inc << endl;
      }
    }
    /* update photon's transmission */
    ph -> transmit(size_t(-1),
		   transmission(freqObs*ggredm1, dsem,coord_ph_hit));
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "NO data requested!" << endl;
#   endif
  }
}

void PageThorneDisk::tell(Hook::Teller* msg) {
  updateSpin();
  updateB();
}

int PageThorneDisk::setParameter(std::string name,
				 std::string content,
				 std::string unit) {
  char* tc = const_cast<char*>(content.c_str());
  if (name=="BlackbodyMdot") {
    blackbody_=1;
    mdot_=atof(tc);
  }
  if (name=="UniFlux") uniflux_=1;
  else return ThinDisk::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void PageThorneDisk::fillElement(FactoryMessenger *fmp) const {
  fmp->setMetric(gg_);
  ThinDisk::fillElement(fmp);
}
#endif
