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

#include "GyotoUtils.h"
#include "GyotoPhoton.h"
#include "GyotoFixedStar.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstring>
#include <float.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

FixedStar::FixedStar() : UniformSphere("FixedStar")
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  for (int i=0;i<3;++i) pos_[i]=0.;
}

FixedStar::FixedStar(SmartPointer<Gyoto::Metric::Generic> gg, double StPsn[3],
		     double rad) :
  UniformSphere("FixedStar", gg, rad)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "(metric, pos, rad)" << endl;
# endif
  for (int i=0;i<3;++i) pos_[i] = StPsn[i]; 
  setRadius(rad);
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done" << endl;
# endif
}

FixedStar::FixedStar(const FixedStar& orig) :
  UniformSphere(orig)
{
  for (int i=0; i<3; ++i) pos_[i] = orig.pos_[i];
}
FixedStar* FixedStar::clone() const { return new FixedStar(*this); }

FixedStar::~FixedStar() {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
}

void FixedStar::getCartesian(double const * const , size_t const n_dates,
			     double * const x, double * const y,
			     double * const z, double * const xprime,
			     double * const yprime, 
			     double * const zprime) {
  double xs, ys, zs;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    xs= pos_[0];
    ys= pos_[1];
    zs= pos_[2];
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    {
      double rs=pos_[0];
      double ths=pos_[1];
      double phs=pos_[2];
      double st, ct, sp, cp;
      sincos(ths, &st, &ct);
      sincos(phs, &sp, &cp);
      xs= rs*st*cp;
      ys= rs*st*sp;
      zs= rs*ct;
    }
    break;
  default:
    throwError("unsupported coordkind");
    xs=ys=zs=0.;
  }
  for (size_t i=0; i<n_dates; ++i) {
    if (x) x[i] = xs;
    if (y) y[i] = ys;
    if (z) z[i] = zs;
    if (xprime) xprime[i] = 0.;
    if (yprime) yprime[i] = 0.;
    if (zprime) zprime[i] = 0.;
  }
}

void FixedStar::getVelocity(double const pos[4], double vel[4]) {
  for (size_t i=0; i<4; ++i) vel[i]=0.;
  vel[0]=gg_->SysPrimeToTdot(pos, vel+1);
}

double const * FixedStar::getPos() const { return pos_; }

void FixedStar::getPos(double dst[3]) const
{ for (int i=0; i<3;++i) dst[i]=pos_[i]; }

void FixedStar::setMetric(SmartPointer<Metric::Generic> gg) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
 Generic::setMetric(gg);
 setRadius(radius_);
}

void FixedStar::setRadius(double r) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(r) ;
# endif
  radius_ = r;
  critical_value_=r*r;
  safety_value_=1.1*critical_value_;
  if (!gg_()) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "metric is not set yet" << endl;
#   endif
    return;
  }
  switch (gg_ -> getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rmax_=3.*(pos_[0]+radius_);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    rmax_=3.*(sqrt(pos_[0]*pos_[0]+pos_[1]*pos_[1]+pos_[2]*pos_[2])+radius_);
    break;
  default:
    throwError("unimplemented coordinate system in FixedStar");
  } 
}

void FixedStar::setPos(const double p[3])
{ for (int i=0; i<3; ++i) pos_[i]=p[i]; setRadius(radius_);}

int FixedStar::setParameter(string name, string content, string unit) {
  if (name=="Position") {
    char* tc = const_cast<char*>(content.c_str());
    double pos[3];
    for (int i=0;i<3;++i) pos[i] = strtod(tc, &tc);
    setPos(pos);
  } else return UniformSphere::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void FixedStar::fillElement(FactoryMessenger *fmp) const {
  fmp -> setParameter ("Position", const_cast<double*>(pos_), 3);
  UniformSphere::fillElement(fmp);
}
#endif
