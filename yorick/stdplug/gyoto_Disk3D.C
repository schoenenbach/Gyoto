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

#include <cstring>

#include <Gyoto.h>
#include <GyotoFactory.h>
#include "../ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

#define OBJ ao

// on_eval worker
void ygyoto_Disk3D_eval(SmartPointer<Astrobj::Generic> *ao_, int argc) {

  static char const * knames[]={
    "unit",
    "fitsread", "repeatphi", "nu0", "dnu",
    "rin", "rout", "zmin", "zmax",
    "phimin", "phimax",
    "copyemissquant", "copyvelocity",
    "fitswrite",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Astrobj, Disk3D, knames, YGYOTO_ASTROBJ_GENERIC_KW_N+14);

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_RUN( fitsRead(ygets_q(iarg)) );
  YGYOTO_WORKER_GETSET_LONG2(repeatPhi);
  YGYOTO_WORKER_GETSET_DOUBLE2(nu0);
  YGYOTO_WORKER_GETSET_DOUBLE2(dnu);
  YGYOTO_WORKER_GETSET_DOUBLE2(rin);
  YGYOTO_WORKER_GETSET_DOUBLE2(rout);
  YGYOTO_WORKER_GETSET_DOUBLE2(zmin);
  YGYOTO_WORKER_GETSET_DOUBLE2(zmax);
  YGYOTO_WORKER_GETSET_DOUBLE2(phimin);
  YGYOTO_WORKER_GETSET_DOUBLE2(phimax);

  /* EMISSQUANT */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyemissquant=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[4];
      (*ao) -> getEmissquantNaxes(ddims);
      long dims[] = {4, ddims[0], ddims[1], ddims[2], ddims[3]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getEmissquant(),
	     dims[1]*dims[2]*dims[3]*dims[4]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyEmissquant(NULL, 0);
      else if (dims[0]==4) {
	size_t ddims[] = {dims[1], dims[2], dims[3], dims[4]};
	(*ao)->copyEmissquant(in, ddims);
      } else
	y_error("COPYEMISSQUANT must be nil, 0, or array(double, nnu, nphi, nz, nr");
    }
  }

  /* VELOCITY */

  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyvelocity=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[4];
      (*ao) -> getEmissquantNaxes(ddims);
      long dims[] = {4, 3, ddims[1], ddims[2], ddims[3]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getVelocity(),
	     3*dims[2]*dims[3]*dims[4]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyVelocity(NULL, 0);
      else if (dims[0]==4 && dims[1]==3) {
	size_t ddims[] = {dims[2], dims[3], dims[4]};
	(*ao)->copyVelocity(in, ddims);
      } else
	y_error("COPYVELOCITY must be nil, 0, or array(double, 3, nphi, nz, nr");
    }
  }

  YGYOTO_WORKER_RUN( fitsWrite(ygets_q(iarg)) );

  YGYOTO_WORKER_CALL_GENERIC(Astrobj);
}

extern "C" {
  void Y__gyoto_Disk3D_register_as_Astrobj(){
    ygyoto_Astrobj_register("Disk3D",&ygyoto_Disk3D_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_Disk3D(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT(Astrobj, Disk3D);
    if ((*ao)->getKind().compare("Disk3D"))
	y_error("Expecting Astrobj of kind Disk3D");
    ygyoto_Disk3D_eval(ao, argc);
  }

}
