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

#include "ygyoto.h"
#include "yapi.h"
#include <iostream>
#include "GyotoFactory.h"
#include "GyotoBlackBodySpectrum.h"

namespace YGyoto {
  namespace Spectrum {
    ygyoto_Spectrum_eval_worker_t BlackBodyYEval;
  }
}

using namespace Gyoto;
using namespace YGyoto;
using namespace Gyoto::Spectrum;
using namespace YGyoto::Spectrum;
using namespace std;

void YGyoto::Spectrum::BlackBodyYEval(SmartPointer<Generic> * OBJ_, int argc) {

  static char const * knames[]={
    "unit",
    "temperature", "scaling",
    YGYOTO_SPECTRUM_GENERIC_KW,
    0
  };

  YGYOTO_WORKER_INIT(Spectrum, BlackBody,
		     knames, YGYOTO_SPECTRUM_GENERIC_KW_N+3);

  YGYOTO_WORKER_SET_UNIT;

  YGYOTO_WORKER_GETSET_DOUBLE(Temperature);
  YGYOTO_WORKER_GETSET_DOUBLE(Scaling);
  YGYOTO_WORKER_CALL_GENERIC(Spectrum);

}

extern "C" {
  void
  Y_gyoto_BlackBodySpectrum(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT(Spectrum, BlackBody);
    if ((*OBJ)->getKind().compare("BlackBody"))
      y_error("Expecting Spectrum of kind BlackBody");
    BlackBodyYEval(OBJ, argc);
  }

  void Y__gyoto_BlackBodySpectrum_register_as_Metric(){
    ygyoto_Spectrum_register("BlackBody",&BlackBodyYEval);
  }
}
