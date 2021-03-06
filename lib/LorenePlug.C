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

// include Metric headers
#include "GyotoRotStar3_1.h"
/*
namespace Gyoto {
  void __StdLibInit();
}
*/

extern "C" void __GyotoloreneInit() {
  Gyoto::Metric::Register("RotStar3_1",
			  &(Gyoto::Metric::Subcontractor
			    <Gyoto::Metric::RotStar3_1>));
}
