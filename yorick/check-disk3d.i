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

#include "gyoto.i"
gyoto_debug, 0;
#include "gyoto_std.i"

if (get_env("GYOTO_CHECK_NODISPLAY")) nodisplay = 1;

emissquant=array(double, 1, 2, 10 ,10);
emissquant(,1,,1:3)=100.;
emissquant(,2,,4:10)=100.;

velocity=array(1., 3, 2, 10, 10);

metric = gyoto_KerrBL(mass=4e6*GYOTO_SUN_MASS);

write, format="%s", "Creating Disk3D...";
pd = gyoto_Disk3D(copyemissquant=emissquant,
                  copyvelocity=velocity,
                  rin=3, rout=28,
                  zmin=1., zmax=10.,
                  phimin=0., phimax=2.*pi,
                  repeatphi=8,
                  metric=metric);
write, format="%s\n", " done.";

write, format="%s\n", "Printing Disk3D:";
pd;
write, format="%s\n", " done.";

screen = gyoto_Screen(metric=metric, resolution=64,
                      time=1000.*metric.unitlength/GYOTO_C,
                      distance=100.*metric.unitlength, fov=30./100.,
                      inclination=110./180.*pi, paln=pi);

write, format="%s", "Attaching Disk3D to scenery...";
sc = gyoto_Scenery(metric=metric, screen=screen, astrobj=pd);
write, format="%s\n", " done.";

write, format="%s", "Saving data to fits file...";
pd, fitswrite="!check-disk3d.fits.gz";
write, format="%s\n", " done.";

write, format="%s", "Saving scenery to XML file...";
sc, xmlwrite="check-disk3d.xml";
write, format="%s\n", " done.";

write, format="%s", "Reading back scenery...";
sc2 = gyoto_Scenery("check-disk3d.xml");
write, format="%s\n", " done.";

write, format="%s", "Removing temporary files...";
remove, "check-disk3d.xml";
remove, "check-disk3d.fits.gz";
write, format="%s\n", " done.";

write, format="%s", "Getting Disk3D...";
pd2 = sc2.astrobj;
write, format="%s\n", " done.";

write, format="%s", "Comparing emissquant array...";
if (anyof(emissquant != pd2.copyemissquant)) error, "CHECK FAILED";
write, format="%s\n", " done.";

write, format="%s", "Comparing velocity array...";
if (anyof(velocity != pd2.copyvelocity)) error, "CHECK FAILED";
write, format="%s\n", " done.";

/*write, format="%s", "Performing raytracing...\n";
im = sc();
write, format="%s\n", "done.";

if (!nodisplay) {
  write, format="%s", "Displaying image...";
  window, style="nobox.gs";
  pli, im;
  write, format="%s\n", " done.";
  pause, 1000;
  if (batch()) winkill;
 }
*/
if (batch()) quit;
