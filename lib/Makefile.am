ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCESCPPFLAGS)
AM_LDFLAGS  = $(XERCESLDFLAGS) -export-dynamic
AM_CXXFLAGS = @use_xerces@ @use_cfitsio@ \
	-DGYOTO_PLUGIN_SFX=\"@plugin_sfx@\" -DGYOTO_PREFIX=\"@prefix@\" \
	-DGYOTO_SOVERS=\"@sovers@\"

# HEADERS: where they are, where to install them
library_includedir=$(includedir)/Gyoto
library_include_HEADERS=$(top_srcdir)/include/*.h

# MAIN LIBRARY
lib_LTLIBRARIES = libgyoto.la
libgyoto_la_SOURCES =  Astrobj.C Factory.C Register.C SmartPointer.C Utils.C Metric.C	\
	Worldline.C FocalPlane.C Photon.C Scenery.C			\
	WorldlineIntegState.C Error.C Screen.C Spectrum.C		\
	Spectrometer.C ThinDisk.C
libgyoto_la_LIBS = $(XERCESLDFLAGS)
libgyoto_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(VERSINFO)

# STANDARD PLUGIN
soverdir=$(pkglibdir)/@sovers@
sover_LTLIBRARIES = libgyoto-stdplug.la
libgyoto_stdplug_la_SOURCES =  KerrBL.C KerrKS.C \
	Star.C FixedStar.C Torus.C \
	PowerLawSpectrum.C BlackBodySpectrum.C \
	ComplexAstrobj.C UniformSphere.C \
	StandardAstrobj.C PageThorneDisk.C \
	ThinDiskPL.C \
	StdPlug.C
# those files are used only if cfitsio is available
cfitsio_stdplug_sources = PatternDisk.C PatternDiskBB.C DynamicalDisk.C \
	 Disk3D.C Disk3D_BB.C
EXTRA_libgyoto_stdplug_la_SOURCES = $(cfitsio_stdplug_sources)
if HAVE_CFITSIO
libgyoto_stdplug_la_SOURCES  +=	$(cfitsio_stdplug_sources)
endif
libgyoto_stdplug_la_LDFLAGS = -module -export-dynamic $(CFITSIOLDFLAGS)

# LORENE PLUGIN
if HAVE_LORENE
sover_LTLIBRARIES += libgyoto-lorene.la
endif
libgyoto_lorene_la_SOURCES = RotStar3_1.C LorenePlug.C
libgyoto_lorene_la_LDFLAGS = -module -export-dynamic $(LORENELDFLAGS)
libgyoto_lorene_la_CPPFLAGS = $(AM_CPPFLAGS) $(LORENECPPFLAGS)

# pkg-config file
pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA=gyoto.pc