AUTOMAKE_OPTIONS=gnu

EXTRA_LTLIBRARIES = libfastquery_nompi.la
libfastquery_nompi_la_SOURCES = fqColumn.cpp fqIndex.cpp fqPart.cpp indexBuilder.cpp fq.cpp fqParser.cpp metadataMgr.cpp queryProcessor.cpp queryParser.cc queryLexer.cc

lib_LTLIBRARIES = libfastquery.la
#libfastquery_la_LDFLAGS = -version-info $(LIB_VERSION_INFO)
libfastquery_la_SOURCES = fqColumn.cpp fqIndex.cpp fqPart.cpp indexBuilder.cpp fq.cpp fqParser.cpp metadataMgr.cpp queryProcessor.cpp queryParser.cc queryLexer.cc

# if _NOMPI_ONLY
# all: nompi
# endif

if BUILD_HDF5
libfastquery_la_SOURCES += hdf5file.cpp
libfastquery_nompi_la_SOURCES += hdf5file.cpp
endif

if BUILD_NETCDF
libfastquery_la_SOURCES += netCDFfile.cpp
libfastquery_nompi_la_SOURCES += netCDFfile.cpp
endif

if BUILD_PNETCDF
libfastquery_la_SOURCES += pnetCDFfile.cpp
endif

if BUILD_BP
libfastquery_la_SOURCES += ADIOS_Wrapper.cpp BPArrayIODriver.cpp
libfastquery_nompi_la_SOURCES += ADIOS_Wrapper.cpp BPArrayIODriver.cpp
endif

include_HEADERS = fastquery-config.h const.h FlexLexer.h fq.h fqParser.h fqVar.h indexBuilder.h queryLexer.h arrayIODriver.h fqColumn.h fqIndex.h fqPart.h metadataMgr.h qToken.h queryProcessor.h queryParser.hh stack.hh position.hh location.hh

if BUILD_HDF5
include_HEADERS += hdf5file.h
endif

if BUILD_NETCDF
include_HEADERS += netCDFfile.h
endif

if BUILD_PNETCDF
libfastquery_la_SOURCES += pnetCDFfile.h
endif

if BUILD_BP
include_HEADERS += ADIOS_Wrapper.h BPArrayIODriver.h
endif

EXTRA_DIST=queryLexer.ll queryParser.yy

# if _NOMPI_ONLY
# libfastquery.la:
# 	rm -f libfastquery.la
# 	echo can NOT create parallel version of FastQuery library libfastquery.la
# else
libfastquery.la: $(libfastquery_la_SOURCES) $(libfastquery_la_DEPENDENCIES) $(EXTRA_libfastquery_la_DEPENDENCIES)
	if test -f libfastquery_nompi.la; then rm -f *.o *.lo; fi
	$(MAKE)	$(libfastquery_la_OBJECTS)
	$(CXXLINK) -rpath $(libdir) $(libfastquery_la_OBJECTS) $(libfastquery_la_LIBADD) $(LIBS)
# endif

libfastquery_nompi.la: $(libfastquery_nompi_la_SOURCES) $(libfastquery_nompi_la_DEPENDENCIES) $(EXTRA_libfastquery_nompi_la_DEPENDENCIES) 
	if test -f libfastquery.la; then rm -f *.o *.lo; fi
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DFQ_NOMPI" CFLAGS="$(CFLAGS) -DFQ_NOMPI" $(libfastquery_nompi_la_OBJECTS)
	$(CXXLINK) -rpath $(libdir) $(libfastquery_nompi_la_OBJECTS) $(libfastquery_nompi_la_LIBADD) $(LIBS)

install-exec-local:
	$(NORMAL_INSTALL)
	test -z "$(libdir)" || $(MKDIR_P) "$(DESTDIR)$(libdir)"
	if test -f libfastquery_nompi.la; then $(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) libfastquery_nompi.la "$(DESTDIR)$(libdir)"; fi

clean-am: clean-generic clean-libLTLIBRARIES clean-libtool  mostlyclean-am
	rm -f $(EXTRA_LTLIBRARIES)

docs: Doxyfile
	perl -pi -e 's/(PROJECT_NUMBER\s*=\s*).+/\1$(PACKAGE_VERSION)/' Doxyfile
	doxygen

nompi: libfastquery_nompi.la

.PHONY: nompi
