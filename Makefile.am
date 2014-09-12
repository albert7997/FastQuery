AUTOMAKE_OPTION=gnu
ACLOCAL_AMFLAGS=-I m4
SUBDIRS=src examples

EXTRA_DIST=runconfig README test-hdf5.sh test-netcdf.sh test-pnetcdf.sh

docs:
	cd src && $(MAKE) docs

examples:
	cd examples && $(MAKE) examples

clean: clean-recursive
	-rm -fr autom4te.cache bin/* lib/* 

nompi:
	cd src && $(MAKE) nompi

.PHONY: docs examples nompi

