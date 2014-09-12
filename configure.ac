#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_COPYRIGHT((c) the Regents of the University of California)
AC_REVISION($Revision$)
AC_PREREQ(2.56)

AC_INIT(FastQuery, 0.8.2.5, fastbit-users@hpcrdm.lbl.gov)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/fq.cpp])
AC_CONFIG_HEADER([src/fastquery-config.h])
AC_DEFINE(PACKAGE_URL, "http://sdm.lbl.gov/fastbit/", [FastQuery web site])
AM_INIT_AUTOMAKE([no-define])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
#AC_PROG_RANLIB # `AC_PROG_RANLIB' is rendered obsolete by `LT_INIT'
#AC_PROG_INSTALL
AC_PROG_LIBTOOL
#LT_INIT([all-static])
#LT_INIT

# Checks for libraries.
AC_CHECK_LIB([m], [log])
#AC_CHECK_LIB([rt], [clock_gettime]) // needed by FastBit
#AC_CHECK_LIB([z], [main]) // needed by HDF5, NetCDF

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h stdint.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
#AC_C_CONST
#AC_C_INLINE
#AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
#AC_CHECK_FUNCS([floor memset])

FQ_PARA_MODE="yes"; # default to include MPI library
FQ_TRY_HDF5="yes"; # default to enable HDF5
FQ_TRY_NETCDF="no";
FQ_TRY_PNETCDF="no";
FQ_TRY_BP="no";
DEFAULT_FB_VERSION="fastbit-ibis1.3.3"; # default fastbit version

AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE(debug,
              [AC_HELP_STRING([--enable-debug], [turn on debugging (default: set by autoconf)])],
              [if test "$enableval" = "yes" ; then debug="yes" ; else debug="no" ; fi],[debug=""])
AC_MSG_RESULT($debug)
AC_MSG_CHECKING(whether to enable extra optimization)
AC_ARG_ENABLE(xopt,
              [AC_HELP_STRING([--enable-xopt], [turn on extra optimization (default: no, leave the default CXXFLAGS alone)])],
              [if test "$enableval" = "yes" ; then xopt="yes" ; else xopt="no" ; fi],[xopt="no"])
AC_MSG_RESULT($xopt)
if test "x$debug" = "xyes" ; then
  if test `echo $CFLAGS | fgrep -- '-g' | wc -l` -eq 0 ; then
    CFLAGS="$CFLAGS -g"
  fi
  if test `echo $CXXFLAGS | fgrep -- '-g' | wc -l` -eq 0 ; then
    CXXFLAGS="$CXXFLAGS -g"
  fi
elif test "x$debug" = "xno" ; then
  if test `echo $CFLAGS | fgrep -- '-g' | wc -l` -gt 0 ; then
    OIFS="$IFS"
    IFS=' '
    newflags=""
    for g in $CFLAGS ; do if test "x$g" != "x-g[0123]" -a "x$g" != "x-g" ; then newflags="$newflags $g"; fi; done
    CFLAGS="$newflags"
    IFS=$OIFS
  fi
  if test `echo $CXXFLAGS | fgrep -- '-g' | wc -l` -gt 0 ; then
    OIFS="$IFS"
    IFS=' '
    newflags=""
    for g in $CXXFLAGS ; do if test "x$g" != "x-g[0123]" -a "x$g" != "x-g" ; then newflags="$newflags $g"; fi; done
    CXXFLAGS="$newflags"
    IFS=$OIFS
  fi
fi
if test "x$xopt" = "xyes" ; then
   # add -O5 if the compiler accepts it, otherwise -O3
   OLDCXXFLAGS="$CXXFLAGS";
   CXXFLAGS="$CXXFLAGS -O5";
   AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([[
        #include <cmath>
        ]], [[
        return (1 == static_cast<int>(std::sqrt(2.0)));
        ]])], [], [CXXFLAGS="$OLDCXXFLAGS -O3"])
fi
if test "x$xopt" = "xno" ; then
  if test `echo $CFLAGS | fgrep -- '-O' | wc -l` -gt 0 ; then
    OIFS=$IFS
    IFS=' '
    newflags=""
    for g in $CFLAGS ; do if test "x$g" != "x-O" -a "x$g" != "x-O2" -a "x$g" != "x-O3" ; then newflags="$newflags $g"; fi; done
    CFLAGS=$newflags
    IFS=$OIFS
  fi
  if test `echo $CXXFLAGS | fgrep -- '-O' | wc -l` -gt 0 ; then
    OIFS=$IFS
    IFS=' '
    newflags=""
    for g in $CXXFLAGS ; do if test "x$g" != "x-O" -a "x$g" != "x-O2" -a "x$g" != "x-O3" ; then newflags="$newflags $g"; fi; done
    CXXFLAGS=$newflags
    IFS=$OIFS
  fi
fi
  
# check for static or dynamic library built
#AC_ARG_ENABLE(static,[AC_HELP_STRING([--enable-static])[build with only static library links]],[CXXFLAGS="-static -fPIC"])

# check for FastBit [required]
AC_ARG_WITH(fastbit,
	AC_HELP_STRING([--with-fastbit=dir], [use dir as the root directory for FastBit include and lib.  If the directory is empty, the configure script will attempt to download a copy with wget or curl.  If a directory is not specified, FastBit will be installed in the directory for this package.]),
	[if test -d "$with_fastbit/include" && test -d "$with_fastbit/lib" && test -x $with_fastbit/bin/fastbit-config; then
	  FASTBIT_HOME="$with_fastbit"
	fi])

# check for default path for FastBit library if no valid path is given by the with-fastbit option
if test "x$FASTBIT_HOME" = "x"; then
   fcpath=`which fastbit-config`
   fcpathlength=`echo ${fcpath} | wc -c`
   if test ${fcpathlength} -gt 20; then
	default_path=`${fcpath} --prefix`
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
    fi
fi
if test "x$FASTBIT_HOME" = "x"; then
	default_path="$FASTBIT_ROOT"
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
fi
if test "x$FASTBIT_HOME" = "x"; then
	default_path="${HOME}/`uname -s`/fastbit"
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
fi
if test "x$FASTBIT_HOME" = "x"; then
	default_path="${prefix}/fastbit"
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
fi
if test "x$FASTBIT_HOME" = "x"; then
	default_path="DEFAULT_FB_VERSION"
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
fi
if test "x$FASTBIT_HOME" = "x"; then
	default_path="${HOME}/DEFAULT_FB_VERSION"
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
fi
if test "x$FASTBIT_HOME" = "x"; then
	default_path="/project/projectdirs/m1248/${NERSC_HOST}/pgi/fastbit" # FastBit location on NERSC machines
	if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib && test -x $default_path/bin/fastbit-config; then
	  FASTBIT_HOME="$default_path";
	fi
fi


AC_MSG_CHECKING([for FastBit])
if test "x$FASTBIT_HOME" != "x"; then
   fbflags="`${FASTBIT_HOME}/bin/fastbit-config --cflags`"
   CFLAGS="$CFLAGS $fbflags"
   CPPFLAGS="$CPPFLAGS $fbflags"
   CXXFLAGS="$CXXFLAGS $fbflags"
   CXXCPPFLAGS="$CXXCPPFLAGS $fbflags"
   LIBS="$LIBS `${FASTBIT_HOME}/bin/fastbit-config --libs`";
   AC_TRY_LINK([#include <ibis.h>],
   	       [ibis::init();],
	       [AC_MSG_RESULT([yes])],
	       [download_fastbit=1;
	        AC_MSG_WARN([FastBit is required, try to download a copy])])
fi
if test "x$download_fastbit" != "x" ; then
   if test -x "/usr/bin/curl"; then
	/usr/bin/curl -O https://codeforge.lbl.gov/frs/download.php/399/fastbit-ibis1.3.3.tar.gz
   elif test -x "/usr/bin/wget"; then
	/usr/bin/wget https://codeforge.lbl.gov/frs/download.php/399/fastbit-ibis1.3.3.tar.gz
   elif test -x "/usr/local/bin/curl"; then
	/usr/local/bin/curl -O https://codeforge.lbl.gov/frs/download.php/399/fastbit-ibis1.3.3.tar.gz
   elif test -x "/usr/local/bin/wget"; then
	/usr/local/bin/wget https://codeforge.lbl.gov/frs/download.php/399/fastbit-ibis1.3.3.tar.gz
    fi
    if test -f fastbit-ibis1.3.3.tar.gz; then
	echo "Configuring and building FastBit, will take about half a hour";
	tar xf fastbit-ibis1.3.3.tar.gz
	if test "x$FASTBIT_HOME" == "x" ; then
	    FASTBIT_HOME="${prefix}"/fastbit
        fi
	cd fastbit-ibis1.3.3
	./configure --prefix=${FASTBIT_HOME} --disable-debug
	make -j 3
	make install
	cd ..
    fi
    if test "x$FASTBIT_HOME" != "x"; then
       fbflags="`${FASTBIT_HOME}/bin/fastbit-config --cflags`"
       CFLAGS="$CFLAGS $fbflags"
       CPPFLAGS="$CPPFLAGS $fbflags"
       CXXFLAGS="$CXXFLAGS $fbflags"
       CXXCPPFLAGS="$CXXCPPFLAGS $fbflags"
       LIBS="$LIBS `${FASTBIT_HOME}/bin/fastbit-config --libs`";
       AC_TRY_LINK([#include <ibis.h>],
     	       [ibis::init();],
	       [AC_MSG_RESULT([yes])],
	       [AC_MSG_ERROR([FastBit is required, but failed to link the library])])
    else
	AC_MSG_ERROR([
FastBit is required, please download the source code from http://sdm.lbl.gov/fastbit/])
    fi
fi

# check for MPI 
AC_ARG_ENABLE(parallel, [AC_HELP_STRING([--enable-parallel], [enable MPI parallel version (default: yes)])],
	[if test "$enableval" != "no"; then FQ_PARA_MODE="yes";
	 else
	    FQ_PARA_MODE="no";
	    CPPFLAGS="${CPPFLAGS} -D_NOMPI -DFQ_NOMPI";
	    CXXFLAGS="${CXXFLAGS} -D_NOMPI -DFQ_NOMPI";
	    CFLAGS="${CFLAGS} -D_NOMPI -DFQ_NOMPI";
	 fi])

AC_ARG_WITH(mpi,
	AC_HELP_STRING([--with-mpi=dir], [use dir as the root directory for MPI include and lib]),
	[if test -d "$with_mpi" && test -d $with_mpi/include && test -d $with_mpi/lib; then
	    	FQ_PARA_MODE="yes";
	  	MPI_HOME="$with_mpi";
	else
		MPI_HOME="";
	fi]
)

if test "x$FQ_PARA_MODE" != "xno"; then
# check for default path for MPI library if no valid path is given from the with-mpi option
if test "x$MPI_HOME" = "x"; then
	default_path="$MPICH_DIR";
	[if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib; then
	  MPI_HOME=$default_path;
	fi]
fi
if test "x$MPI_HOME" = "x"; then
	default_path="$MPI_ROOT";
	[if test -d "$default_path" && test -d $default_path/include && test -d $default_path/lib; then
	  MPI_HOME=$default_path;
	fi]
fi
if test "x$MPI_HOME" = "x"; then
	default_path="mpich2-gnu"
	[if test -d $default_path && test -d $default_path/include && test -d $default_path/lib; then
	  MPI_HOME=$default_path;
	fi]
fi
if test "x$MPI_HOME" = "x"; then
	default_path="/usr/common/mpich2-gnu"
	[if test -d $default_path && test -d $default_path/include && test -d $default_path/lib; then
	  MPI_HOME=$default_path;
	fi]
fi
if test "x$MPI_HOME" = "x"; then
	default_path="/opt/cray/mpt/5.1.4/xt/gemini/mpich2-gnu" # MPI location on NERSC Hopper
	[if test -d $default_path && test -d $default_path/include && test -d $default_path/lib; then
	  MPI_HOME=$default_path;
	fi]
fi
if test "x$MPI_HOME" = "x"; then
	default_path="/opt/cray/mpt/default/xt/seastar/mpich2-gnu" # MPI location on NERSC Franklin
	[if test -d $default_path && test -d $default_path/include && test -d $default_path/lib; then
	  MPI_HOME=$default_path;
	fi]
fi

AC_MSG_CHECKING([for MPI])
if test -d "$MPI_HOME" -a -d "$MPI_HOME/include" -a -d "$MPI_HOME/lib"; then 
	OLDLIBS="$LIBS";
   	OLDCFLAGS="$CFLAGS";
	OLDCPPFLAGS="$CPPFLAGS";
	OLDCXXFLAGS="$CXXFLAGS";
	OLDCXXCPPFLAGS="$CXXCPPFLAGS";
	LIBS="$LIBS -L$MPI_HOME/lib -lmpich -lm"; # try mpich
	CFLAGS="$CFLAGS -I$MPI_HOME/include";
	CPPFLAGS="$CPPFLAGS -I$MPI_HOME/include";
	CXXFLAGS="$CXXFLAGS -I$MPI_HOME/include";
	CXXCPPFLAGS="$CXXCPPFLAGS -I$MPI_HOME/include";
  	AC_TRY_LINK([#include <mpi.h>], [MPI_Finalize();],
		    [AC_MSG_RESULT([yes])],
		    [LIBS="$OLDLIBS -L$MPI_HOME/lib -lmpi -lm"; # try openmpi
		     CFLAGS="$OLDCFLAGS -I$MPI_HOME/include/openmpi";
		     CPPFLAGS="$OLDCPPFLAGS -I$MPI_HOME/include/openmpi";
		     CXXFLAGS="$OLDCXXFLAGS -I$MPI_HOME/include/openmpi";
		     CXXCPPFLAGS="$OLDCXXCPPFLAGS -I$MPI_HOME/include/openmpi";
		     AC_TRY_LINK([#include <mpi.h>], [MPI_Finalize();],
		     	         [AC_MSG_RESULT([yes])],
		     		 [LIBS="$OLDLIBS"; # try no additional flags
		     		  CFLAGS="$OLDCFLAGS";
		     		  CPPFLAGS="$OLDCPPFLAGS";
		     		  CXXFLAGS="$OLDCXXFLAGS";
		     		  CXXCPPFLAGS="$OLDCXXCPPFLAGS";
		     		  AC_TRY_LINK([#include <mpi.h>], [MPI_Finalize();],
		     	           	      [AC_MSG_RESULT([yes])],
		     		 	      [FQ_PARA_MODE="no";
						AC_MSG_WARN([no: failed to link with the MPI library! 
Can only produce NOMPI version of the library, make sure MPI is not needed anywhere.])
					      ])])])
else
	# try if it works without any additional flags
	AC_TRY_LINK([#include <mpi.h>], [MPI_Finalize();],
		    [AC_MSG_RESULT([yes])],
		    [OLDLIBS="$LIBS";
		     LIBS="$LIBS -lmpi";
		     AC_TRY_LINK([#include <mpi.h>], [MPI_Finalize();],
		     		 [AC_MSG_RESULT([yes])],
				 [LIBS="$OLDLIBS -lmpich";
		     		  AC_TRY_LINK([#include <mpi.h>], [MPI_Finalize();],
				  [AC_MSG_RESULT([yes])],
				  [FQ_PARA_MODE="no";
                                   AC_MSG_WARN([no: failed to link with the MPI library!
Can only produce NOMPI version of the library, make sure MPI is not needed anywhere.])
				  ])])])
fi
fi

# check for HDF5
AC_ARG_ENABLE(hdf5, [AC_HELP_STRING([--enable-hdf5], [enable HDF5 file format interface (default: enable)])],
	[if test "$enableval" = "yes"; then FQ_TRY_HDF5="yes"; else FQ_TRY_HDF5="no"; fi])

AC_ARG_WITH(hdf5,
	AC_HELP_STRING([--with-hdf5=dir], [use dir as the root directory for HDF5 include and lib]),
	[if test -d "$with_hdf5" -a -d "$with_hdf5/include" -a -d "$with_hdf5/lib"; then
		HDF5_HOME="$with_hdf5";
		FQ_TRY_HDF5="yes";
	elif test "$with_hdf5" = "no"; then
	     	FQ_TRY_HDF5="no";
		HDF5_HOME=""
	else
		HDF5_HOME=""
	fi]
)
AC_ARG_WITH(szlib,
	[AC_HELP_STRING([--with-szlib], [the SZip library used by HDF5])],
	[if test -f "${with_szlib}/include/szlib.h" -a -d "${with_szlib}/lib"; then SZLIB_HOME="$with_szlib";
	 else AC_MSG_WARN([Option "--with-szlib=${with_szlib}" does not specify a valid directory]); fi])

if test "x$FQ_TRY_HDF5" = "xyes"; then
# check for default path for HDF5 library if no valid path is given from the with-hdf5 option
if test "x$HDF5_HOME" = "x"; then
   	if test -d "$HDF5_ROOT"; then
	   default_path="$HDF5_ROOT"
	elif test -d "$HDF5_DIR"; then
	   default_path="$HDF5_DIR"
	else
	   default_path="$CRAY_HDF5_DIR"
	fi
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  HDF5_HOME="$default_path";
	fi]
fi
if test "x$HDF5_HOME" = "x"; then
	default_path="hdf5-1.8.5.0" 
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  HDF5_HOME="$default_path";
	fi]
fi
if test "x$HDF5_HOME" = "x"; then
	default_path="/usr/common/hdf5-1.8.5.0" 
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  HDF5_HOME=$default_path;
	fi]
fi
if test "x$HDF5_HOME" = "x"; then
	default_path="/opt/cray/hdf5-parallel/1.8.5.0/hdf5-parallel-gnu" # HDF5 location on NERSC Hopper
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  HDF5_HOME=$default_path;
	fi]
fi
if test "x$HDF5_HOME" = "x"; then
	default_path="/opt/cray/hdf5-parallel/default" # HDF5 location on NERSC Franklin
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  HDF5_HOME=$default_path;
	fi]
fi
if test "x$HDF5_HOME" = "x"; then
	default_path="/opt/cray/hdf5/default" # HDF5 location on NERSC Franklin
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  HDF5_HOME=$default_path;
	fi]
fi

AC_MSG_CHECKING([for HDF5])
if test -d "$HDF5_HOME" && test -d "$HDF5_HOME/include" && test -d "$HDF5_HOME/lib"; then
	OLDLIBS="$LIBS"
   	OLDCFLAGS="$CFLAGS";
	OLDCPPFLAGS="$CPPFLAGS";
	OLDCXXFLAGS="$CXXFLAGS";
	OLDCXXCPPFLAGS="$CXXCPPFLAGS";
    	CFLAGS="$CFLAGS -I$HDF5_HOME/include"
    	CPPFLAGS="$CPPFLAGS -I$HDF5_HOME/include"
    	CXXFLAGS="$CXXFLAGS -I$HDF5_HOME/include"
    	CXXCPPFLAGS="$CXXCPPFLAGS -I$HDF5_HOME/include"
    	LIBS="$LIBS -L$HDF5_HOME/lib -lhdf5 -lhdf5_hl -lz"
	if test "x$SZLIB_HOME" != "x"; then
	    CFLAGS="$CFLAGS -I$SZLIB_HOME/include"
	    CPPFLAGS="$CPPFLAGS -I$SZLIB_HOME/include"
	    CXXFLAGS="$CXXFLAGS -I$SZLIB_HOME/include"
	    CXXCPPFLAGS="$CXXCPPFLAGS -I$SZLIB_HOME/include"
	    LIBS="$LIBS -L$HDF5_HOME/lib -lsz"
	fi
    	AC_TRY_LINK([#include <hdf5.h>], [H5Fclose(0);],
		    [AC_MSG_RESULT([yes]); FQ_VERIFIED_HDF5="yes"],
		    [AC_MSG_RESULT([no: failed to link the HDF5 library!]);
			LIBS="$OLDLIBS";
			CFLAGS="$OLDCFLAGS";
			CPPFLAGS="$OLDCPPFLAGS";
			CXXFLAGS="$OLDCXXFLAGS";
			CXXCPPFLAGS="$OLDCXXCPPFLAGS";])
else
	AC_MSG_RESULT([no: HDF5 directory not found!])
fi
fi

# check for NETCDF
AC_ARG_ENABLE(netcdf,
	[AC_HELP_STRING([--enable-netcdf], [enable NETCDF file format interface (default: disabled)])],
	[if test "$enableval" = "yes"; then FQ_TRY_NETCDF="yes"; else FQ_TRY_NETCDF="no"; fi])

AC_ARG_WITH(netcdf,
	AC_HELP_STRING([--with-netcdf=dir], [use dir as the root directory for NETCDF include and lib]),
	[if test -d "$with_netcdf" && test -d "$with_netcdf/include" && test -d "$with_netcdf/lib"; then
		NETCDF_HOME="$with_netcdf"
		FQ_TRY_NETCDF="yes";
	elif test "$with_netcdf" = "no"; then
	     	FQ_TRY_NETCDF="no";
		NETCDF_HOME=""
	else
		NETCDF_HOME=""
	fi]
)
if test "x$FQ_TRY_NETCDF" = "xyes"; then
# check for default path for NETCDF library if no valid path is given from the with-netcdf option
if test "x$NETCDF_HOME" = "x"; then
   	if test -d "$NETCDF_ROOT"; then
	    default_path="$NETCDF_ROOT"
	elif test -d "$NETCDF_DIR"; then
	    default_path="$NETCDF_DIR"
	else
	    default_path="$CRAY_NETCDF_ROOT"
	fi
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  NETCDF_HOME="$default_path";
	fi]
fi
if test "x$NETCDF_HOME" = "x"; then
	default_path="netcdf-4.1.1.0" 
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  NETCDF_HOME=$default_path;
	fi]
fi
if test "x$NETCDF_HOME" = "x"; then
	default_path="/usr/common/netcdf-4.1.1.0" 
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  NETCDF_HOME=$default_path;
	fi]
fi
if test "x$NETCDF_HOME" = "x"; then
	default_path="/opt/cray/netcdf/4.1.1.0/netcdf-gnu" # NETCDF location on NERSC Hopper
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  NETCDF_HOME=$default_path;
	fi]
fi
if test "x$NETCDF_HOME" = "x"; then
	default_path="/opt/cray/netcdf-hdf5parallel/default" # netCDF location on NERSC Franklin
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  NETCDF_HOME=$default_path;
	fi]
fi
if test "x$NETCDF_HOME" = "x"; then
	default_path="/opt/cray/netcdf/default" # netCDF location on NERSC Franklin
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  NETCDF_HOME=$default_path;
	fi]
fi

AC_MSG_CHECKING([for NETCDF])
if test -d "$NETCDF_HOME" && test -d "$NETCDF_HOME/include" && test -d "$NETCDF_HOME/lib"; then
	OLDLIBS="$LIBS"
   	OLDCFLAGS="$CFLAGS";
	OLDCPPFLAGS="$CPPFLAGS";
	OLDCXXFLAGS="$CXXFLAGS";
	OLDCXXCPPFLAGS="$CXXCPPFLAGS";
    	CFLAGS="$CFLAGS -I$NETCDF_HOME/include"
    	CPPFLAGS="$CPPFLAGS -I$NETCDF_HOME/include"
    	CXXFLAGS="$CXXFLAGS -I$NETCDF_HOME/include"
    	CXXCPPFLAGS="$CXXCPPFLAGS -I$NETCDF_HOME/include"
    	LIBS="-L$NETCDF_HOME/lib -lnetcdf $LIBS"
    	AC_TRY_LINK([#include "netcdf.h"], [nc_close(0);],
		    [AC_MSG_RESULT([yes]); FQ_VERIFIED_NETCDF="yes"; check="yes"],
		    [AC_MSG_RESULT([no: failed to link the NETCDF library!])
		    LIBS="$OLDLIBS"
   		    CFLAGS="$OLDCFLAGS";
		    CPPFLAGS="$OLDCPPFLAGS";
		    CXXFLAGS="$OLDCXXFLAGS";
		    CXXCPPFLAGS="$OLDCXXCPPFLAGS";])
else
	AC_MSG_RESULT([no: NETCDF directory not found!])
fi
fi

# check for PNETCDF
AC_ARG_ENABLE(pnetcdf,
	[AC_HELP_STRING([--enable-pnetcdf], [enable PNETCDF file format interface (default: disabled)])],
	[if test "$enableval" = "yes"; then FQ_TRY_PNETCDF="yes"; else FQ_TRY_PNETCDF="no"; fi])

AC_ARG_WITH(pnetcdf,
	AC_HELP_STRING([--with-pnetcdf=dir], [use dir as the root directory for PNETCDF include and lib]),
	[if test -d "$with_pnetcdf" && test -d "$with_pnetcdf/include" && test -d "$with_pnetcdf/lib"; then
		PNETCDF_HOME="$with_pnetcdf"
		FQ_TRY_PNETCDF="yes";
	elif test "$with_pnetcdf" = "no"; then
	     	FQ_TRY_PNETCDF="no";
		PNETCDF_HOME=""
	else
		PNETCDF_HOME=""
	fi]
)

if test "x$FQ_TRY_PNETCDF" = "xyes"; then
# check for default path for PNETCDF library if no valid path is given from the with-netcdf option
if test "x$PNETCDF_HOME" = "x"; then
   	if test -d "$PNETCDF_ROOT"; then
	    default_path="$PNETCDF_ROOT"
	elif test -d "$PNETCDF_DIR"; then
	    default_path="$PNETCDF_DIR"
	else
	    default_path="$CRAY_PNETCDF_ROOT"
	fi
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  PNETCDF_HOME="$default_path";
	fi]
fi
if test "x$PNETCDF_HOME" = "x"; then
	default_path="pnetcdf-1.2.0" 
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  PNETCDF_HOME=$default_path;
	fi]
fi

AC_MSG_CHECKING([for PNETCDF])
if test -d "$PNETCDF_HOME" && test -d "$PNETCDF_HOME/include" && test -d "$PNETCDF_HOME/lib"; then
	OLDLIBS="$LIBS"
   	OLDCFLAGS="$CFLAGS";
	OLDCPPFLAGS="$CPPFLAGS";
	OLDCXXFLAGS="$CXXFLAGS";
	OLDCXXCPPFLAGS="$CXXCPPFLAGS";
    	CFLAGS="$CFLAGS -I$PNETCDF_HOME/include"
    	CPPFLAGS="$CPPFLAGS -I$PNETCDF_HOME/include"
    	CXXFLAGS="$CXXFLAGS -I$PNETCDF_HOME/include"
    	CXXCPPFLAGS="$CXXCPPFLAGS -I$PNETCDF_HOME/include"
    	LIBS="-L$PNETCDF_HOME/lib -lpnetcdf $LIBS"
    	AC_TRY_LINK([#include "pnetcdf.h"], [ncmpi_close(0);],
		    [AC_MSG_RESULT([yes]); FQ_VERIFIED_PNETCDF="yes"; check="yes"],
		    [AC_MSG_RESULT([no: failed to link the PNETCDF library!])
		    LIBS="$OLDLIBS"
   		    CFLAGS="$OLDCFLAGS";
		    CPPFLAGS="$OLDCPPFLAGS";
		    CXXFLAGS="$OLDCXXFLAGS";
		    CXXCPPFLAGS="$OLDCXXCPPFLAGS";])
else
	AC_MSG_RESULT([no: PNETCDF directory not found!])
fi
fi

# check for BP
AC_ARG_ENABLE(adios,
	[AC_HELP_STRING([--enable-bp], [enable ADIOS BP file format interface (default: disable)])],
	[if test "$enableval" = "yes"; then FQ_TRY_BP="yes"; else FQ_TRY_BP="no"; fi])

AC_ARG_ENABLE(bp,
	[AC_HELP_STRING([--enable-bp], [enable ADIOS BP file format interface (default: disable)])],
	[if test "$enableval" = "yes"; then FQ_TRY_BP="yes"; else FQ_TRY_BP="no"; fi])

AC_ARG_WITH(adios,
	AC_HELP_STRING([--with-adios=dir], [use dir as the root directory for ADIOS BP include and lib]),
	[if test -d "$with_adios" && test -d "$with_adios/include" && test -d "$with_adios/lib"; then
		BP_HOME="$with_adios"
		FQ_TRY_BP="yes"
	elif test "$with_adios" = "no"; then
	     	FQ_TRY_BP="no";
		BP_HOME=""
	else
		BP_HOME=""
	fi]
)

AC_ARG_WITH(bp,
	AC_HELP_STRING([--with-bp=dir], [use dir as the root directory for ADIOS BP include and lib]),
	[if test -d "$with_bp" && test -d "$with_bp/include" && test -d "$with_bp/lib"; then
		BP_HOME="$with_bp"
		FQ_TRY_BP="yes"
	else
		BP_HOME=""
	fi]
)

AC_ARG_WITH(mxml,
	AC_HELP_STRING([--with-mxml=dir], [use dir as the root directory for MXML include and lib, assumed to be the same as ADIOS directory if no specified]),
	[if test -d "$with_mxml" && test -d "$with_mxml/include" && test -d "$with_mxml/lib"; then
		MXML_HOME="$with_mxml"
		#echo "MXML_HOME=${MXML_HOME}"
	fi]
)

if test "x$FQ_TRY_BP" = "xyes"; then		
# check for default path for BP library if no valid path is given from the with-bp option
if test "x$BP_HOME" = "x" && test "x$ADIOS_DIR" != "x"; then
	default_path="$ADIOS_DIR"
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  BP_HOME="$default_path";
	fi]
fi
if test "x$BP_HOME" = "x" && test "xBP_ROOT" != "x"; then
	default_path="$BP_ROOT"
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  BP_HOME="$default_path";
	fi]
fi
if test "x$BP_HOME" = "x"; then
	default_path="/global/homes/j/jinohkim/sys/local/adios" # BP location on NERSC Hopper
	[if test -d "$default_path" && test -d "$default_path/include" && test -d "$default_path/lib"; then
	  BP_HOME=$default_path;
	fi]
fi

AC_MSG_CHECKING([for ADIOS BP])
if test -x "$BP_HOME/bin/adios_config" && test "x$FQ_PARA_MODE" != "xno"; then
    BP_INC=`$BP_HOME/bin/adios_config -c`
    BP_LIB=`$BP_HOME/bin/adios_config -l`
elif test -d "$BP_HOME/include" && test -d "$BP_HOME/lib"; then
    BP_INC="-I$BP_HOME/include"
    if test "x$FQ_PARA_MODE" = "xno" ; then
    	BP_LIB="-L$BP_HOME/lib -ladios_nompi -ladiosread_nompi";
    else
        BP_LIB="-L$BP_HOME/lib -ladios -ladiosread"
    fi
    if test -d "$MXML_HOME/include" && test -d "$MXML_HOME/lib"; then
	BP_INC="$BP_INC -I$MXML_HOME/include"
	BP_LIB="$BP_LIB -L$MXML_HOME/lib -lmxml";
    else
	BP_LIB="$BP_LIB -lmxml";
    fi
else
    if test "x$FQ_PARA_MODE" = "xno" ; then
    	BP_LIB="-ladios_nompi -ladiosread_nompi -lmxml";
    else
        BP_LIB="-ladios -ladiosread -mxml"
    fi
fi

OLDLIBS="$LIBS"
OLDCFLAGS="$CFLAGS";
OLDCPPFLAGS="$CPPFLAGS";
OLDCXXFLAGS="$CXXFLAGS";
OLDCXXCPPFLAGS="$CXXCPPFLAGS";
CFLAGS="$CFLAGS $BP_INC"
CPPFLAGS="$CPPFLAGS $BP_INC"
CXXFLAGS="$CXXFLAGS $BP_INC"
CXXCPPFLAGS="$CXXCPPFLAGS $BP_INC"
LIBS="$LIBS $BP_LIB";
AC_TRY_LINK([#include "adios.h"], [adios_close(0);],
	    [AC_MSG_RESULT([yes]); FQ_VERIFIED_BP="yes"],
	    [AC_MSG_RESULT([no: failed to link the BP library!])
	    LIBS="$OLDLIBS";
	    CFLAGS="$OLDCFLAGS";
	    CPPFLAGS="$OLDCPPFLAGS";
	    CXXFLAGS="$OLDCXXFLAGS";
	    CXXCPPFLAGS="$OLDCXXCPPFLAGS";])
fi

if test "x$FQ_VERIFIED_NETCDF" = "xyes"; then
    AC_DEFINE(FQ_HAVE_NETCDF, 1, [Define to 1 if NETCDF format is supported.])
fi
if test "x$FQ_VERIFIED_PNETCDF" = "xyes"; then
    AC_DEFINE(FQ_HAVE_PNETCDF, 1, [Define to 1 if PNETCDF format is supported.])
fi
if test "x$FQ_VERIFIED_HDF5" = "xyes"; then
    AC_DEFINE(FQ_HAVE_HDF5, 1, [Define to 1 if HDF5 format is supported.])
fi
if test "x$FQ_VERIFIED_BP" = "xyes"; then
    AC_DEFINE(FQ_HAVE_BP, 1, [Define to 1 if ADIOS BP format is supported.])
fi
if test "x$FQ_VERIFIED_HDF5" != "xyes" -a "x$FQ_VERIFIED_NETCDF" != "xyes" -a "x$FQ_VERIFIED_BP" != "xyes" -a "x$FQ_VERIFIED_PNETCDF" != "xyes"; then
    AC_MSG_ERROR([FastQuery needs to have at least one of HDF5, netCDF, PnetCDF or ADIOS BP, none found!])
fi

AM_CONDITIONAL(_NOMPI_ONLY, test "x$FQ_PARA_MODE" = "xno")
AM_CONDITIONAL(BUILD_NETCDF, test "x$FQ_VERIFIED_NETCDF" = "xyes")
AM_CONDITIONAL(BUILD_PNETCDF, test "x$FQ_VERIFIED_PNETCDF" = "xyes")
AM_CONDITIONAL(BUILD_HDF5, test "x$FQ_VERIFIED_HDF5" = "xyes")
AM_CONDITIONAL(BUILD_BP, test "x$FQ_VERIFIED_BP" = "xyes")

OBJDIR="$objdir"
AC_SUBST(_NOMPI_ONLY)
AC_SUBST(BUILD_HDF5)
AC_SUBST(BUILD_NETCDF)
AC_SUBST(BUILD_PNETCDF)
AC_SUBST(BUILD_BP)
AC_SUBST(LIB_VERSION_INFO)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(CXXCPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(DEFS)
AC_SUBST(OBJDIR)
AC_SUBST([LIBTOOL_DEPS])
#end stuff not generated by autoscan

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])

AC_OUTPUT
