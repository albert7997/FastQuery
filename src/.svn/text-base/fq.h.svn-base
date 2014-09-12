/*! \mainpage FastQuery

  \author Jerry Chou, Kesheng Wu, Surendra Byna

  \copyright FastQuery Copyright (c) 2011 - 2012, The Regents of the
             University of California, through Lawrence Berkeley National
             Laboratory (subject to receipt of any required approvals from
             the U.S. Dept. of Energy).  This software was developed under
             funding from Dept. of Energy, Office of Science, Office of
             Advanced Scientific Computing Research (ASCR).  All rights
             reserved.

  FastQuery provides a simple API for indexing and querying datasets using
  the FastBit bitmap index technology. Datasets can be retrieved using
  complex compound range queries such as "(energy > 100) && (70 < pressure
  < 90)". The bitmap index technology only retrieves the data elements that
  satisfy the query condition and shows significant speedup compared with
  reading the entire datasets.  Dataset can be managed by a list of file
  system or data model including HDF5, NetCDF, etc.

  Given a file structure:
  @verbatim
                  "/"
             /     |      \
        "time0"   "c"    "time1"
         /  |  \  /     /   |  \
       "a" "b" "d"     "a" "b" "c"
  @endverbatim

  \section definition Variable definition

  Every unique dataset location in a file is treated as an individual
  varaiable in FastQuery.  Since a dataset in a file is allowed to be
  connected by multiple links, a dataset could match to multiple variables,
  but a variable only matches to a single dataset.  E.g: variables in the
  file are {"/time0/a", "/time0/b", "/time0/d", "/c/d", "/time1/a",
  "/time1/b", "/time1/c"}

  \section naming Variable naming schema

  Each variable is identified by its file location which can be described
  by a tuple(varPathStr, varNameStr).  varNameStr is a file location string
  that ends with the name of the dataset of a variable.  If varNameStr is
  not specified, it means any variable name.  varPathStr is a file location
  string that contains the sub-path to the dataset of a variable.  If
  varNameStr is not specified, it means any variable path.  If multiple
  variables match to the tuple specification, the best fit variable is the
  one with the shortest path length.  E.g:

@verbatim
        ("","")           => all variables.
        ("", "a")         => "/time0/a", "/time1/a".
        ("time0", "a")    => "/time0/a", "/time0/b", "/time0/d"
        ("", "/time0/a")  => "/time0/a".
        ("c", "")         => "/c/d".
        ("", "c")         => "/time1/c".
@endverbatim
*/

#ifndef _FASTQUERY_H
#define _FASTQUERY_H

#include "const.h"
#include "metadataMgr.h"
#include "arrayIODriver.h"
#include <ibis.h>
#include <vector>
#include <string>

#ifndef FQ_NOMPI
#include <mpi.h>
#endif

/*******************************************
* Basic FastQuery APIs
********************************************/
class FastQuery {
public:
    /*!
      Open a file and provide basic variable information retriving functions.
       The file must be already exist.

      \param dataFileName IN: Name of a file containing data.
      \param ffmt IN: Enumerate data format model for the file.
      \param indexFileName IN: Name of a file containing indexes.
      \param v: verbosenes for print out message.
      \param rcfile: Name of the runtime control file for FastQuery and
             FastBit.
      \param logfile: Name of the log file for FastQuery and FastBit
      \param readOnly IN: Whether to open the index file in read only
             mode.  The data file is always opened in read-only mode.
      \param extra IN: extra options to be passed to the concrete array IO
             drivers.  Currently, only ADIOS driver makes use of this
             argument.
      \param comm IN: MPI communication channel for the file.
             It is only used for MPI parallel mode.
    */
#ifdef FQ_NOMPI
    FastQuery(const std::string &dataFileName,
              const FQ::FileFormat ffmt,
              const std::string &indexFileName="",
              const int v=ibis::gVerbose,
              const char *rcfile=0,
              const char *logfile=0,
              bool readOnly=true,
              void *extra=0);
#else
    FastQuery(const std::string &dataFileName,
              const FQ::FileFormat ffmt,
              const std::string &indexFileName="",
              const int v=ibis::gVerbose,
              const char *rcfile=0,
              const char *logfile=0,
              bool readOnly=true,
              MPI_Comm comm=MPI_COMM_WORLD,
              void *extra=0);
#endif
    /*!
     \brief Close files and clean out memory usage.
    */
    ~FastQuery();

    /*!
      Get the full-path of all variables which contains a specific path
      string or name string.

      \param variables OUT: A vector containing the full-path of the
      variables found.
      \param varPathStr IN: A variable path string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.
      \param varNameStr IN: A variable name string in the file
      structure. (optional)
      If varNameStr is not given, it means any variable name.
      \return the number of variables found.

      Given following file structure,
      @verbatim
                                "/"
                             /   |    \
                      "time0"   "c"    "time1"
                       /  \      |     /  |  \
                     "a"  "b"   "d"  "a" "b" "c"
        @endverbatim
       - getAllVariables(variables): find all variables.
       - getAllVariables(variables,"time0"): find "/time0/a", "/time0/b".
       - getAllVariables(variables,"","a"): find "/time0/a", "/time1/a".
       - getAllVariables(variables,"","time/a"): find "/time0/a".
       - getAllVariables(variables,"","c"): find "/time1/c".
       - getAllVariables(variables,"c"): find "/c/d".
       - getAllVariables(variables,"/time0","a"): find "/time0/a".

    */
    unsigned int getAllVariables(std::vector<std::string> &variables,
                                 const std::string &varPathStr="",
                                 const std::string &varNameStr="");

    /*!
      \brief Check if the FastQuery object is initiated successfully.
    */
    bool isValid();

    /*!
      \brief Retrieve information of a particular variable.
        If multiple variables match to the given variable specifications,
        pick the best fit variable with the shortest path-name.

      \param varNameStr IN: A variable name string in the file structure.
      \param dims OUT: Vector containing the dimensions for a given variable.
      \param type OUT: Enumerated type of variable
      \param varPathStr IN: A variable name string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.

      \return false if no variables match to the given varNameStr and
        varPathStr, or failed to retrieve the variable information.
        Otherwise, return true.
    */
    bool getVariableInfo(const std::string &varNameStr,
                         std::string &variable,
                         std::vector <uint64_t> &dims,
                         FQ::DataType *type,
                         const std::string &varPathStr="");

    /*!
      \brief Query the existence of a particular variable

      \param varNameStr IN: A variable name string in the file structure.
      \param varPathStr IN: A variable path string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.

      \return false if no variables match to the given varNameStr and
        varPathStr.  Otherwise return false.
    */
    bool checkForVariable(const std::string &varNameStr,
                          const std::string &varPathStr="");

    /*!
      \brief Get all the data values of a variable.
        If multiple variables match to the given variable specifications,
        pick the best fit variable with the shortest path-name.

      \param varNameStr IN: A variable name string in the file structure.
      \param data OUT: Data values. The memory size of the data must be
      reserved.
      \param varPathStr IN: A variable path string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.
      \param collective IN: Enable collective call from all MPI tasks (default: true).
      \param bcast IN: Broadcast data to all MPI tasks (default: true).

      \return True if success, otherwise return False.
    */
    bool getData(const std::string &varNameStr,
                 void *data,
                 const std::string &varPathStr="",
                 const bool collective=true,
                 const bool bcast=true);

    /*!
      \brief Get all the data values of a variable attribute.

      \param varNameStr IN: A variable name string in the file structure.
      \param attrNameStr IN: A attribute name of the given variable.
      \param values OUT: Attribute data values. The memory size of the data
      must be reserved.
      \param varPathStr IN: A variable path string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.

      \return True if success, otherwise return False.
    */

    bool getAttribute(const std::string &varNameStr,
                      const std::string &attrName,
                      void *values,
                      const std::string &varPathStr="");

    /*!
      \brief Retrieve information of a variable attribute.
        If multiple variables match to the given variable specifications,
        pick the best fit variable with the shortest path-name.

      \param varNameStr IN: A variable name string in the file structure.
      \param attrNameStr IN: A attribute name of the given variable.
      \param length OUT: The number of values for the variable attribute.
      \param type OUT: Enumerated type of variable
      \param varPathStr IN: A variable path string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.

      \return True if success, otherwise return False.
    */
    bool getAttributeInfo(const std::string &varNameStr,
                          const std::string &attrName,
                          uint64_t *length,
                          FQ::DataType *type,
                          const std::string &varPathStr="");

    static bool reportTiming() {return report_timing;}

protected:
/*******************************************
* Members objects
********************************************/
    // ArrayIODriver implements the BaseArrayIODriver to enable file I/O
    ArrayIODriver* dataFile;
    ArrayIODriver* indexFile;
    MetadataMgr* metadataMgr;
    bool valid;
    bool isValid(const std::string &func) const;

    /// Shall we turn on the internal timers?  This is set to true if one of
    /// the following condition is satisfied.  The default value is false.
    /// - the macro DEBUG is specified and it greater than 0,
    /// - the macro FQ_REPORT_STATISTIC is defined and the run-time
    /// configuration parameter with the name is true,
    /// - when a constructor is called, the verboseness level is greater
    /// than 3,
    /// - the parameter named "FastQuery.reportTiming" is true.
    static bool report_timing;

#ifndef FQ_NOMPI
    MPI_Comm mpi_comm;
    int mpi_rank;
    int mpi_size;
#endif
};

namespace util {
    double compute_median(double* arr, int num, bool allow_modify=false);
};
#endif
