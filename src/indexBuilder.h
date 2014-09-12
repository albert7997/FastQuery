#ifndef _FQ_INDEXBUILDER_H
#define _FQ_INDEXBUILDER_H

#include <vector>
#include <string>
#include "ibis.h"
#include "fq.h"
#include "fqPart.h"

/*!
  \brief Open a file and provide indexing and data modification
  functions to it.
  The data file must be already exist, and it will be opend with write
  permission.

  \param dataFileName (IN): Name of a file containing data.

  \param dataModel (IN): Enumerate data format model for the file.

  \param indexFileName (IN): Name of a file containing indexes.

  \param v (IN): verbosenes for print out message.

  \param rcfile (IN): IN Name of the runtime control file for FastQuery and
  FastBit.

  \param logfile (IN): Name of the log file for FastQuery and FastBit.

  \param strm (IN): whether to use the streaming mode in ADIOS.  This
  parameter only affects operations of ADIOS files.

  \param comm (IN): MPI communication channel for the file.  It is only
  used for MPI parallel mode.

*/

class IndexBuilder : public FastQuery {
public:
    /*******************************************
     * Constructor and Destructor
     ********************************************/
#ifdef FQ_NOMPI
    IndexBuilder(const std::string &dataFileName,
                 const FQ::FileFormat dataModel,
                 const std::string &indexFileName,
                 const int v=ibis::gVerbose,
                 const char *rcfile=0,
                 const char *logfile=0,
                 void *extra = 0);
#else
    IndexBuilder(const std::string &dataFileName,
                 const FQ::FileFormat dataModel,
                 const std::string &indexFileName,
                 const int v=ibis::gVerbose,
                 const char *rcfile=0,
                 const char *logfile=0,
                 const MPI_Comm comm=MPI_COMM_WORLD,
                 void *extra = 0);
#endif

    ~IndexBuilder(){};
    /*******************************************
     * Building Index APIs
     ********************************************/
/*!
  \brief Build bitmap indices for all variables who have a particular name
  or path string.
        This must be a collective call from all MPI tasks.

  \param binning IN: Binning option (optional)
  The format for the binning option is '\<binning nbins=xxxx /\>'
  where xxxx is the number of bins.
  If no binning is specified, a bitmap index is built for
  each distinct data value.
  Details on the binning option can be found at:
  http://sdm.lbl.gov/fastbit/doc/indexSpec.html

  \param varPathStr IN: A variable path string in the file structure. (optional)
  If varPathStr is not given, it means any variable path.
  \param varNameStr IN: A variable name string in the file structure. (optional)
  If varNameStr is not given, it means any variable name.
  \param mpi_dim IN: The dimension to split the data into subarray in the MPI mode.
  \param mpi_len IN: The subarray size for splitting the data in the MPI mode.
  \return The number of variables whose indexes are available in the file.

  Given a file structure:
\pre
  "/"
  /   |    \
  "time0"   "c"    "time1"
  /  \      |     /  |  \
  "a"  "b"   "d"  "a" "b" "c"
\endpre

  buildIndexes(): builds indexes for all variables for default binning option.
  buildIndexes(binning, "", "/time0/a"): builds indexes for variables "/time0/a".
  buildIndexes(binning, "/time0/a"): builds indexes for variables "/time0/a".
  buildIndexes(binning, "time0"): builds indexes for variables "/time0/a" and "/time0/b".
  buildIndexes(binning, "", "a"): builds indexes for variables "/time0/a" and "/time1/a".
  buildIndexes(binning, "c"): builds indexes for variables "/c/d" and "/time1/c".
*/

    int buildIndexes(const char *binning=0,
                     const std::string &varPathStr="",
                     const std::string &varNameStr="",
                     uint64_t mpi_dim=FQ_DEFAULT_MPI_DIM,
                     uint64_t mpi_len=FQ_DEFAULT_MPI_LEN);

    /*******************************************
     * Extend File I/O APIs handled by ArrayIODriver
     ********************************************/
/*!
  \brief Create a new variable
        This must be a collective call from all MPI tasks

  \param variable IN: The full-path name of the dataset for the variable.

  \param dims IN: Vector containing the dimensions of the variable.

  \param type IN: Enumerated type of the variable.

  \return True if success, otherwise return False.
*/

    bool createNewVariable(const std::string &variable,
                           std::vector<uint64_t> dims,
                           FQ::DataType type);

/*!
  \brief Insert data values for a prticular variable

  \param varNameStr IN: A variable name string in the file structure.

  \param data IN: Data values.

  \param varPathStr IN: A variable path string in the file structure. (optional)
  If varPathStr is not given, it means any variable path.

  \param collective IN: Enable collective call from all MPI tasks (default: true).

  \return True if success, otherwise return False.
*/

    bool setData(const std::string &varNameStr,
                 const void* data,
                 const std::string &varPathStr="",
                 const bool collective=true);

/*!
  \brief Create a vraiable attribute and set its values.
        This must be a collective call from all MPI tasks.

  \param varNameStr IN: A variable name string in the file structure.

  \param attrNameStr IN: A attribute name of the given variable.

  \param values IN: Attribute data values. The memory size of the data must
  be reserved.

  \param length IN: The number of values for the attribute.

  \param fqType IN:  Enumerated type of variable.

  \param varPathStr IN: A variable path string in the file
  structure. (optional) If varPathStr is not given, it means any variable
  path.

  \return True if success, otherwise return False.
*/

    bool setAttribute(const std::string &varNameStr,
                      const std::string &attrName,
                      const void *values,
                      const uint64_t length,
                      const FQ::DataType fqType,
                      const std::string &varPathStr="");
};
#endif
