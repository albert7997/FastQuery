#ifndef _FQ_QUERYPROCESSOR_H
#define _FQ_QUERYPROCESSOR_H

#include <vector>
#include <string>
#include "ibis.h"
#include "fq.h"
#include "fqPart.h"

/// Main FQ class with functionality for users to perform queries on files
/// and get varaible data.
class QueryProcessor : public FastQuery {
public:
    /*******************************************
     * Constructor and Destructor
     ********************************************/
#ifdef FQ_NOMPI
    QueryProcessor(const std::string &dataFileName,
                   const FQ::FileFormat dataModel,
                   const std::string &indexFileName="",
                   const int v=ibis::gVerbose,
                   const char *rcfile=0,
                   const char *logfile=0,
                   void *extra = 0);
#else
    QueryProcessor(const std::string &dataFileName,
                   const FQ::FileFormat dataModel,
                   const std::string &indexFileName="",
                   const int v=ibis::gVerbose,
                   const char *rcfile=0,
                   const char *logfile=0,
                   const MPI_Comm comm=MPI_COMM_WORLD,
                   void *extra = 0);
#endif

    ~QueryProcessor(){};
    /*******************************************
     * Querying Index APIs
     ********************************************/
    /*!
      \brief Get data values from a variable by a list of selected coordinates.
        If multiple variables match to the given variable specifications,
        pick the best fit variable with the shortest path-name.
             This is only supported as an independent call for all MPI tasks.

      \param varNameStr IN: A variable name string in the file structure.
      \param coords IN: A list of selected coordinates.
      \param data OUT: Data values. The memory size of the data must be
      reserved.
      \param varPathStr IN: A variable path string in the file
      structure. (optional)
       If varPathStr is not given, it means any variable path.
      \param selectForm IN: indicate the coordinates are points or boxes (default: points)

      \return True if success, otherwise return False.
    */

    bool getSelectedData(const std::string &varNameStr,
                         const std::vector<uint64_t> &coords,
                         void *data,
                         const std::string &varPathStr="",
                         const FQ::SelectForm selectForm=FQ::POINTS_SELECTION);

    /*!
      \brief Get the max/min value of a variable.

      \param variable IN: The location of a variable in the file structure.
      The location can be givin by either the name or the full-path of the
      variable.  E.g. variable="a" refers to the variable locating at "/a"
      in the file.  variable="/time0/a" refers to the variable locating at
      "/time0/a".

      \return the max/min value.  The type double is used to represent all
      types of data value.  Users are responsible to cast the value into
      the correct type of the variable.  There is a potential prcision lose
      for unsigned int 64
    */
    //    double getDataMax(const std::string &variable);
    //    double getDataMin(const std::string &variable);


    /* Core FastBit Search & Indexing APIs */
/*!
  \brief Get the number of hits that satisfy a query condition.
        This must be a collective call from all MPI tasks.

  \param query IN: SQL-like query string.

  \param varPathStr IN: A variable path string in the file
  structure (optional).
  If varPathStr is not given, it means any variable path.
  The same varPathStr is applied to every varNameStr used in the query.

  \param mpi_dim IN: The dimension to split the data into subarray in
  the MPI mode.

  \param mpi_len IN: The subarray size for splitting the data in the
  MPI mode.

  \return the number of hits.
*/

    uint64_t getNumHits(const std::string &query,
                        const std::string &varPathStr="",
                        uint64_t mpi_dim=FQ_DEFAULT_MPI_DIM,
                        uint64_t mpi_len= FQ_DEFAULT_MPI_LEN);

/*!
  \brief Get the vector coordinates of values that satisfy a query
  condition.
        This must be a collective call from all MPI tasks.

  \param query IN: SQL-like query string.

  \param coords OUT: Vector containing the coordinates of the hit values
  that satisfy the query condition.

  \param varPathStr IN: A variable path string in the file structure
  (optional).
  If varPathStr is not given, it means any variable path.
  The same varPathStr is applied to every varNameStr used in the query.

  \param selectForm: The selection form of data selection.
  - If selectForm is POINTS_SELECTION, the data is retrieved by a vector
  of points.
  - If selectForm is LINES_SELECTION, the data is retrieved by a vector
  of bounding lines.
  - If selectForm is BOXES_SELECTION, the data is retrieved by a vector
  of bounding boxes.
  If selectForm is not specified, the default is POINTS_SELECTION.

  \param mpi_dim IN: The dimension to split the data into subarray in
  the MPI mode.

  \param mpi_len IN: The subarray size for splitting the data in the
  MPI mode.

  \param bcast IN: Broadcast data to all MPI tasks (default: true).

  \return the number of hit points/boxes.
*/

    uint64_t executeQuery(const std::string &query,
                          std::vector<uint64_t> &coords,
                          const std::string &varPathStr="",
                          const FQ::SelectForm selectForm=FQ::POINTS_SELECTION,
                          uint64_t mpi_dim=FQ_DEFAULT_MPI_DIM,
                          uint64_t mpi_len= FQ_DEFAULT_MPI_LEN,
                          const bool bcast=true);

/*!
  \brief Get the point coordinates of values that satisfy the equality
  condition form a variable.  The order of return coordinates may be
  different from the given order of identifier.
        This must be a collective call from all MPI tasks.

  \param varNameStr IN: A variable name string in the file structure.

  \param identifiers IN: Vector containing the equality selection values.
  Values is specified in double to avoid string creation/parsing overhead.

  \param coords OUT: Vector containing the coordinates of the data
  that satisfy the equality condition.

  \param varPathStr IN: A variable path string in the file structure
  (optional).
  If varPathStr is not given, it means any variable path.
  The same varPathStr is applied to every varNameStr used in the query.

  \param mpi_dim IN: The dimension to split the data into subarray in
  the MPI mode.

  \param mpi_len IN: The subarray size for splitting the data in the
  MPI mode.

  \param bcast IN: Broadcast data to all MPI tasks (default: true).

  \return the number of hits
*/

    uint64_t executeEqualitySelectionQuery
		    (const std::string &varNameStr,
		     const std::vector<double> &identifiers,
		     std::vector<uint64_t> &coords,
		     const std::string &varPathStr="",
		     const FQ::SelectForm selectForm=FQ::POINTS_SELECTION,
		     uint64_t mpi_dim=FQ_DEFAULT_MPI_DIM,
		     uint64_t mpi_len=FQ_DEFAULT_MPI_LEN,
                     const bool bcast=true);

    uint64_t get1DHistogram (const char *query,
                             const std::string& varNameStr,
                             const std::string& varPathStr,
                             double begin, double end, double stride,
                             std::vector<uint32_t> &counts,
                             uint64_t mpi_dim, uint64_t mpi_len);

    // Compute weighted conditional 1D histogram with regularly spaced bins.
    uint64_t get1DHistogram (const char *query,
                             const std::string& varNameStr,
                             const std::string& varPathStr,
                             double begin, double end, double stride,
                             const std::string& wtNameStr,
                             std::vector<double> &weights);

    // Compute 1D histogram with adaptive bins.
    uint64_t get1DHistogram (const std::string& varNameStr,
                             uint32_t nbin,
                             const std::string& varPathStr,
                             std::vector<double> &bounds,
                             std::vector<uint32_t> &counts);

    // Compute conditional 1D histogram with adaptive bins.
    uint64_t get1DHistogram (const char *query,
                             const std::string& varNameStr,
                             uint32_t nbin,
                             const std::string& varPathStr,
                             std::vector<double> &bounds,
                             std::vector<uint32_t> &counts);

    // Partition values of the named variable into regularly spaced bins.
    uint64_t get1DBins  (const char *query,
                         const std::string& varNameStr,
                         const std::string& varPathStr,
                         double begin, double end, double stride,
                         std::vector<ibis::bitvector> &bins);

    uint64_t get1DBins  (const char *query,
                         const std::string& varNameStr,
                         const std::string& varPathStr,
                         double begin, double end, double stride,
                         std::vector<ibis::bitvector *> &bins);

    uint64_t get1DBins  (const char *query,
                         const std::string& varNameStr,
                         const std::string& varPathStr,
                         double begin, double end, double stride,
                         const std::string& wtNameStr,
                         std::vector<double> &weights,
                         std::vector<ibis::bitvector *> &bins);

    uint64_t get1DBins  (const char *query,
                         const std::string& varNameStr,
                         const std::string& varPathStr,
                         uint32_t nbin,
                         double begin, double end, double stride,
                         std::vector<double> &bounds,
                         std::vector<ibis::bitvector> &bins);


    // 2D Histogram functions
    uint64_t get2DHistogram (const char *query,
                             const std::string& varPathStr,
                             const std::string& varNameStr1,
                             double begin1, double end1, double stride1,
                             const std::string& varNameStr2,
                             double begin2, double end2, double stride2,
                             std::vector<uint32_t> &counts,
                             uint64_t mpi_dim, uint64_t mpi_len);

    uint64_t get2DHistogram (const char *query,
                             const std::string& varPathStr,
                             const std::string& varNameStr1,
                             double begin1, double end1, double stride1,
                             const std::string& varNameStr2,
                             double begin2, double end2, double stride2,
                             const std::string& wtNameStr,
                             std::vector<double> &weights);

    // Compute 2D histogram with adaptive bins.
    uint64_t get2DHistogram (const std::string& varNameStr1,
                             const std::string& varNameStr2,
                             const std::string& varPathStr,
                             uint32_t nbin1,
                             uint32_t nbin2,
                             std::vector<double> &bounds1,
                             std::vector<double> &bounds2,
                             std::vector<uint32_t> &counts,
                             const char* const option);

    // Compute conditional 2D histogram with adaptive bins.
    uint64_t get2DHistogram (const char *query,
                             const std::string& varNameStr1,
                             const std::string& varNameStr2,
                             const std::string& varPathStr,
                             uint32_t nbin1,
                             uint32_t nbin2,
                             std::vector<double> &bounds1,
                             std::vector<double> &bounds2,
                             std::vector<uint32_t> &counts);

    uint64_t get2DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         double begin1, double end1, double stride1,
                         const std::string& varNameStr2,
                         double begin2, double end2, double stride2,
                         std::vector<ibis::bitvector> &bins);

    uint64_t get2DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         double begin1, double end1, double stride1,
                         const std::string& varNameStr2,
                         double begin2, double end2, double stride2,
                         std::vector<ibis::bitvector *> &bins);

    uint64_t get2DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         double begin1, double end1, double stride1,
                         const std::string& varNameStr2,
                         double begin2, double end2, double stride2,
                         const std::string& wtNameStr,
                         std::vector<double> &weights,
                         std::vector<ibis::bitvector *> &bins);

    uint64_t get2DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         const std::string& varNameStr2,
                         uint32_t nbin1,
                         uint32_t nbin2,
                         std::vector<double> &bounds1,
                         std::vector<double> &bounds2,
                         std::vector<ibis::bitvector> &bins);


    // 3D Histogram functions
    uint64_t get3DHistogram (const char *query,
                             const std::string& varPathStr,
                             const std::string& varNameStr1,
                             double begin1, double end1, double stride1,
                             const std::string& varNameStr2,
                             double begin2, double end2, double stride2,
                             const std::string& varNameStr3,
                             double begin3, double end3, double stride3,
                             std::vector<uint32_t> &counts,
                             uint64_t mpi_dim, uint64_t mpi_len);

    uint64_t get3DHistogram (const char *query,
                             const std::string& varPathStr,
                             const std::string& varNameStr1,
                             double begin1, double end1, double stride1,
                             const std::string& varNameStr2,
                             double begin2, double end2, double stride2,
                             const std::string& varNameStr3,
                             double begin3, double end3, double stride3,
                             const std::string& wtNameStr,
                             std::vector<double> &weights);

    uint64_t get3DHistogram (const std::string& varNameStr1,
                             const std::string& varNameStr2,
                             const std::string& varNameStr3,
                             const std::string& varPathStr,
                             uint32_t nbin1,
                             uint32_t nbin2,
                             uint32_t nbin3,
                             std::vector<double> &bounds1,
                             std::vector<double> &bounds2,
                             std::vector<double> &bounds3,
                             std::vector<uint32_t> &counts);

    uint64_t get3DHistogram (const char *query,
                             const std::string& varNameStr1,
                             const std::string& varNameStr2,
                             const std::string& varNameStr3,
                             const std::string& varPathStr,
                             uint32_t nbin1,
                             uint32_t nbin2,
                             uint32_t nbin3,
                             std::vector<double> &bounds1,
                             std::vector<double> &bounds2,
                             std::vector<double> &bounds3,
                             std::vector<uint32_t> &counts);

    uint64_t get3DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         double begin1, double end1, double stride1,
                         const std::string& varNameStr2,
                         double begin2, double end2, double stride2,
                         const std::string& varNameStr3,
                         double begin3, double end3, double stride3,
                         std::vector<ibis::bitvector> &bins);

    uint64_t get3DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         double begin1, double end1, double stride1,
                         const std::string& varNameStr2,
                         double begin2, double end2, double stride2,
                         const std::string& varNameStr3,
                         double begin3, double end3, double stride3,
                         std::vector<ibis::bitvector *> &bins);


    uint64_t get3DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         double begin1, double end1, double stride1,
                         const std::string& varNameStr2,
                         double begin2, double end2, double stride2,
                         const std::string& varNameStr3,
                         double begin3, double end3, double stride3,
                         const std::string& wtNameStr,
                         std::vector<double> &weights,
                         std::vector<ibis::bitvector *> &bins);


    uint64_t get3DBins  (const char *query,
                         const std::string& varPathStr,
                         const std::string& varNameStr1,
                         const std::string& varNameStr2,
                         const std::string& varNameStr3,
                         uint32_t nbin1,
                         uint32_t nbin2,
                         uint32_t nbin3,
                         std::vector<double> &bounds1,
                         std::vector<double> &bounds2,
                         std::vector<double> &bounds3,
                         std::vector<ibis::bitvector> &bins);


    int recordRegions(const std::string &outputfile,
                      const std::string &path,
                      const std::string &cond,
                      const std::string &col,
                      const std::vector<double> &thr) const;

};

#endif
