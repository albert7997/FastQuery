#ifndef _FQ_PART_H
#define _FQ_PART_H

#include <part.h>       // ibis::part
#include <meshQuery.h>  // ibis::meshQuery
#include <util.h>
#include <vector>
#include <string>
#include "arrayIODriver.h"
#include "fqVar.h"

/**
   The class FQ_Part is a container of the variables in data file.
   It mirrors a FastBit data partition.
*/
class FQ_Part : public ibis::part {
public:
    /*******************************************
     * Constructor and Destructor
     ********************************************/
    /*!
      \brief Initiate a FQ_Part object for a set of variables.

      \param variables IN: Vector containing a list of variables given by
      their full-paths.
      \param dataFile IN: A reference to a ArrayIODriver object which
      implements all file I/O APIs needed by the FQ_Part object.
    */
    FQ_Part(const std::vector<VarInfo> &varInfoLists,
            const std::vector<VarSpace> &varSpaceLists,
            const ArrayIODriver &dataFile,
            const ArrayIODriver &indexFile,
            const int mpi_rank=-1);

    virtual ~FQ_Part() {releaseAllQueries();}

    /*******************************************
     * Core APIs for FastbitHelper
     ********************************************/
    bool isValid (const std::string &func) const;

    void destroyQuery(const char *token);
    const char* createQuery(const std::string &cond);
    const char* createEqualitySelectionQuery
    (const std::string &variable, const std::vector<double> &identifiers);

    /// Fully evaluate the query. This will return the exact number of
    /// hits.
    int64_t submitQuery(const char *token);
    const ibis::meshQuery* getQuery(const char*) const;

    int64_t getNumHits(const char *token) const;
    int64_t getHitsAsPoints(const char *token,
                            std::vector<uint64_t>&) const;
    int64_t getHitsAsLines(const char *token,
                           std::vector<uint64_t>&) const;
    int64_t getHitsAsBoxes(const char *token,
                           std::vector<uint64_t>&) const;

    using ibis::part::buildIndexes;
    int buildIndexes(ArrayIODriver &indexFile, const char* opt=0);

    /// Return an reference to an array storing the extend of the
    /// dimensions of the regular mesh that defines the data sets.
    const std::vector<uint32_t>& getMeshDims() const { return shapeSize;}

    //Added by Suren
    // get conditional 1D histogram with regularly spaced bins.
    int get1DHist_counts(const std::string &cond, const std::string &varName,
                         double begin, double end, double stride,
                         std::vector<uint32_t> &counts);

    // get weighted conditional 1D histogram with regularly spaced bins.
    int get1DHist_weights(const std::string &cond, const std::string &varName,
                          double begin, double end, double stride,
                          const std::string &wtName,
                          std::vector<double> &weights);

    // get 1D histogram with adaptive bins
    int get1DAdaptiveHist_counts(const std::string &varName,
                                 uint32_t nbin,
                                 std::vector<double> &bounds,
                                 std::vector<uint32_t> &counts);

    // get 1D conditional histogram with adaptive bins
    int get1DAdaptiveHist_cond_counts(const std::string &cond,
                                      const std::string &varName,
                                      uint32_t nbin,
                                      std::vector<double> &bounds,
                                      std::vector<uint32_t> &counts);

    // get partition values of the named variable into regularly spaced bins.
    int get1DBins(const std::string &cond, const std::string &varName,
                  double begin, double end, double stride,
                  std::vector<ibis::bitvector> &bins);

    // get partition values of the named variable into regularly spaced bins
    int get1DBins(const std::string &cond, const std::string &varName,
                  double begin, double end, double stride,
                  std::vector<ibis::bitvector *> &bins);

    // get partition values of the named variable into regularly spaced bins.
    int get1DBins_weights(const std::string &cond, const std::string &varName,
                          double begin, double end, double stride,
                          const std::string &wtName,
                          std::vector<double> &weights,
                          std::vector<ibis::bitvector *> &bins);


    //2D Histogram functions
    // get conditional 2D histogram with regularly spaced bins.
    int get2DHist_counts(const std::string &cond,
                         const std::string &varName1,
                         double begin1, double end1, double stride1,
                         const std::string &varName2,
                         double begin2, double end2, double stride2,
                         std::vector<uint32_t> &counts);

    // get weighted conditional 2D histogram with regularly spaced bins.
    int get2DHist_weights(const std::string &cond,
                          const std::string &varName1,
                          double begin1, double end1, double stride1,
                          const std::string &varName2,
                          double begin2, double end2, double stride2,
                          const std::string &wtName,
                          std::vector<double> &weights);

    // get 2D histogram with adaptive bins - Without condition
    int get2DAdaptiveHist_counts(const std::string &varName1,
                                 const std::string &varName2,
                                 uint32_t num_bins1,
                                 uint32_t num_bins2,
                                 std::vector<double> &bounds1,
                                 std::vector<double> &bounds2,
                                 std::vector<uint32_t> &counts);

    // get 2D histogram with adaptive bins - With condition
    int get2DAdaptiveHist_cond_counts(const std::string &cond,
                                      const std::string &varName1,
                                      const std::string &varName2,
                                      uint32_t num_bins1,
                                      uint32_t num_bins2,
                                      std::vector<double> &bounds1,
                                      std::vector<double> &bounds2,
                                      std::vector<uint32_t> &counts);

    // get partition values of the named variable into regularly spaced bins.
    int get2DBins(const std::string &cond,
                  const std::string &varName1,
                  double begin1, double end1, double stride1,
                  const std::string &varName2,
                  double begin2, double end2, double stride2,
                  std::vector<ibis::bitvector> &bins);

    // get partition values of the named variable into regularly spaced bins
    int get2DBins(const std::string &cond,
                  const std::string &varName1,
                  double begin1, double end1, double stride1,
                  const std::string &varName2,
                  double begin2, double end2, double stride2,
                  std::vector<ibis::bitvector *> &bins);

    // get partition values of the named variable into regularly spaced bins.
    int get2DBins_weights(const std::string &cond,
                          const std::string &varName1,
                          double begin1, double end1, double stride1,
                          const std::string &varName2,
                          double begin2, double end2, double stride2,
                          const std::string &wtName,
                          std::vector<double> &weights,
                          std::vector<ibis::bitvector *> &bins);

    //3D Histogram functions
    // get conditional 2D histogram with regularly spaced bins.
    int get3DHist_counts(const std::string &cond,
                         const std::string &varName1,
                         double begin1, double end1, double stride1,
                         const std::string &varName2,
                         double begin2, double end2, double stride2,
                         const std::string &varName3,
                         double begin3, double end3, double stride3,
                         std::vector<uint32_t> &counts);

    // get weighted conditional 2D histogram with regularly spaced bins.
    int get3DHist_weights(const std::string &cond,
                          const std::string &varName1,
                          double begin1, double end1, double stride1,
                          const std::string &varName2,
                          double begin2, double end2, double stride2,
                          const std::string &varName3,
                          double begin3, double end3, double stride3,
                          const std::string &wtName,
                          std::vector<double> &weights);

    //end additions by Suren

    /*******************************************
     * Internal APIs used by Fastbit
     ********************************************/
    using ibis::part::doScan;
    // ---- required by FastBit to perform candidate check ----
    /// Evaluate query by reading the values from the data file.
    virtual long doScan(const ibis::qRange& cmp,
                        const ibis::bitvector& mask,
                        ibis::bitvector& hits) const;

private:
    // member variables
    typedef std::map< const char*, ibis::meshQuery*,
                      std::less< const char* > > queryList;
    queryList qlist; // list of querys on this set of variables
    bool valid;
    int mpi_rank;

    FQ_Part(); // no default constructor
    FQ_Part(const FQ_Part&); // no copy constructor
    const FQ_Part& operator=(const FQ_Part&); // no assignment
    void releaseAllQueries();
}; // class FQ_Part
#endif
