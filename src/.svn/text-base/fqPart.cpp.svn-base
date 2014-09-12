#include "fqPart.h"
#include "fqIndex.h"
#include "fq.h"         // FastQuery::reportTiming()

/*******************************************
 * Constructor and Destructor
 ********************************************/
FQ_Part::FQ_Part(const std::vector<VarInfo> &varInfoList,
                 const std::vector<VarSpace> &varSpaceList,
                 const ArrayIODriver &dataFile,
                 const ArrayIODriver &indexFile, const int rank)
    : ibis::part(static_cast<const char*>(0), static_cast<const char*>(0)) {
    // give a name to the partition object
    std::ostringstream ostr;
    ostr << dataFile.getFileName();
    m_name = ibis::util::strnewdup(ostr.str().c_str());
    mpi_rank = rank;

    if (varSpaceList.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_Part: must have at least one variable to "
            "create a partition";
        valid = false;
        return;
    }

    // initialize column object for each variable
    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part: initializing variable partition";

    std::vector<uint64_t> dims = varSpaceList[0].getCounts();
    for (unsigned int i = 0; i < varInfoList.size(); ++ i) {
        std::ostringstream varFBName("");
        varFBName << "Var";
        varFBName << i;
        varFBName << "\0";
        std::string varPath = varInfoList[i].getPath();
        LOGGER(ibis::gVerbose > 2)
            << "FQ_Part: new column for variable "
            << varSpaceList[i].getText().c_str()
            << " as " << varFBName.str().c_str();
        /*
          std::vector<uint64_t> newDims = varSpaceList[i].getCounts();
          if (newDims.size() != dims.size()) {
          LOGGER(ibis::gVerbose > 0)
          << "Warning -- FQ_Part: variables have different number "
          "of dimensions";
          valid = false;
          return;
          }
          for (unsigned int j=0; j<dims.size(); j++) {
          if (newDims[j] != dims[j]) {
          LOGGER(ibis::gVerbose > 0)
          << "Warning -- FQ_Part variables must have matching"
          " dimensions";
          valid = false;
          return;
          }
          }
        */
        FQ_Variable *var =
            new FQ_Variable(this, varFBName.str(), varPath,
                            varInfoList[i], varSpaceList[i],
                            dataFile, indexFile);
        columns[var->name()] = var;
    }
    // set dimension size
    shapeSize.resize(dims.size());
    nEvents = 1;
    for (unsigned j = 0; j < dims.size(); ++ j) {
        shapeSize[j] = static_cast<unsigned>(dims[j]);
        nEvents *= dims[j];
    }
    if (nEvents > 0 && ! columns.empty()) {// define a name for the table
        amask.set(1, nEvents); // every record is valid

        // TODO: need to allocate a cleaner object, actually define a new
        // cleaner class!
    }
    valid = true;
} // FQ_Part::FQ_Part

/*******************************************
 * Public Functions
 ********************************************/
/// Create a query object for the given set of query conditions.
/// All query conditions are specified in a single string.
const char *FQ_Part::createQuery(const std::string &cond) {
    const char *token = 0;
    if (cond.empty()) return token;
    if (! isValid("FQ_Part::createQuery")) return token;
    { // is the query expression already in the list of queries
        readLock lock(this, "createQuery");
        for (queryList::const_iterator it = qlist.begin();
             it != qlist.end(); ++ it) {
            if (cond.compare(it->second->getWhereClause()) == 0)
                return it->first;
        }
    }

    ibis::meshQuery *sel = new ibis::meshQuery(ibis::util::userName(), this);
    sel->setWhereClause(cond.c_str());
    //sel->evaluate(); // evaluate the query to determine hits
    token = sel->id();
    { // need a mutex lock to modify the query list
        mutexLock lock(this, "createQuery");
        qlist[token] = sel;
        LOGGER(ibis::gVerbose > 2)
            << "FQ_Part::createQuery: has " << qlist.size()
            << " quer" << (qlist.size() > 1 ? "ies" : "y") << " in memory";
    }
    return token;
} // FQ_Part::createQuery

const char *
FQ_Part::createEqualitySelectionQuery(const std::string &variable,
                                      const std::vector<double>& identifiers) {
    if (! isValid("FQ_Part::createEqualitySelectionQuery")) return "";

    const char *token = 0;
    ibis::meshQuery *sel = new ibis::meshQuery(ibis::util::userName(), this);
    ibis::qDiscreteRange* dRange =
        new ibis::qDiscreteRange(variable.c_str(), identifiers);
    sel->setWhereClause(dRange);

    token = sel->id();
    if (token != 0) {
        mutexLock lock(this, "createQuery");
        qlist[token] = sel;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::createEqualitySelectionQuery:"
        << " has " << qlist.size()
        << " quer" << (qlist.size() > 1 ? "ies" : "y") << " in memory";
    return token;
} // FQ_Part::createEqualitySelectionQuery

int64_t FQ_Part::submitQuery(const char *token) {
    if (! isValid("FQ_Part::submitQuery")) return 0;

    ibis::horometer timer;
    timer.start();
    int64_t nhits = 0;
    readLock lock(this, "submitQuery");
    queryList::const_iterator it = qlist.find(token);
    if (it != qlist.end()) {
        ibis::query::QUERY_STATE qst = (*it).second->getState();
        timer.stop();
        switch (qst) {
        default:
            logWarning("submitQuery", "query not fully specified");
            break;
        case ibis::query::QUICK_ESTIMATE:
        case ibis::query::SET_PREDICATE:
        case ibis::query::SPECIFIED:
        case ibis::query::SET_RIDS:
            (*it).second->evaluate();
        case ibis::query::FULL_EVALUATE:
            nhits = (*it).second->getNumHits();
            break;
        }
    }
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::submitQuery\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;
    return nhits;
} // FQ_Part::submitQuery

//Added by Suren
// Compute 1D histogram counts
int FQ_Part::get1DHist_counts(const std::string &cond,
                              const std::string &varName,
                              double begin, double end, double stride,
                              std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get1DHist_counts")) return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DHist_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DHist_counts: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DDistribution(cond.c_str(), varName.c_str (),
                                  begin, end, stride, counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DHist_counts\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get1DHist_counts

//Compute weighted 1D Histogram
int FQ_Part::get1DHist_weights(const std::string &cond,
                               const std::string &varName,
                               double begin, double end, double stride,
                               const std::string &wtName,
                               std::vector<double> &weights) {
    if (! isValid("FQ_Part::get1DHist_weights")) return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DHist_weights");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DHist_weights: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DDistribution(cond.c_str(), varName.c_str (),
                                  begin, end, stride, wtName.c_str(), weights);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DHist_weights\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
}  //FQ_Part::get1DHist_weights

// get 1D histogram with adaptive bins
int FQ_Part::get1DAdaptiveHist_counts(const std::string &varName,
                                      uint32_t nbin,
                                      std::vector<double> &bounds,
                                      std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get1DAdaptiveHist_counts"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DAdaptiveHist_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DHist_weights: Condition string: \t"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DDistribution(varName.c_str (), nbin, bounds, counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DAdaptiveHist_counts\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;

} // FQ_Part::get1DAdaptiveHist_counts

// get 1D conditional histogram with adaptive bins
int FQ_Part::get1DAdaptiveHist_cond_counts(const std::string &cond,
                                           const std::string &varName,
                                           uint32_t nbin,
                                           std::vector<double> &bounds,
                                           std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get1DHist_cond_counts"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DAdaptiveHist_cond_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DHist_weights: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DDistribution(cond.c_str(), varName.c_str (), nbin,
                                  bounds, counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DAdaptiveHist_cond_counts\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;

} // FQ_Part::get1DAdaptiveHist_cond_counts

// Compute 1D histogram bins
int FQ_Part::get1DBins (const std::string &cond, const std::string &varName,
                        double begin, double end, double stride,
                        std::vector<ibis::bitvector> &bins) {
    if (! isValid("FQ_Part::get1DBins"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DBins");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DBins: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DBins (cond.c_str(), varName.c_str (),
                           begin, end, stride, bins);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DBins\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get1DBins

// Compute 1D histogram bins --> passing <ibis::bitvector *> &bins
int FQ_Part::get1DBins (const std::string &cond, const std::string &varName,
                        double begin, double end, double stride,
                        std::vector<ibis::bitvector *> &bins) {
    if (! isValid("FQ_Part::get1DBins"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DBins");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DBins: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DBins(cond.c_str(), varName.c_str (),
                          begin, end, stride, bins);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DBins\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get1DBins

// Compute 1D histogram bins with weights
int FQ_Part::get1DBins_weights (const std::string &cond,
                                const std::string &varName,
                                double begin, double end, double stride,
                                const std::string &wtName,
                                std::vector<double> &weights,
                                std::vector<ibis::bitvector *> &bins) {
    if (! isValid("FQ_Part::get1DBins_weights"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get1DBins_weights");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get1DHist_bins: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name: " << varName.c_str ();

    ibis::part::get1DBins(cond.c_str(), varName.c_str (),
                          begin, end, stride, wtName.c_str (), weights, bins);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get1DBins\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get1DBins

// 2D Histogram function definitions
int FQ_Part::get2DHist_counts (const std::string &cond,
                               const std::string &varName1,
                               double begin1, double end1, double stride1,
                               const std::string &varName2,
                               double begin2, double end2, double stride2,
                               std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get2DHist_counts"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DHist_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DHist_counts: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\nVariable name 2: " << varName2.c_str () ;

    ibis::part::get2DDistribution(cond.c_str(),
                                  varName1.c_str (), begin1, end1, stride1,
                                  varName2.c_str (), begin2, end2, stride2,
                                  counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DHist_counts and counts.size\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t"
        << counts.size() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get2DHist_counts

//Compute weighted 2D Histogram
int FQ_Part::get2DHist_weights(const std::string &cond,
                               const std::string &varName1,
                               double begin1, double end1, double stride1,
                               const std::string &varName2,
                               double begin2, double end2, double stride2,
                               const std::string &wtName,
                               std::vector<double> &weights) {
    if (! isValid("FQ_Part::get2DHist_weights"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DHist_weights");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DHist_weights: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\nVariable name 2: " << varName2.c_str ();

    ibis::part::get2DDistribution(cond.c_str(),
                                  varName1.c_str (), begin1, end1, stride1,
                                  varName2.c_str (), begin2, end2, stride2,
                                  wtName.c_str(), weights);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DHist_weights\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
}  //FQ_Part::get2DHist_weights

int FQ_Part::get2DAdaptiveHist_counts(const std::string &varName1,
                                      const std::string &varName2,
                                      uint32_t num_bins1,
                                      uint32_t num_bins2,
                                      std::vector<double> &bounds1,
                                      std::vector<double> &bounds2,
                                      std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get2DAdaptiveHist_counts"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DAdaptiveHist_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DHist_counts: Condition string: \t"
        << "Variable name 1: " << varName1.c_str ()
        << "\nVariable name 2: " << varName2.c_str () ;

    ibis::part::get2DDistribution(varName1.c_str (), varName2.c_str (),
                                  num_bins1, num_bins2,
                                  bounds1, bounds2,
                                  counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DAdaptiveHist_counts\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get2DAdaptiveHist_counts

int FQ_Part::get2DAdaptiveHist_cond_counts(const std::string &cond,
                                           const std::string &varName1,
                                           const std::string &varName2,
                                           uint32_t num_bins1,
                                           uint32_t num_bins2,
                                           std::vector<double> &bounds1,
                                           std::vector<double> &bounds2,
                                           std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get2DAdaptiveHist_cond_counts"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DAdaptiveHist_cond_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DHist_counts: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\nVariable name 2: " << varName2.c_str () ;

    ibis::part::get2DDistribution(cond.c_str (),
                                  varName1.c_str (), varName2.c_str (),
                                  num_bins1, num_bins2,
                                  bounds1, bounds2,
                                  counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DAdaptiveHist_cond_counts\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get2DAdaptiveHist_cond_counts

// get2DBins
int FQ_Part::get2DBins (const std::string &cond,
                        const std::string &varName1,
                        double begin1, double end1, double stride1,
                        const std::string &varName2,
                        double begin2, double end2, double stride2,
                        std::vector<ibis::bitvector> &bins) {
    if (! isValid("FQ_Part::get2DBins"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DBins");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DBins: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\tVariable name 2: " << varName2.c_str () ;

    ibis::part::get2DBins(cond.c_str(),
                          varName1.c_str (), begin1, end1, stride1,
                          varName2.c_str (), begin2, end2, stride2,
                          bins);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DBins \t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get2DBins

// get2DBins ->> <ibis::bitvector *> &bins
int FQ_Part::get2DBins(const std::string &cond,
                       const std::string &varName1,
                       double begin1, double end1, double stride1,
                       const std::string &varName2,
                       double begin2, double end2, double stride2,
                       std::vector<ibis::bitvector *> &bins) {
    if (! isValid("FQ_Part::get2DBins")) {
        return 0;
    }

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DBins");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DBins: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\tVariable name 2: " << varName2.c_str () ;

    ibis::part::get2DBins(cond.c_str(),
                          varName1.c_str (), begin1, end1, stride1,
                          varName2.c_str (), begin2, end2, stride2,
                          bins);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DBins \t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get2DBins

// get2DBins_weights
int FQ_Part::get2DBins_weights(const std::string &cond,
                               const std::string &varName1,
                               double begin1, double end1, double stride1,
                               const std::string &varName2,
                               double begin2, double end2, double stride2,
                               const std::string &wtName,
                               std::vector<double> &weights,
                               std::vector<ibis::bitvector *> &bins) {
    if (! isValid("FQ_Part::get2DBins_weights"))  {
        return 0;
    }

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get2DBins_weights");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get2DBins_weights: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\tVariable name 2: " << varName2.c_str () ;

    ibis::part::get2DBins(cond.c_str(),
                          varName1.c_str (), begin1, end1, stride1,
                          varName2.c_str (), begin2, end2, stride2,
                          wtName.c_str(), weights,
                          bins);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get2DBins_weights \t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get2DBins


// 3D Histogram function definitions
int FQ_Part::get3DHist_counts(const std::string &cond,
                              const std::string &varName1,
                              double begin1, double end1, double stride1,
                              const std::string &varName2,
                              double begin2, double end2, double stride2,
                              const std::string &varName3,
                              double begin3, double end3, double stride3,
                              std::vector<uint32_t> &counts) {
    if (! isValid("FQ_Part::get3DHist_counts")) {
        return 0;
    }

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get3DHist_counts");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get3DHist_counts: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\nVariable name 2: " << varName2.c_str ()
        << "\nVariable name 3: " << varName3.c_str () ;

    ibis::part::get3DDistribution(cond.c_str(),
                                  varName1.c_str (), begin1, end1, stride1,
                                  varName2.c_str (), begin2, end2, stride2,
                                  varName3.c_str (), begin3, end3, stride3,
                                  counts);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get3DHist_counts\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
} //FQ_Part::get3DHist_counts

//Compute weighted 3D Histogram
int FQ_Part::get3DHist_weights(const std::string &cond,
                               const std::string &varName1,
                               double begin1, double end1, double stride1,
                               const std::string &varName2,
                               double begin2, double end2, double stride2,
                               const std::string &varName3,
                               double begin3, double end3, double stride3,
                               const std::string &wtName,
                               std::vector<double> &weights) {
    if (! isValid("FQ_Part::get3DHist_weights"))
        return 0;

    ibis::horometer timer;
    timer.start();
    readLock lock (this, "get3DHist_weights");

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::get3DHist_weights: Condition string: \t"
        << cond.c_str () << "\n"
        << "Variable name 1: " << varName1.c_str ()
        << "\nVariable name 2: " << varName2.c_str ()
        << "\nVariable name 3: " << varName3.c_str () ;

    ibis::part::get3DDistribution(cond.c_str(),
                                  varName1.c_str (), begin1, end1, stride1,
                                  varName2.c_str (), begin2, end2, stride2,
                                  varName3.c_str (), begin3, end3, stride3,
                                  wtName.c_str(), weights);

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::get3DHist_weights\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;

    return 1;
}  //FQ_Part::get3DHist_weights

// Suren: End histogram functions

/// Return the query object associated with the given query token.
const ibis::meshQuery* FQ_Part::getQuery(const char *tok) const {
    queryList::const_iterator it = qlist.find(tok);
    if (it != qlist.end()) {
        return it->second;
    }
    return 0;
} // FQ_Part::getQuery

/// Return the number of records selected.  It will return -1 before
/// the function submitQuery is called.
int64_t FQ_Part::getNumHits(const char *token) const {
    if (! isValid("FQ_Part::getNumHits")) return 0;

    int64_t nhits = 0;
    readLock lock(this, "getNumHits");
    queryList::const_iterator it = qlist.find(token);
    if (it != qlist.end() &&
        ibis::query::FULL_EVALUATE == (*it).second->getState())
        nhits = (*it).second->getNumHits();
    return nhits;
} // FQ_Part::getNumHits

/// Destroy the query specified by the token.  Should be done to
/// release internal resources if the query objection is no longer
/// needed.  Once this function is called, the query token can
/// no longer be used because the token was part of the query
/// object this function destroys.
void FQ_Part::destroyQuery(const char *token) {
    if (! isValid("FQ_Part::destroyQuery")) return;

    // need an exclusive lock to prevent variable it being modified between
    // the read and write lock.
    mutexLock mtx(this, "destroyQuery");
    queryList::iterator it = qlist.find(token);
    if (it != qlist.end()) {
        delete (*it).second;
        qlist.erase(it);
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::destoryQuery: " << qlist.size()
        << " quer" << (qlist.size() > 1 ? "ies" : "y") << " in memory";
} // FQ_Part::destroyQuery

/// Remove all stored selections.
void FQ_Part::releaseAllQueries() {
    if (! isValid("FQ_Part::releaseAllQueries")) return;

    // need to an exclusive lock to avoid interfering with destroyQuery
    mutexLock mtx(this, "releaseAllQueries");
    const unsigned nq = qlist.size();
    for (queryList::iterator it = qlist.begin(); it != qlist.end(); ++ it) {
        delete (*it).second;
    }
    qlist.clear();

    LOGGER(ibis::gVerbose > 2)
        << "FQ_Part::releaseAllQueries: removed " << nq
        << " quer" << (nq>1?"ies":"y")
        << ".\n\tBytes under FastBit management: "
        << ibis::fileManager::bytesInUse();
} // FQ_Part::releaseAllQueries

/// Build indexes for all variables in the data partition.  If binning is a
/// nil pointer, it attempts to build the binned indexes for the
/// floating-point valued columns and unbinned indexes for the integer
/// valued columns.
int FQ_Part::buildIndexes(ArrayIODriver &indexFile, const char* binning) {
    int ret = 0;
    int cnt = 0;
    if (! isValid("FQ_Part::buildIndexes")) return ret;
    if (nRows() < 1024) {
        LOGGER(ibis::gVerbose > 2)
            << "FQ_Part[" << name() << "]::buildIndexes will skip partitions "
            "with less than 1024 rows";
        return ret;
    }

    bool binned = false;
    if (binning != 0 && *binning != 0) {
        binned = (strncmp(binning, "<binning null", 13) != 0 &&
                  strncmp(binning, "<binning none", 13) != 0 &&
                  strncmp(binning, "<binning no ", 12) != 0 &&
                  strncmp(binning, "unbinned", 8) != 0);
    }

    ibis::horometer timer;
    timer.start();
    bool readOnly = false;
    for (ibis::part::columnList::const_iterator it = columns.begin();
         it != columns.end(); ++ it) {
        // call the constructor to either read an index or build a new one
        ibis::horometer splitTimer;
        splitTimer.start();
        const FQ_Variable* variable =
            reinterpret_cast<const FQ_Variable*>((*it).second);
        std::string varName = variable->getVarInfo().getPath();

        if (binned ||
            ((binning == 0 || *binning == 0) && variable->isFloat())) {
            FQ_IndexBinned tmp(variable, binning, readOnly);
            if (tmp.getStatus() < 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_Part::buildIndexes:"
                    << " failed to build the binned index for variable \""
                    << varName.c_str() << "\"" << std::endl;
            }
            else if (tmp.getStatus() > 0) { // write file
                if (tmp.write(indexFile) == false) {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- FQ_Part::buildIndexes:"
                        << " failed to write the binned index for variable \""
                        << varName.c_str() << "\"" << std::endl;
                }
                else {
                    LOGGER(ibis::gVerbose > 2)
                        << "FQ_Part::buildIndexes:"
                        << " successfully built and wrote the binned index "
                        "for variable \""
                        << varName.c_str() << "\"" << std::endl;
                    ret++;
                }
            }
            else {
                LOGGER(ibis::gVerbose > 2)
                    << "FQ_Part::buildIndexes:"
                    << " successfully read the binned index for variable \""
                    << varName.c_str() << "\"" << std::endl;
                ret++;
                cnt++;
            }
        }
        else {
            FQ_IndexUnbinned tmp(variable, readOnly);
            if (tmp.getStatus() < 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_Part::buildIndexes:"
                    << " failed to build the unbinned index for variable \""
                    << varName.c_str() << "\"" << std::endl;
            }
            else if (tmp.getStatus() > 0) { // write file
                if (tmp.write(indexFile) == false) {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- FQ_Part::buildIndexes:"
                        << " failed to write the unbinned index for variable \""
                        << varName.c_str() << "\"" << std::endl;
                }
                else {
                    LOGGER(ibis::gVerbose > 2)
                        << "FQ_Part::buildIndexes:"
                        << " successfully built and wrote the unbinned index "
                        "for variable \""
                        << varName.c_str() << "\"" << std::endl;
                    ret++;
                }
            }
            else {
                LOGGER(ibis::gVerbose > 2)
                    << "FQ_Part::buildIndexes:"
                    << " successfully read the binned unindex for variable \""
                    << varName.c_str() << "\"" << std::endl;
                ret++;
                cnt++;
            }
        }
    }
    if (ibis::gVerbose > 2) {
        timer.stop();
        logMessage("buildIndexes", "building %d index%s and reading %d index%s "
                   "used %G sec CPU time and %G sec elapsed time",
                   ret-cnt, ((ret-cnt)>1?"es":""), cnt, (cnt>1?"es":""),
                   timer.CPUTime(), timer.realTime());
    }
    return ret;
} // FQ_Part::buildIndexes

/// Retrieve the coordinates of the records selected by the query.  The
/// variable @c coord contains a list of bounding boxes for all mesh points
/// that satisfy the specified query conditions.  Let @c d be the number of
/// dimension of the mesh, the size of @c coord is expected to be @code 2 *
/// d * number-of-hits. @endcode Each bounding box is specified by the
/// position of the lower-left corner and the upper-right corder.  This
/// function return the actual number of bounding boxes it placed in @c
/// coord.  It returns a negative number if the @c token is not a known
/// token.  Make sure to call @c submitQuery before trying to call this
/// function.
int64_t FQ_Part::getHitsAsBoxes(const char *token,
                                std::vector<uint64_t>& coord) const {
    if (! isValid("FQ_Part::getHitsAsBoxes")) return 0;

    int64_t nhits = 0;
    readLock lock(this, "getHitsAsBoxes");
    queryList::const_iterator it = qlist.find(token);
    if (it != qlist.end()) {
        std::vector< std::vector<unsigned> > boxes;
        (void) (*it).second->getHitsAsBlocks(boxes);
        const unsigned int nDims = getMeshDims().size();
        coord.reserve(boxes.size()*(nDims+nDims));
        for (unsigned i = 0; i < boxes.size(); i++) {
            std::vector<unsigned> box = boxes[i];
            if (box.size() == nDims+nDims) {
                for (unsigned j = 0; j < nDims+nDims; j++) {
                    coord.push_back(box[j]);
                }
                ++ nhits;
            }
            else {
                LOGGER(ibis::gVerbose > 2)
                    << "Warning -- FQ_Part[" << name() << "]::getHitsAsBoxes("
                    << token << ") encountered bounding box " << i << " with "
                    << box.size() << " elements, but " << nDims+nDims
                    << " were expected";
            }
        }
    }
    return nhits;
} // FQ_Part::getHitsAsBoxes

/// Retrieve the coordinates of the records selected by the query.  The
/// variable @c coord contains a list of line segments for all mesh points
/// that satisfy the specified query conditions.  Let @c d be the number of
/// dimension of the mesh, the size of @c coord is expected to be @code
/// (d+1) * number-of-hits. @endcode The slow varying d-1 dimensions of a
/// line segment are specified by their coordiantes while the fastest
/// varying dimension is specified by a pair of (begin, end) pair.  This
/// function return the actual number of bounding boxes it placed in @c
/// coord.  It returns a negative number if the @c token is not a known
/// token.  Make sure to call @c submitQuery before trying to call this
/// function.
int64_t FQ_Part::getHitsAsLines(const char *token,
                                std::vector<uint64_t>& coord) const {
    if (! isValid("FQ_Part::getHitsAsLines")) return 0;

    int64_t nhits = 0;
    readLock lock(this, "getHitsAsLines");
    queryList::const_iterator it = qlist.find(token);
    if (it != qlist.end()) {
        std::vector<unsigned> lines;
        nhits = (*it).second->getHitsAsLines(lines);
        coord.reserve(lines.size());
        for (unsigned i = 0; i < lines.size(); i++) {
            coord.push_back(lines[i]);
        }
    }
    return nhits;
} // FQ_Part::getHitsAsLines

/// Retrieve the coordinates of the records selected by the query.  The
/// variable @c coord contains a list of coordinates of each mesh point
/// that satisfies the specified query conditions.  Let @c d be the
/// number of dimension of the mesh, the size of @c coord is expected
/// to be @c d * number-of-hits.  The coordinates of each mesh point is
/// given in the same order as return by the function @c getMeshDims.
/// This function return the actual number of point it placed in @c
/// coord.  It returns a negative number if the @c token is not a known
/// token.  Make sure to call @c submitQuery before trying to call this
/// function.
int64_t FQ_Part::getHitsAsPoints(const char *token,
                                 std::vector<uint64_t>& coord) const {
    if (! isValid("FQ_Part::getHitsAsPoints")) return 0;

    int64_t nhits = 0;
    readLock lock(this, "getHitsAsPoints");
    queryList::const_iterator it = qlist.find(token);
    if (it != qlist.end()) {
        // need to copy the coordinates
        std::vector<uint32_t> tmp;
        nhits = ibis::meshQuery::bitvectorToCoordinates
            (*((*it).second->getHitVector()), getMeshDims(), tmp);
        for (unsigned i = 0; i < tmp.size(); ++ i)
            coord.push_back(tmp[i]);
    }
    return nhits;
} // FQ_Part::getHitsAsPoints

// This implementation overrides the default one that relys on the
// vertically partitioned data files.
long FQ_Part::doScan(const ibis::qRange& cmp,
                     const ibis::bitvector& mask,
                     ibis::bitvector& hits) const {
    if (! isValid("FQ_Part::FQ_Part::doScan")) return 0;

    ibis::horometer readTimer;
    ibis::horometer compareTimer;

    LOGGER(ibis::gVerbose >= 5)
        << "... entering  to FQ_Part::doScan to resolve " << cmp;
    std::string colname(cmp.colName());
    columnList::const_iterator it = columns.find(cmp.colName());
    if (it == columns.end()) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_Part::doScan:"
            << " unable to find named column " << cmp.colName()
            << " in the data table";
        return 0;//-1;
    }
    int ierr = 0;
    const float mulfactor = 2;
    std::string variablePath =
        reinterpret_cast<const FQ_Variable*>((*it).second)->description();

    switch ((*it).second->type()) {
        //case ibis::column::KEY:
    case ibis::UINT: {
        readTimer.start();
        ibis::array_t<uint32_t> arr;
        if (mask.cnt() * mulfactor >= mask.size()) {
            // get all values
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getValuesArray(&arr);
        }
        else { // get only the selected values
            std::vector<uint32_t> tmp;
            std::vector<uint64_t> coords;
            ibis::meshQuery::bitvectorToCoordinates(mask, getMeshShape(),
                                                    tmp);
            coords.reserve(tmp.size());
            for (unsigned i = 0; i < tmp.size(); ++ i)
                coords.push_back(tmp[i]);
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getPointValues(arr, coords);
        }
        readTimer.stop();
        compareTimer.start();
        if (ierr >= 0)
            ierr = doCompare(arr, cmp, mask, hits);
        else
            hits.set(0, mask.size());
        break;
    }
    case ibis::INT: {
        readTimer.start();
        ibis::array_t<int32_t> arr;
        if (mask.cnt() * mulfactor >= mask.size()) {
            // get all values
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getValuesArray(&arr);
        }
        else { // get only the selected values
            std::vector<uint32_t> tmp;
            std::vector<uint64_t> coords;
            ibis::meshQuery::bitvectorToCoordinates(mask, getMeshShape(),
                                                    tmp);
            coords.reserve(tmp.size());
            for (unsigned i = 0; i < tmp.size(); ++ i)
                coords.push_back(tmp[i]);
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getPointValues(arr, coords);
        }
        readTimer.stop();
        compareTimer.start();
        if (ierr >= 0)
            ierr = doCompare(arr, cmp, mask, hits);
        else
            hits.set(0, mask.size());
        break;
    }
    case ibis::LONG: {
        readTimer.start();
        ibis::array_t<uint64_t> arr;
        if (mask.cnt() * mulfactor >= mask.size()) {
            // get all values
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getValuesArray(&arr);
        }
        else { // get only the selected values
            std::vector<uint32_t> tmp;
            std::vector<uint64_t> coords;
            ibis::meshQuery::bitvectorToCoordinates(mask, getMeshShape(),
                                                    tmp);
            coords.reserve(tmp.size());
            for (unsigned i = 0; i < tmp.size(); ++ i)
                coords.push_back(tmp[i]);
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getPointValues(arr, coords);
        }
        readTimer.stop();
        compareTimer.start();
        if (ierr >= 0)
            ierr = doCompare(arr, cmp, mask, hits);
        else
            hits.set(0, mask.size());
        break;
    }
    case ibis::FLOAT: {
        readTimer.start();
        ibis::array_t<float> arr;
        if (mask.cnt() * mulfactor >= mask.size()) {
            // get all values
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getValuesArray(&arr);
        }
        else { // get only the selected values
            std::vector<uint32_t> tmp;
            std::vector<uint64_t> coords;
            ibis::meshQuery::bitvectorToCoordinates(mask, getMeshShape(),
                                                    tmp);
            coords.reserve(tmp.size());
            for (unsigned i = 0; i < tmp.size(); ++ i)
                coords.push_back(tmp[i]);
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getPointValues(arr, coords);
        }
        readTimer.stop();
        compareTimer.start();
        if (ierr >= 0)
            ierr = doCompare(arr, cmp, mask, hits);
        else
            hits.set(0, mask.size());
        break;
    }
    case ibis::DOUBLE: {
        readTimer.start();
        ibis::array_t<double> arr;
        if (mask.cnt() * mulfactor >= mask.size()) {
            // get all values
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getValuesArray(&arr);
        }
        else { // get only the selected values
            std::vector<uint32_t> tmp;
            std::vector<uint64_t> coords;
            ibis::meshQuery::bitvectorToCoordinates(mask, getMeshShape(),
                                                    tmp);
            coords.reserve(tmp.size());
            for (unsigned i = 0; i < tmp.size(); ++ i)
                coords.push_back(tmp[i]);
            ierr = reinterpret_cast<const FQ_Variable*>((*it).second)
                ->getPointValues(arr, coords);
        }
        readTimer.stop();
        compareTimer.start();
        if (ierr >= 0)
            ierr = doCompare(arr, cmp, mask, hits);
        else
            hits.set(0, mask.size());
        break;
    }
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_Part::doScan:"
            << " unable to process column " << cmp.colName()
            << " with unsupport data type "
            << ibis::TYPESTRING[(int)(*it).second->type()];
        hits.set(0, mask.size());
    }
    }
    compareTimer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_Part::getValues\t"
        << readTimer.CPUTime() << "\t" << readTimer.realTime() << "\t"
        << mpi_rank << "\nStatistic\tFQ_Part::doCompare\t"
        << compareTimer.CPUTime() << "\t" << compareTimer.realTime()
        << "\t" << mpi_rank;
    return hits.cnt();
} // FQ_Part::doScan

bool FQ_Part::isValid (const std::string &func) const {
    if (! valid) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << func.c_str()
            << ": partition was not initialized correctly";
        return false;
    } else {
        return true;
    }
} // FQ_Part::isValid
