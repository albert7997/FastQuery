/**
   Print out the number of records selected by the specified conditions.
   If more than one time step is present in the data file, the same set of
   conditions will be applied to each time step separatedly.

   command line arguments:
   -d name-of-data-file
   -i name-of-index-file
   -q query-conditions-in-a-single-string
   -p path-of-the-variables
   -n name-of-variable-to-print-or-histogram
   -t time-step-to-query

   Each option -f and -q should be specified only once.  If multiple of
   them are specified, the last one will override all previous ones.
   Options -v and -t are cumulative.  The arguments to -v will be added
   together and the arguments to -t will be cumulated and sorted in a
   set.

*/

#include "BPArrayIODriver.h"
#include "queryProcessor.h"
#include <iostream>

static const char *options="d:D:f:F:n:N:q:Q:p:P:i:I:v:V:xXl:L:bBg:G:";
static char *condstring = 0;
static std::string datafile;
static std::string indexfile;
static std::ostringstream logfile;
static char *varPath = 0;
static char *varName = 0;
static bool useBoxSelection = false;
static bool xport = false;
static int mpi_len = 100000;
static int mpi_dim = 0;
static int mpi_size=1, mpi_rank=0;

static FQ::FileFormat model = FQ::FQ_BP;
static QueryProcessor* queryProcessor = 0;

static void
parseArgs(int argc, char **argv)
{
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
        case 'd':
        case 'D':
        case 'f':
        case 'F': datafile = optarg; break;
        case 'i':
        case 'I': indexfile = optarg; break;
        case 'g':
        case 'G': logfile << optarg; break;
        case 'n':
        case 'N': varName = optarg; break;
        case 'p':
        case 'P': varPath = optarg; break;
        case 'c':
        case 'q':
        case 'Q': condstring = optarg; break;
        case 'v':
        case 'V': ibis::gVerbose = atoi(optarg); break;
        case 'b':
        case 'B': useBoxSelection = true; break;
        case 'x':
        case 'X': xport = true; break;
        case 'l':
        case 'L': mpi_len = atoi(optarg);break;
        default: break;
        } // switch
    } // while
} // parseArgs

static int
do_query(const std::string& varNameStr, const std::string& varPathStr)
{
    unsigned int hits = 0;
    // getNumHits
    hits = queryProcessor->getNumHits(condstring, varPathStr, mpi_dim, mpi_len);
    LOGGER(ibis::gVerbose > 1)
        << "do_query(" << condstring << ") produced " << hits << " hit"
        << (hits>1?"s":"");

    if (hits <= 0) {
        LOGGER(ibis::gVerbose > 1)
            << "  No element is seleteced ==> the rest of the test is skipped!";
        return hits;
    }
    if (varNameStr.empty())
        return hits;

    std::string variable;
    std::vector<uint64_t> dims;
    FQ::DataType type;
    if (!queryProcessor->getVariableInfo
        (varNameStr, variable, dims, &type, varPathStr)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning: Failed to get the information for variable \""
            << variable.c_str()  << "\" from file \"" << datafile.c_str()
            << "\"";

        return -1;
    }

    // executeQuery
    uint64_t hits1=0;
    std::vector<uint64_t> coords;
    if (useBoxSelection) {
        coords.reserve(hits*2*dims.size());
        hits1 = queryProcessor->executeQuery
            ((char*)condstring, coords, varPathStr, FQ::BOXES_SELECTION,
             mpi_dim, mpi_len);
    } else if (varPath != 0) {
        coords.reserve(hits*dims.size());
        hits1 = queryProcessor->executeQuery
            ((char*)condstring, coords, varPathStr, FQ::POINTS_SELECTION,
             mpi_dim, mpi_len);
    } else {
        coords.reserve(hits*dims.size());
        hits1 = queryProcessor->executeQuery
            ((char*)condstring, coords, "", FQ::POINTS_SELECTION,
             mpi_dim, mpi_len);
    }
    LOGGER(ibis::gVerbose > 1)
        << "do_query(" << condstring << ") retrieved " << hits1 << " point"
        << (hits1>1?"s":"");
    if (hits != hits1) {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- number of hits does not match";
        return -1;
    }
    if (xport) {
        ibis::util::logger lg;
        lg() << "do_query(" << condstring << ") selected coordinates:";
        for(int i = 0; i < coords.size();) {
            const int bnd = i + dims.size();
            lg() << "\n";
            while (i < bnd) {
                lg() << ' ' << coords[i];
                ++ i;
            }
        }
    }

    // getSelectedData
    bool berr = true;
    std::vector<double> values;
    switch(type) {
    case FQ::FQT_BYTE: {
        std::vector<char> data;
        data.resize(hits1);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    case FQ::FQT_FLOAT: {
        std::vector<float> data;
        data.resize(hits);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    case FQ::FQT_DOUBLE: {
        std::vector<double> data;
        data.resize(hits);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    case FQ::FQT_INT: {
        std::vector<int32_t> data;
        data.resize(hits);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    case FQ::FQT_LONG: {
        std::vector<int64_t> data;
        data.resize(hits);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    case FQ::FQT_UINT: {
        std::vector<uint32_t> data;
        data.resize(hits);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    case FQ::FQT_ULONG: {
        std::vector<uint64_t> data;
        data.resize(hits);
        if (useBoxSelection) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr,
                 FQ::BOXES_SELECTION);
        } else if (varPath != 0) {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0], varPathStr);
        } else {
            berr = queryProcessor->getSelectedData
                (varNameStr, coords, &data[0]);
        }
        values.resize(data.size());
        for(unsigned int i=0; i<data.size(); i++) {
            values[i] = (double)data[i];
        }
        break;}
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- Data type is not supported";
        return -1;
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- Failed to get selected data";
        return -1;
    }
    LOGGER(ibis::gVerbose > 1)
        << "do_query(" << condstring << ") retrieved " << values.size()
        << " value" << (values.size()>1?"s":"") << " for " << varNameStr;
    if (xport) {
        ibis::util::logger lg;
        lg() << "do_query(" << condstring << ") selected data ["
             << values.size() << ']';
        for (unsigned int i = 0; i < values.size(); ++ i) {
            lg() << "\n" << values[i];
        }
    }

    // executeEqualitySelectionQuery
    std::vector<uint64_t> eqCoords;
    eqCoords.reserve(coords.size());
    if (coords.size() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "  No element is seleteced ==>"
            << " executeEqualitySelectionQuery is skipped!";
    } else {
        if (useBoxSelection) {
            queryProcessor->executeEqualitySelectionQuery
                (varNameStr, values, eqCoords, varPathStr,
                 FQ::BOXES_SELECTION, mpi_dim, mpi_len);
        } else if (varPath !=0) {
            queryProcessor->executeEqualitySelectionQuery
                (varNameStr, values, eqCoords, varPathStr,
                 FQ::POINTS_SELECTION, mpi_dim, mpi_len);
        } else {
            queryProcessor->executeEqualitySelectionQuery
                (varNameStr, values, eqCoords, "", FQ::POINTS_SELECTION,
                 mpi_dim, mpi_len);
        }
        bool match = (eqCoords.size() == coords.size());
        for (unsigned int i = 0; match && i < eqCoords.size(); i++)
            match = (eqCoords[i] == coords[i]);
        if (! match) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- the equality selection (as a list of doubles)"
                << " produced a different set of coordinates as input "
                << eqCoords.size() << " != " << coords.size();
            return -1;
        } else {
            LOGGER(ibis::gVerbose > 2)
                << " the equality selection (as a list of doubles)"
                << " produced the same set of coordinates";
        }
    }

    LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
        << "do_query successfully complete the query with "
        << hits << " hits" << std::endl;

    return hits;
} // do_query

int
main(int argc, char **argv)
{
    parseArgs(argc, argv);

    if (datafile.empty() || condstring == 0) {
        std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
            " -q query-conditions-in-a-single-string"
            " [-i index-file-name]"
            " [-g log-file-name]"
            " [-n name-of-variable]"
            " [-p path-of-variable]"
            " [-b use-boundingbox-data-selection]"
            " [-v verboseness]"
            " [-l mpi-subarray-length]"
            "\n e.g:   ./queryIndex -f h5uc-data-index.h5 -q 'px < 0.3' -n y -p TimeStep2\n\n"
            "\tFor More detailed usage description and examples, please see file GUIDE"
                  << std::endl;
        return -1;
    }

#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
    adios_init_noxml();
    adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, 50);

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
        << *argv << " starts with mpi_size= " << mpi_size << " ...";

    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");

    std::string varNameStr;
    std::string varPathStr;

    if (! indexfile.empty()) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
            << *argv << " uses indexfile \"" << indexfile.c_str() << "\"";
    } else {
        std::ostringstream oss;
        oss << datafile;
        long pos = oss.tellp();
        oss.seekp(pos-3);
        oss << "_idx.bp";
        indexfile = oss.str();
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
            << *argv << " uses the default indexfile \"" << indexfile.c_str()
            << "\" ...";
    }

    if (varPath != 0) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
            << *argv << " use variable path \"" << varPath << "\"";
        varPathStr = varPath;
    }
    if (varName != 0) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
            << " get data for variable \"" << varName << "\"";
        varNameStr = varName;
    }
    if (logfile.str().empty() != true){
        logfile << "-" << mpi_rank << ".log";
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
            << " using logfile \"" << logfile.str().c_str() << "\"";
    }

    // the file handler
    queryProcessor = new QueryProcessor
        (datafile, model, indexfile, ibis::gVerbose, "", logfile.str().c_str());
    // open the named file
    if (queryProcessor->isValid() == false) {
        LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
            << "Warning: failed to initiate the QueryProcessor on file \""
            << datafile.c_str() << "\"";

        return -1;
    }

    ibis::horometer timer;
    timer.start();

    do_query(varNameStr, varPathStr);

    timer.stop();

    int idx = 0;
    const int IDX_TOT_TIME = idx++;
    const int IDX_CPU_TIME = idx++;
    const int IDX_READ_DATA_TIME = idx++;
    const int IDX_WRITE_INDEX_TIME = idx++;

    double timeArray[idx];
    timeArray[IDX_TOT_TIME] = timer.realTime();
    timeArray[IDX_CPU_TIME] = timer.CPUTime();
    timeArray[IDX_READ_DATA_TIME] =
        BPArrayIODriver::getReadDataTimer().realTime();
    timeArray[IDX_WRITE_INDEX_TIME] =
        BPArrayIODriver::getWriteIndexTimer().realTime();

#ifndef FQ_NOMPI
    double timeArrayArray[idx*mpi_size];
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(timeArray, idx, MPI_DOUBLE, timeArrayArray, idx, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
#else
    double *timeArrayArray = timeArray;
#endif
    if (mpi_rank == 0) {
        double totTimeArray[mpi_size];
        double cpuTimeArray[mpi_size];
        double readDataTimeArray[mpi_size];
        double writeIndexTimeArray[mpi_size];
        double sumTimeArray[idx];
        for (int i=0; i<mpi_size; i++) {
            totTimeArray[i] = 0;
            cpuTimeArray[i] = 0;
            readDataTimeArray[i] = 0;
            writeIndexTimeArray[i] = 0;
        }
        for (int i=0; i<idx; i++) {
            sumTimeArray[i] = 0;
        }
        for (int i=0; i<mpi_size; i++) {
            sumTimeArray[IDX_TOT_TIME] +=
                timeArrayArray[i*idx+IDX_TOT_TIME];
            sumTimeArray[IDX_CPU_TIME] +=
                timeArrayArray[i*idx+IDX_CPU_TIME];
            sumTimeArray[IDX_READ_DATA_TIME] +=
                timeArrayArray[i*idx+IDX_READ_DATA_TIME];
            sumTimeArray[IDX_WRITE_INDEX_TIME] +=
                timeArrayArray[i*idx+IDX_WRITE_INDEX_TIME];

            totTimeArray[i] = timeArrayArray[i*idx+IDX_TOT_TIME];
            cpuTimeArray[i] = timeArrayArray[i*idx+IDX_CPU_TIME];
            readDataTimeArray[i] = timeArrayArray[i*idx+IDX_READ_DATA_TIME];
            writeIndexTimeArray[i] = timeArrayArray[i*idx+IDX_WRITE_INDEX_TIME];

            printf("Running seconds for %d: Total=%.3f, CPU=%.3f, "
                   "ReadData=%.3f, WriteIndex=%.3f\n",
                   i, timeArrayArray[i*idx+IDX_TOT_TIME],
                   timeArrayArray[i*idx+IDX_CPU_TIME],
                   timeArrayArray[i*idx+IDX_READ_DATA_TIME],
                   timeArrayArray[i*idx+IDX_WRITE_INDEX_TIME]);
        }

        printf("MPI size=%d, Mean Total=%.3f, CPU=%.3f, ReadData=%.3f, "
               "WriteIndex=%.3f, Median Total=%.3f, CPU=%.3f, ReadData=%.3f, "
	       "WriteIndex=%.3f\n",
               mpi_size, sumTimeArray[IDX_TOT_TIME]/mpi_size,
               sumTimeArray[IDX_CPU_TIME]/mpi_size,
               sumTimeArray[IDX_READ_DATA_TIME]/mpi_size,
               sumTimeArray[IDX_WRITE_INDEX_TIME]/mpi_size,
               util::compute_median(totTimeArray, mpi_size),
               util::compute_median(cpuTimeArray, mpi_size),
               util::compute_median(readDataTimeArray, mpi_size),
               util::compute_median(writeIndexTimeArray, mpi_size));
    }

    delete queryProcessor;

    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    return 0;
} // main
