/**
   A simple program to test the capability of building and storing
   indexes to a file.

   binning option is specified in the form of "<binning ... />".  
   NOTE that on most systems, the binning option needs to be quoted 
   because it involved characters that have special meaning to most shells.
*/

#include "indexBuilder.h"
#include <iostream>	
#include "BPArrayIODriver.h"

static const char *options="b:B:f:F:d:D:p:P:n:N:i:I:v:V:m:M:l:L:g:G:rR";

static std::string datafile;
static std::string indexfile;
static std::ostringstream logfile;
static char *binning = 0; //"<binning precision=1 />";
static bool forcerebuild = false;
static char* varPath = 0;
static char* varName = 0;
static int mpi_len = 100000;
static int mpi_dim = 0;
static int mpi_size = 1, mpi_rank = 0;

static FQ::FileFormat model = FQ::FQ_BP;
static IndexBuilder* indexBuilder=0;

static void 
parseArgs(int argc, char **argv) 
{
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	    //case 'm':
	    //case 'M': fileModel = optarg; break;
	case 'd':
	case 'D':
	case 'f':
	case 'F': datafile = optarg; break;
	case 'i':
	case 'I': indexfile = optarg; break;
	case 'g':
	case 'G': logfile << optarg; break;
	case 'p':
	case 'P': varPath = optarg; break;
	case 'n':
	case 'N': varName = optarg; break;
	case 'v':
	case 'V': ibis::gVerbose = atoi(optarg); break;
	case 'r':
	case 'R': forcerebuild = true; break;
	case 'l':
	case 'L': mpi_len = atoi(optarg);break; 
	case 'b':
	    //case 'B': binning = optarg; break;
	case 'B':
	    binning = new char[128];
	    sprintf(binning, "<binning %s />", optarg);
	    break;
	default: break;
        } // switch
    } // while
} // parseArgs


static int
do_build(const std::string& varNameStr, const std::string& varPathStr)
{
    std::string pathstr = varPathStr;
    std::vector<std::string> variables;
    int numVar = indexBuilder->getAllVariables
	(variables, pathstr, varNameStr);
    if (numVar < 1) {
	// another attempt with an empty path
	pathstr.clear();
	numVar = indexBuilder->getAllVariables
	    (variables, pathstr, varNameStr);
	if (numVar < 1) {
	    LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
		<< "Warning: no variable matching ("
		<< varPathStr << "," << varNameStr << ")";
	    return -2;
	}
    }

    int builtVar = 0;
    if (mpi_len > 0) {
	builtVar = indexBuilder->buildIndexes
	    (binning, pathstr, varNameStr, mpi_dim, mpi_len);
    }
    else {
	if (varName != 0) {
	    builtVar = indexBuilder->buildIndexes
		(binning, pathstr, varNameStr);
	}
	else if (varPath != 0) {
	    builtVar = indexBuilder->buildIndexes(binning, pathstr);
	}
	else {
	    builtVar = indexBuilder->buildIndexes(binning);
	}
    }
    if (builtVar != numVar) {
	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << "Warning: Failed to build the index of " << numVar-builtVar 
	    << " variable" << (numVar-builtVar>1?"s":"");
    }
    else {
	LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	    << "Successfully complete building index for " 
	    << numVar << " variable" << (numVar>1?"s":"");
    }
    return (pathstr.empty() ? -3 : builtVar);
}


int 
main(int argc, char** argv) 
{
    parseArgs(argc, argv);
    if (datafile.empty()) {
	std::cout
	    << "Usage:\n" << *argv 
	    << " -f data-file-name [-i index-file-name] [-g log-file]" 
	    << " [-n variable-name] [-p variable-path]"
	    << " [-b '<binning nbins=1000 />' (default <binning none>)]"
	    << " [-r (force-rebuild-index)] [-v verboseness]"
	    << " [-l mpi_subarray_size(default=100000)]\n"
	    << "\tIt builds index for a set of variables whose dataset location has the prefix\n"
	    << "\tvariable-path and postfix variable-name.\n\n"
	    << "\tUse option \"-i\" to specify the output file for storing indexes.\n"
	    << "\tOtherwise, the indexes are written back to data file \"data-file-name\".\n\n"
	    << "\tUse option \"-r\" to enforce rebuild and replace the existing index.\n\n"
	    << "\tUnder parallel mode, use \"-l\" to set the subarray size for spitting dataset.\n\n"
	    << "\tUse option \"-b\" to specify the binning option to build the index.\n"
	    << "\tThe available binning option is defined and provided by the FastBit.\n"
	    << "\tMore information can be found at http://crd.lbl.gov/~kewu/fastbit/doc/indexSpec.html.\n"
	    << "\tBinning option is suggested to be used with large dataset to reduce the size of built index.\n"
	    << "\tPrecision option is suggested to be used when the query involves floating point numbers.\n\n"
	    << "\tFor More detailed usage description and examples, please see file GUIDE"
	    << std::endl;
        return -1;
    }

#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

    adios_init_noxml();
    adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, 20);

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	<< "BPbuildIndex is running with mpi_size=" << mpi_size << "...";

    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    ibis::gParameters().add("fileManager.maxBytes", "4GB");
    if (forcerebuild) {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "true");
    } else {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "false");
    }


    std::string varPathStr;
    std::string varNameStr;

    if (! indexfile.empty()) {
	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << *argv << " using indexfile \"" << indexfile.c_str() << "\" ... ";
    } else {
        std::ostringstream oss;
        oss << datafile;
        long pos = oss.tellp();
        oss.seekp(pos-3);
        oss << "_idx.bp";
        indexfile = oss.str();
        LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << *argv << " using default indexfile \"" << indexfile.c_str()
	    << "\" ... ";
    }

    if (binning != 0) {
      	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << *argv << " using binning option \"" << binning << "\" ...";
    }
    if (varPath != 0) {
        LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << *argv << " specifying variable path \"" << varPath << "\" ...";
	varPathStr = varPath;
    }
    if (varName != 0) {
        LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << *argv << " specifying variable name \"" << varName << "\" ...";
	varNameStr = varName;
    }

    if (logfile.str().empty() != true){
        LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << *argv << " using logfile \"" << logfile.str().c_str()
	    << "\" ...";
    	logfile << "-" << mpi_rank << ".log";
    }

    indexBuilder = new IndexBuilder
	(datafile, model, indexfile, ibis::gVerbose, NULL,
	 logfile.str().c_str());
    if (indexBuilder->isValid() == false) {
	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << "Warning: Failed to initialize the IndexBuilder object for \"" 
	    << datafile.c_str() << "\"";
        return -1;
    } else {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << *argv << " initialize the IndexBuilder object for file \"" 
	    << datafile.c_str() << "\"";
    }

    ibis::horometer timer;
    timer.start();

    if (varPath != 0) {
        do_build(varNameStr, varPathStr);
    }
    else {
	do_build(varNameStr, varPathStr);
    }

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
	    sumTimeArray[IDX_TOT_TIME] += timeArrayArray[i*idx+IDX_TOT_TIME];
	    sumTimeArray[IDX_CPU_TIME] += timeArrayArray[i*idx+IDX_CPU_TIME];
	    sumTimeArray[IDX_READ_DATA_TIME] +=
		timeArrayArray[i*idx+IDX_READ_DATA_TIME];
	    sumTimeArray[IDX_WRITE_INDEX_TIME] +=
		timeArrayArray[i*idx+IDX_WRITE_INDEX_TIME];

	    totTimeArray[i] = timeArrayArray[i*idx+IDX_TOT_TIME];
	    cpuTimeArray[i] = timeArrayArray[i*idx+IDX_CPU_TIME];
	    readDataTimeArray[i] = timeArrayArray[i*idx+IDX_READ_DATA_TIME];
	    writeIndexTimeArray[i] = timeArrayArray[i*idx+IDX_WRITE_INDEX_TIME];
        
	    printf("Running seconds for %d: Total=%.3f, CPU=%.3f, "
		   "ReadData=%.3f, WriteIndex=%.3f\n", i,
		   timeArrayArray[i*idx+IDX_TOT_TIME],
		   timeArrayArray[i*idx+IDX_CPU_TIME],
		   timeArrayArray[i*idx+IDX_READ_DATA_TIME],
		   timeArrayArray[i*idx+IDX_WRITE_INDEX_TIME]);
	}

	printf("MPI size=%d, Mean Total=%.3f, CPU=%.3f, "
	       "ReadData=%.3f, WriteIndex=%.3f, Median Total=%.3f, "
	       "CPU=%.3f, ReadData=%.3f, WriteIndex=%.3f\n", mpi_size,
	       sumTimeArray[IDX_TOT_TIME]/mpi_size,
	       sumTimeArray[IDX_CPU_TIME]/mpi_size,
	       sumTimeArray[IDX_READ_DATA_TIME]/mpi_size,
	       sumTimeArray[IDX_WRITE_INDEX_TIME]/mpi_size,
	       util::compute_median(totTimeArray, mpi_size),
	       util::compute_median(cpuTimeArray, mpi_size),
	       util::compute_median(readDataTimeArray, mpi_size),
	       util::compute_median(writeIndexTimeArray, mpi_size));
    }

    delete indexBuilder;

    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif

    return 0;
} // main
