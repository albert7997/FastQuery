/**
   A simple program to test the capability of building and storing
   indexes to a file.

   binning option is specified in the form of "<binning ... />".  
   NOTE that on most systems, the binning option needs to be quoted 
   because it involved characters that have special meaning to most shells.
*/

#include "BPArrayIODriver.h"
#include "indexBuilder.h"
#include <iostream>	

#define MAX_MPI_LEN 10000
static const char *options="b:B:f:F:p:P:n:N:i:I:v:V:m:M:l:L:g:G:rR:t:T:";
static int verboseness = 1;

static std::string datafile;
static std::string indexfile;
static std::ostringstream logfile;
static char *binning = 0;
static bool forcerebuild = false;
static const char* varPath = 0;
static const char* varName = "B";
static int mpi_len = MAX_MPI_LEN;
static double threshold = 0.1;
static int mpi_size=0, mpi_rank=0;

static std::string varPathStr;
static std::string varNameStr;
static int firstTimestep = 0;
static int lastTimestep = 0;

static ADIOS_File* adiosDataFile = 0;


// forward declaration
static int do_query();


static void 
parseArgs(int argc, char **argv) 
{
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'f':
	case 'F': datafile = optarg; break;
	    //case 'i':
	    //case 'I': indexfile = optarg; break;
	case 'g':
	case 'G': logfile << optarg; break;
	    //case 'p':
	    //case 'P': varPath = optarg; break;
	    //case 'n':
	    //case 'N': varName = optarg; break;
	case 'v':
	case 'V': verboseness = atoi(optarg); break;
	case 'r':
	case 'R': forcerebuild = true; break;
	case 'l':
	case 'L': mpi_len = atoi(optarg);break; 
	case 'b':
	    //case 'B': binning = optarg; break;
	case 'B': binning = new char[128]; sprintf(binning, "<binning %s />", optarg); break;
	case 't': threshold = atof(optarg);break;
	case 'T': lastTimestep = atoi(optarg);break;
	default: break;
        } // switch
    } // while
} // parseArgs


int 
main(int argc, char** argv) 
{
    int ret = 0;
    parseArgs(argc, argv);
    if (datafile.empty()) {
	std::cerr << "Usage:\n" << *argv 
                  << " -f data-file-name [-i index-file-name] [-g log-file]" 
		  << " [-r (force-rebuild-index)] [-v verboseness]"
                  << " [-t threshold for /var/B < threshold]"
                  << " [-T last timestep]"
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
    adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, 1000);

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	<< "MagNullQueryScan is running with mpi_size = " << mpi_size << " ...";

    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    ibis::gParameters().add("fileManager.maxBytes", "6GB");
    if (forcerebuild) {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "true");
    } else {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "false");
    }

    if (! indexfile.empty()) {
	LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: using indexfile \"" << indexfile.c_str() << "\" ...";
    } else {
        std::ostringstream oss;
        oss << datafile;
        long pos = oss.tellp();
        oss.seekp(pos-3);
        oss << "_idx" << mpi_size << ".bp";
        indexfile = oss.str();
	LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: using auto indexfile \"" << indexfile.c_str() << "\" ...";
    }

    if (binning !=0) {
      	LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: using binning option \"" << binning << "\" ...";
    }

    if (varPath !=0) {
        LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: specifying variable path \"" << varPath << "\" ...";
	varPathStr = varPath;
    }

    if (varName !=0) {
        LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: specifying variable name \"" << varName << "\" ...";
	varNameStr = varName;
    }

    if (logfile.str().empty() != true){
        LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: using logfile \"" << logfile.str().c_str() << "\" ...";
    	logfile << "-" << mpi_rank << ".log";
    }

    if (mpi_len != MAX_MPI_LEN && mpi_size == 1) {
        LOGGER(verboseness > 1 && mpi_rank == 0)
	    << "DEBUG: ignore mpi_len for a single process ...";
        mpi_len = MAX_MPI_LEN;
    }

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	<< "Threshold = " << threshold << ", Timestep = " << lastTimestep;


    ibis::horometer timer;
    timer.start();

    ret = do_query();
    LOGGER(ret < 0 && ibis::gVerbose >= 0)
	<< "DEBUG: do_query return value is " << ret;

    timer.stop();

    LOGGER(ibis::gVerbose > 0)
	<< "MPI size=" << mpi_size << ", Running seconds: Total="
	<< timer.realTime() << ", CPU = " << timer.CPUTime();

    if (adiosDataFile)
	delete adiosDataFile;

    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    return ret;
} // main


static int
do_query()
{
    adiosDataFile = new ADIOS_File(datafile);
    if (adiosDataFile == 0) {
	LOGGER(verboseness > 0)
	    << "ERROR: Failed to new ADIOS_File() \""
	    << datafile.c_str() << "\"";
	return -1;
    }
    if (0 == adiosDataFile->getHandle()) {
	LOGGER(verboseness > 0)
	    << "ERROR: Failed to open the file \""
	    << datafile.c_str() << "\"";
	return -2;
    }

    firstTimestep = 0;

    ADIOS_Var* v5 = adiosDataFile->getVariable("/var/v5");
    ADIOS_Var* v6 = adiosDataFile->getVariable("/var/v6");
    ADIOS_Var* v7 = adiosDataFile->getVariable("/var/v7");
    if (v5==0 || v6==0 || v7==0) {
	LOGGER(verboseness > 0)
	    << "ERROR: No such variable v5,v6,v7 in record";
	return -3;
    }

    std::vector<uint64_t> dims(v5->getNumDimension());
    v5->getDimensionArray(&dims);

    int nelem = 1;
    std::vector<uint64_t> start(dims.size());
    std::vector<uint64_t> count(dims.size());
    start[0] = 0;
    count[0] = 1;
    for (int i=1; i<dims.size(); i++) {
	nelem *= dims[i];
	start[i] = 0;
	count[i] = dims[i];
    }
    start[1] = mpi_len*mpi_rank;
    if (start[1]+mpi_len <= dims[1])
	count[1] = mpi_len;
    else
	count[1] = dims[1]-start[1];

    int nlocalelem = 1;
    for (int i=0; i<count.size(); i++) {
	nlocalelem *= count[i];
    }
    LOGGER(ibis::gVerbose > 0)
	<< "Number of elements is " << nelem << ", local elements is "
	<< nlocalelem;

    double* v5_data = new double[nlocalelem];
    double* v6_data = new double[nlocalelem];
    double* v7_data = new double[nlocalelem];

    for (unsigned time=firstTimestep; time<=lastTimestep; time++) {
	start[0] = time-firstTimestep;

	int64_t nbytes = v5->readData(v5_data, -1, start, count);
	if (nbytes <= 0) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- Can't get data for v5";
	    return -4;
	}

	nbytes = v6->readData(v6_data, -1, start, count);
	if (nbytes <= 0) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- Can't get data for v6";
	    return -5;
	}

	nbytes = v7->readData(v7_data, -1, start, count);
	if (nbytes <= 0) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- Can't get data for v7";
	    return -6;
	}

	int cnt=0;
	int nhits = 0;
	for (int k=0; k<nlocalelem; k++) {
	    double val = sqrt(v5_data[k]*v5_data[k] + v6_data[k]*v6_data[k] +
			      v7_data[k]*v7_data[k]);
	    if (val < threshold) ++cnt;
	}
	LOGGER(ibis::gVerbose > 0)
	    << "[" << mpi_rank << "] Local counter for B < " << threshold
	    << " = " << cnt;
#ifndef FQ_NOMPI
	int nhitsArray[mpi_size];
	MPI_Gather(&cnt, 1, MPI_INT, nhitsArray, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	for (int i=0; i<mpi_size; i++)
	    nhits += nhitsArray[i];
#else
	nhits = cnt;
#endif
	LOGGER(ibis::gVerbose >=0)
	    << "Number of hits for timestep " << time << " = " << nhits;
    }

    delete[] v5_data;
    delete[] v6_data;
    delete[] v7_data;

    LOGGER(ibis::gVerbose > 5)
	<< "do_query(): operation completed";

    return 0;
}

