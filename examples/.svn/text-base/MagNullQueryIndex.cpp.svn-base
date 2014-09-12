/**
   Print out the number of records selected by the specified conditions.
   If more than one time step is present in the data file, the same set of
   conditions will be applied to each time step separatedly.

   command line arguments:
   -f name-of-hdf5-file
   -q query-conditions-in-a-single-string
   -p name-of-variable-to-print-or-histogram
   -h use-h5part-data
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

#define MAX_MPI_LEN 10000

static const char *options="f:F:n:N:m:M:q:Q:p:P:i:I:v:V:dDl:L:bBg:G:t:T:";
static const char *condstring = 0;
static std::string datafile;
static std::string indexfile;
static std::ostringstream logfile;
static const char *varPath = 0;
static const char *varName = "B";
static bool useBoxSelection = false;
static bool xport = false;
static int mpi_len = MAX_MPI_LEN;
static int mpi_dim = 0;
static int mpi_size=0, mpi_rank=0;
static double threshold = 0.1;

static std::string varPathStr;
static std::string varNameStr;
static int firstTimestep = 0;
static int lastTimestep = 0;

static FQ::FileFormat model = FQ::FQ_BP;
static ADIOS_File* adiosDataFile = 0;

// forward declaration
static int do_query(QueryProcessor*);


static void 
parseArgs(int argc, char **argv) {
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
	    //case 'p':
	    //case 'P': varPath = optarg; break;
	case 'c':
	    //case 'q':
	    //case 'Q': condstring = optarg; break;
	case 'v':
	case 'V': ibis::gVerbose = atoi(optarg); break;
	case 'b':
	case 'B': useBoxSelection = true; break;
	case 'x':
	case 'X': xport = true; break;
	    //case 'm':
	    //case 'M': fileModel = optarg; break;
	case 'l':
	case 'L': mpi_len = atoi(optarg);break; 
	case 't': threshold = atof(optarg);break; 
	case 'T': lastTimestep = atoi(optarg);break; 
	default: break;
        } // switch
    } // while
} // parseArgs

int 
main(int argc, char **argv) {
    int ret=0;

    parseArgs(argc, argv);
    
    if (datafile.empty()) {
	std::cerr << "Usage:\n" << *argv 
                  << " -f data-file-name" 
	    " [-i index-file-name]"
	    " [-t threshold for /var/B < threshold]"
	    " [-T last timestep]"
	    " [-g log-file-name]"
	    " [-b use-boundingbox-data-selection]"
	    " [-v verboseness]"
	    " [-l mpi-subarray-length]"
	    "\n\n e.g:   ./queryIndex -f h5uc-data-index.h5 -q 'px < 0.3' "
	    "-n y -p TimeStep2\n"
		  << std::endl;
	return -1;
    }
#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	<< "MagNullQueryIndex is running with mpi_size = " << mpi_size << " ...";
    adios_init_noxml ();
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 1000);
    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");

    if (! indexfile.empty()) {
	LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: using indexfile \"" << indexfile.c_str() << "\" ...";
    }
    if (varPath != 0) {
	LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "Debug: use variable path \"" << varPath << "\"";
        varPathStr = varPath;
    }
    if (varName != 0) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "Debug: get data for variable \"" << varName << "\"";
        varNameStr = varName;
    }
    if (logfile.str().empty() != true){
        logfile << "-" << mpi_rank << ".log";
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "Debug: using logfile \"" << logfile.str().c_str() << "\"";
    }

    if (mpi_len != MAX_MPI_LEN && mpi_size == 1) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: ignore mpi_len for a single process ...";
        mpi_len = MAX_MPI_LEN;
    }

    LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	<< "Threshold = " << threshold << ", Timestep = " << lastTimestep;


    adiosDataFile = new ADIOS_File(datafile);
    if (adiosDataFile == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to new ADIOS_File() \""
	    << datafile.c_str() << "\"";
	return -1;
    }
    if (0 == adiosDataFile->getHandle()) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to open the file \""
	    << datafile.c_str() << "\"";
	return -2;
    }

    firstTimestep = 0;

    ibis::horometer timer;
    timer.start();

    QueryProcessor* queryProcessor = new QueryProcessor
	(datafile, model, indexfile, ibis::gVerbose, "", logfile.str().c_str());
    ret = do_query(queryProcessor);
    LOGGER(ret < 0 && ibis::gVerbose >= 0)
	<< "Warning -- do_query return value is " << ret;

    timer.stop();

    LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	<< "Running seconds: Total = " << timer.realTime() << ", CPU = "
	<< timer.CPUTime();

    if (queryProcessor)
        delete(queryProcessor);

    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    return ret;
}

static int
do_query(QueryProcessor* queryProcessor) 
{
    // open the named file
    if (queryProcessor->isValid() == false) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: failed to initiate the QueryProcessor object for file \"" 
	    << datafile.c_str() << "\" ...\n"
	    << "REPORT: failed to complete processing query";
	return -1;
    }


    std::string variable;
    std::vector<uint64_t> dims;
    FQ::DataType type;    
    if (!queryProcessor->getVariableInfo(varNameStr, variable, dims,
					 &type, varPathStr)) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to get the information for variable \""
	    << variable.c_str() 
	    << "\" from file \"" << datafile.c_str() << "\"";
	return -2;
    }

    std::ostringstream oss;
    oss << "B < " << threshold;
    std::string cond = oss.str();
    condstring = cond.c_str();
    LOGGER(ibis::gVerbose > 0)
	<< "Query condition: " << condstring;

    for (unsigned time=firstTimestep; time<=lastTimestep; time++) {

	std::ostringstream oss2;
	oss2 << "/time" << time;
	varPathStr = oss2.str();

	unsigned int hits = 0;
	// getNumHits
	hits = queryProcessor->getNumHits(condstring, varPathStr, mpi_dim,
					  mpi_len);
	LOGGER(ibis::gVerbose > 0)
	    << "Debug: conditions \"" << condstring 
	    << "\" number of hits " << hits;;

	if (hits == 0) {
	    LOGGER(ibis::gVerbose > 1)
		<< "Warning -- No element is seleteced ==>"
		<< " the rest of the test is skipped!";
	    LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
		<< "REPORT: successfully complete processing query with " 
		<< hits << " hits";
	    return hits;
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
	if (xport) {
	    ibis::util::logger lg;
	    lg() << "selected coordinates:";
	    for(int i=0; i<coords.size(); i++) {
		lg() << "\n" << coords[i];
	    }
	}
	if (hits != hits1) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Error -- number of hits does not match!"
		<< "\nREPORT: failed to complete processing query";
	    return -1;
	}
    
	if (varNameStr != "") {
	    // getSelectedData
	    bool berr = true;
	    std::vector <double> values;
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
		    << "Error -- Data type is not supported!"
		    << "\nREPORT: failed to complete processing query";

		delete(queryProcessor);
		adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
		MPI_Finalize();
#endif
		return -1;
	    }
	    if (! berr) {
		LOGGER(ibis::gVerbose > 0)
		    << "Error -- Failed to get selected data"
		    << "\nREPORT: failed to complete processing query";

		delete(queryProcessor);
		adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
		MPI_Finalize();
#endif
		return -1;
	    }

	    if (xport) {
		ibis::util::logger lg;
		lg() << "Selected data:";
		for (unsigned int i=0; i<values.size(); i++) {
		    lg() << "\n" << values[i];
		}
	    }

	    // executeEqualitySelectionQuery
	    std::vector<uint64_t> eqCoords;
	    eqCoords.reserve(coords.size());
	    if (coords.size() == 0) {
		LOGGER(ibis::gVerbose > 1)
		    << "Warning -- No element is seleteced ==>"
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
			(varNameStr, values, eqCoords, "",
			 FQ::POINTS_SELECTION, mpi_dim, mpi_len);
		}
		bool match = (eqCoords.size() == coords.size());
		for (unsigned int i = 0; match && i < eqCoords.size(); i++)
		    match = (eqCoords[i] == coords[i]);
		if (! match) {
		    LOGGER(ibis::gVerbose > 0)
			<< "Warning -- the equality selection (as a list of "
			"doubles) produced a different set of coordinates as "
			"input " << eqCoords.size() << " != " << coords.size()
			<< "\nREPORT: failed to complete processing query";

		    delete(queryProcessor);
		    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
		    MPI_Finalize();
#endif
		    return -1;
		} else {
		    LOGGER(ibis::gVerbose > 1)
			<< "Debug: the equality selection (as a list of doubles)"
			<< " produced the same set of coordinates"; 
		}
	    }
	}
	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << "REPORT: successfully complete processing query with " 
	    << hits << " hits";

	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << "query: timestep " << time << " is processed...\n";
    }

    return 0;
}
