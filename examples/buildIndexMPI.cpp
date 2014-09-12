/**
   A simple program to test the capability of building and storing
   indexes to a file.

   binning option is specified in the form of "<binning ... />".  
   NOTE that on most systems, the binning option needs to be quoted 
   because it involved characters that have special meaning to most shells.
*/

#include "indexBuilder.h"
#include <iostream>	

static const char *options="b:B:f:F:p:P:n:N:i:I:v:V:m:M:l:L:k:K:g:G:s:S:rR";
static int verboseness = 1;

static char *datafile = 0;
static char *indexfile = 0;
static std::ostringstream logfile;
static char *fileModel = 0;
static char *binning = 0;
static bool forcerebuild = false;
static char* varPath = 0;
static char* varName = 0;
static int mpi_len = 100000;
static int mpi_dim = 0;
static int startIdx = 0;
static int endIdx = -1;
static int groupSize = -1;

void parseArgs(int argc, char **argv) {
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'm':
	case 'M': fileModel = optarg; break;
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
	case 'V': verboseness = atoi(optarg); break;
	case 'r':
	case 'R': forcerebuild = true; break;
	case 'l':
	case 'L': mpi_len = atoi(optarg);break; 
	case 'k':
	case 'K': {
		std::string str = optarg;
		size_t pos = str.find(':');
		if (pos == std::string::npos) {
		    endIdx = atoi(optarg); 
		} else {
		    startIdx = atoi(str.substr(0, pos).c_str());
		    endIdx = atoi(str.substr(pos+1).c_str());
		}
		break;
	}
	case 's':
	case 'S': groupSize = atoi(optarg); break;
	case 'b':
	case 'B': binning = optarg; break;
	default: break;
        } // switch
    } // while
} // parseArgs


int main(int argc, char** argv) {
    std::string varPathStr;
    std::string varNameStr;
    int ret = 0;
    parseArgs(argc, argv);
    if (datafile == 0) {
        std::cerr << "Usage:\n" << *argv 
                  << " -f data-file-name [-i index-file-name] [-g log-file]\n" 
		  << " [-p varPath] [-n varName]\n"
		  << " [-b '<binning nbins=1000 />' (default unbinned)]"
		  << " [-r force-rebuild-index] [-v verboseness]\n"
		  << " [-m fileModel [HDF5(default), H5PART, NETCDF, PNETCDF]]\n"
		  << " [-k [startFildIdx:]endFileIdx]\n"
		  << " [-l mpi_subarray_size(default=100000)]\n"
		  << " [-s mpi_group_size(defaul=evenly splitted)]\n"
		  << std::endl;
        return -1;
    }

#ifndef FQ_NOMPI 
    MPI_Init(&argc, &argv);
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int group_rank = mpi_rank;
#endif

    ibis::gParameters().add(FQ_REPORT_STATISTIC, "true");
    ibis::gParameters().add("fileManager.maxBytes", "6GB");
    if (forcerebuild) {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "true");
    } else {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "false");
    }

    FQ::FileFormat model = FQ::FQ_HDF5;
    if (fileModel != 0) {
	std::string format = fileModel;
    	if (format.compare("HDF5") == 0) {
	    model = FQ::FQ_HDF5;
    	} else if (format.compare("H5PART") == 0) {
	    model = FQ::FQ_H5Part;
    	} else if (format.compare("NETCDF") == 0) {
	    model = FQ::FQ_NetCDF;
    	} else if (format.compare("PNETCDF") == 0) {
	    model = FQ::FQ_pnetCDF;
	}
    }


    if (indexfile != 0) {
	if (verboseness > 1) std::cout << "DEBUG: using indexfile \"" << indexfile << "\" ... \n";
    }
    if (binning !=0) {
      	if (verboseness > 1) std::cout << "DEBUG: using binning option \"" << binning << "\" ...\n";
    }
    if (varPath !=0) {
        if (verboseness > 1) std::cout << "DEBUG: specifying variable path \"" << varPath << "\" ...\n";
	varPathStr = varPath;
    }
    if (varName !=0) {
        if (verboseness > 1) std::cout << "DEBUG: specifying variable name \"" << varName << "\" ...\n";
	varNameStr = varName;
    }

    if (logfile.str().empty() != true){
        if (verboseness > 1) std::cout << "DEBUG: using logfile \"" << logfile.str().c_str() << "\" ...\n";
#ifndef FQ_NOMPI
    	logfile << "/" << mpi_rank << ".log";
#endif
    }



    bool berr = true;
    ibis::horometer totTimer;
    //std::cout << "start:end " << startIdx << ":" << endIdx << std::endl;
    if (endIdx == -1) {
        totTimer.start();
	IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile, 
		verboseness, NULL, logfile.str().c_str()); // the file handler
	if (indexBuilder->isValid() == false) {
	    berr = false;
	    if (verboseness > 0) {
		std::cout << "ERROR: Failed to initiate the IndexBuilder object for file \"" 
			<< datafile << "\"" << std::endl;
		std::cout << "REPORT: Failed to complete building index";
		delete(indexBuilder);
#ifndef FQ_NOMPI
    		MPI_Finalize();
#endif
		return -1;
	    }
	    
	} else {
	    if (verboseness > 1) {
		std::cout << "DEBUG: initiate the IndexBuilder object for file \"" 
		    << datafile << "\"" << std::endl;
	    }
	}
 	std::vector<std::string> variables;
	int numVar = indexBuilder->getAllVariables(variables, varPathStr, varNameStr);
	int builtVar = 0;
	if (mpi_len > 0) {
	    builtVar = indexBuilder->buildIndexes(binning, varPathStr, varNameStr, mpi_dim, mpi_len);
	} else {
    	    if (varName!=0) {
	    	builtVar = indexBuilder->buildIndexes(binning, varPathStr, varNameStr);
    	    } else if (varPath!=0) {
	    	builtVar = indexBuilder->buildIndexes(binning, varPathStr);
	    } else {
		builtVar = indexBuilder->buildIndexes(binning);
	    }
    	}
	if (builtVar != numVar) {
	    if (verboseness > 0) {
		std::cout << "Warning: Failed to build the index of " 
		    	<< numVar-builtVar << " variables" << std::endl;
	    }
	    berr = false;
	}
	ret = builtVar;
	delete(indexBuilder);
    } else {
#ifdef FQ_NOMPI
    	totTimer.start();
    	for (int i=startIdx; i<endIdx; i++) {
	    ibis::horometer openTimer;
	    ibis::horometer closeTimer;
	    std::ostringstream dataFileName("");; 
	    std::ostringstream indexFileName("");; 
	    dataFileName << datafile << i;
	    if (indexfile != 0) {
	    	indexFileName << indexfile << i;
	    }
	    //open file
	    openTimer.start();
	    IndexBuilder* indexBuilder = new IndexBuilder(dataFileName.str(), model, indexFileName.str(), 
			verboseness, NULL, logfile.str().c_str()); // the file handler
	    if (indexBuilder->isValid() == false) {
		if (verboseness > 0) {
		    std::cout << "ERROR: Failed to initiate the IndexBuilder object for file \"" 
			<< dataFileName.str().c_str() << "\"" << std::endl;
		}
		berr = false;
		continue;
	    } else {
		if (verboseness > 1) {
		    std::cout << "DEBUG: initiate the IndexBuilder object for file \"" 
			<< dataFileName.str().c_str() << "\"" << std::endl;
		}
	    }
	    openTimer.stop();
	    // build indices with bins
 	    std::vector<std::string> variables;
	    int numVar = indexBuilder->getAllVariables(variables, varPathStr, varNameStr);
	    int builtVar = 0;
	    if (varName!=0) {
		builtVar = indexBuilder->buildIndexes(binning, varPathStr, varNameStr);
	    } else if (varPath!=0) {
		builtVar = indexBuilder->buildIndexes(binning, varPathStr);
	    } else {
		builtVar = indexBuilder->buildIndexes(binning);
	    }
	    if (builtVar != numVar) {
	    	if (verboseness > 0) {
		    std::cout << "Warning: Failed to build the index of " 
			<< numVar-builtVar << " variables for file \""
			<< dataFileName.str().c_str() << "\"" << std::endl;
		}
		berr = false;
	    } 
	    ret += builtVar;
	    // close file
	    closeTimer.start(); 
	    delete(indexBuilder);
	    closeTimer.stop();

    	    if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC)) {
		LOGGER(true) << "Statistic\tFile_open_time\t" 
		    << openTimer.CPUTime() << "\t" << openTimer.realTime();
		LOGGER(true) << "Statistic\tFile_close_time\t" 
		    << closeTimer.CPUTime() << "\t" << closeTimer.realTime();
	    }
	    if (verboseness > 0) {
	        std::cout << "REPORT: built " << builtVar << " variables for file \"" 
			<< dataFileName.str().c_str() << "\"" << std::endl;
	    }
	}
#else
	int num_file = endIdx-startIdx;
    	int num_group;
    	int groupId = -1;

	std::vector<int*> group_ranks;
    	MPI_Group org_group, new_group;
	if (groupSize <0) {
	    //evenly split group size based on the number of cores
	    num_group = num_file;
	    if (num_group > mpi_size) {
	    	num_group = mpi_size;
	    }
	} else {
	    //use fix group size given from user
	    num_group = ceil(mpi_size/groupSize);
	}
	
        // determine the group size and ranks
	int group_size[num_group];
	int idx = 0;
	for (int i=0; i<num_group; i++) {
	    if (groupSize < 0) {
	    	group_size[i] = floor((float)mpi_size/(float)num_group);
		if (i < mpi_size%num_group) group_size[i]++;
	    } else {
		if (mpi_size%groupSize != 0 && i == (num_group-1)) {
		    group_size[i] = mpi_size%groupSize;
		} else {
		    group_size[i] = groupSize;
		}
	    }
	    int* ranks = new int[group_size[i]];
	    for (int j=0; j<group_size[i]; j++) {
		ranks[j] = idx;
		if (idx == mpi_rank) {
		    groupId = i;
		}
		idx++;
	    }
	    group_ranks.push_back(ranks);
	}
	//construct groups
        MPI_Comm_group(MPI_COMM_WORLD, &org_group); 
	MPI_Group_incl(org_group, group_size[groupId], group_ranks[groupId], &new_group);
	MPI_Group_rank(new_group, &group_rank);
    	MPI_Comm mpi_comm;
	MPI_Comm_create(MPI_COMM_WORLD, new_group, &mpi_comm);
        MPI_Barrier(MPI_COMM_WORLD);
        totTimer.start();
	for (int fileIdx=groupId+startIdx; fileIdx<endIdx; fileIdx+=num_group) {
	    ibis::horometer openTimer;
	    ibis::horometer closeTimer;
	    std::ostringstream dataFileName("");; 
	    std::ostringstream indexFileName("");; 
	    dataFileName << datafile << fileIdx;
	    if (indexfile != 0) {
	    	indexFileName << indexfile << fileIdx;
	    }
	    // open file
	    openTimer.start();
	    IndexBuilder* indexBuilder = new IndexBuilder(dataFileName.str(), model, indexFileName.str(), 
			verboseness, NULL, logfile.str().c_str(), mpi_comm ); // the file handler
	    if (indexBuilder->isValid() == false) {
		if (verboseness > 0) {
		    std::cout << "ERROR: Failed to initiate the IndexBuilder object for file \"" 
			    << dataFileName.str().c_str() << "\"" << std::endl;
		}
		berr = false;
		delete(indexBuilder);
		continue;
	    } else {
		if (verboseness > 1) {
		    std::cout << "DEBUG: initiate the IndexBuilder object for file \"" 
			<< dataFileName.str().c_str() << "\"" << std::endl;
		}
	    }
	    openTimer.stop(); 
	    // build indices with bins
 	    std::vector<std::string> variables;
	    int numVar = indexBuilder->getAllVariables(variables, varPathStr, varNameStr);
	    int builtVar = 0;
	    if (mpi_len > 0) {
		builtVar = indexBuilder->buildIndexes(binning, varPathStr, varNameStr, mpi_dim, mpi_len);
	    } else {
	        if (varName!=0) {
		    builtVar = indexBuilder->buildIndexes(binning, varPathStr, varNameStr);
		} else if (varPath!=0) {
		    builtVar = indexBuilder->buildIndexes(binning, varPathStr);
		} else {
		    builtVar = indexBuilder->buildIndexes(binning);
		}
	    }
	    if (builtVar != numVar) {
		if (verboseness > 0) {
		    std::cout << "Warning: Failed to build the index of " 
			<< numVar-builtVar << " variables for file\"" 
			<< dataFileName.str().c_str() << "\"" << std::endl;
		}
		berr = false;
	    }
	    ret += builtVar;
	    closeTimer.start(); 
	    delete(indexBuilder);
	    closeTimer.stop();
    	    if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC)) {
		LOGGER(true) << "Statistic\tFile_open_time\t" 
		    << openTimer.CPUTime() << "\t" << openTimer.realTime();
		LOGGER(true) << "Statistic\tFile_close_time\t" 
		    << closeTimer.CPUTime() << "\t" << closeTimer.realTime();
	    }
	    if (verboseness > 0) {
	        std::cout << "REPORT: built " << builtVar << " variables for file \"" 
			<< dataFileName.str().c_str() << "\"" << std::endl;
	    }
	}
#endif
    }
#ifndef FQ_NOMPI
    ibis::horometer idleTimer;
    idleTimer.start();
    MPI_Barrier(MPI_COMM_WORLD);
    idleTimer.stop();
    MPI_Finalize();
    if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC)) {
    	LOGGER(true) << "Statistic\tFinal_wait_time\t" 
	    << idleTimer.CPUTime() << "\t" << idleTimer.realTime();
    }
#endif
    totTimer.stop();
    if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC)) {
    	LOGGER(true) << "Statistic\tTotal_time\t" 
	    << totTimer.CPUTime() << "\t" << totTimer.realTime();
    }
    if (verboseness > 0) {
    	if (berr == false) {
	    std::cout << "REPORT: Failed to complete building index\n";
    	} else {
	    std::cout << "REPORT: Successfully complete building index for " 
		<< ret << " variables" <<std::endl;;
	}
    }
    return ret;
} // main
