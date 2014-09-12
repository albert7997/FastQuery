/**
   Print out the number of records selected by the specified conditions.

   command line arguments:
   -m file-format
   -f name-of-data-file
   -i name-of-index-file
   -g name-of-log-file
   -q query-conditions-in-a-single-string
   -p name-of-file-path
   -k file-name-range
   -v verboness
   -l subarray-size (for mpi only)
   -c cores-per-file (for mpi only)
*/

#include "queryProcessor.h"
#include <fstream>
#include <iostream>

#define DEBUG 0

static const char *options="m:M:f:F:i:I:g:G:q:Q:p:P:v:V:l:L:k:K:c:C:";
static int verboseness = 1;

static const char *condstring;
static char *datafile = 0;
static char *indexfile = 0;
static char *logfile = 0;
static char *fileModel = 0;
static char *varPath = 0;
static char *varName = 0;

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
        case 'G': logfile = optarg; break;
	case 'q':
	case 'Q': condstring = optarg; break;
	case 'p':
	case 'P': varPath = optarg; break;
	case 'v':
	case 'V': verboseness = atoi(optarg); break;
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
        case 'c':
        case 'C': groupSize = atoi(optarg); break;
	default: break;
        } // switch
    } // while
} // parseArgs

int main(int argc, char **argv) {
    std::string varPathStr;
    std::string varNameStr;
    std::ostringstream logFile("");
    parseArgs(argc, argv);
    
    if (datafile == 0 || condstring == 0) {
        std::cerr << "Usage:\n" << *argv
                  << " -f data-file [-i index-file] "
                  << " [-g log-file] [-m fileFormat]\n"
                  << "  -q query-conditions-in-a-single-string\n"
                  << " [-p varPath]"
                  << " [-k [startIdx:]endIdx] [-v verboseness]\n"
                  << " [-l mpi_subarray_size(default=100000)]"
                  << " [-c mpi_cores_per_file(defaul=evenly splitted)]\n"
                  << "\t-m -- file model [HDF5(default), H5PART, NETCDF]\n"
                  << "\t-k -- specify the start and end file index(default is single file)\n"
                  << "e.g:   ./queryIndex -f data-index.h5 -q 'px < 0.3' -p TimeStep2\n"
		  << " It processes the query on file data-index.h5\n"
                  << "e.g:   ./queryIndex -f data-index -q 'px < 0.3' -p TimeStep2 -k 10\n"
		  << " It processes the query on files data-index0, data-index1, ..., data-index9\n"
                  << "e.g:   ./queryIndex -f data-index -q 'px < 0.3' -p TimeStep2 -k 3:5\n"
		  << " It processes the query on files data-index3, data-index4, data-index5\n";
	return -1;
    }
    ibis::gParameters().add(FQ_REPORT_STATISTIC, "true");
    ibis::gParameters().add("fileManager.maxBytes", "2GB");

    FQ::FileFormat model = FQ::FQ_HDF5;
    if (fileModel != 0) {
        std::string format = fileModel;
        if (format.compare("HDF5") == 0) {
            model = FQ::FQ_HDF5;
        } else if (format.compare("H5PART") == 0) {
            model = FQ::FQ_H5Part;
        } else if (format.compare("NETCDF") == 0) {
            model = FQ::FQ_NetCDF;
        }
    }

    if (varPath != 0) {
        std::cout << "REPORT -- use variable path \"" << varPath << "\" ...\n";
        varPathStr = varPath;
    }
    if (logfile != 0) {
        std::cout << "REPORT -- use log file \"" << logFile << "\" ...\n";
	logFile << logfile;
    }

#ifndef FQ_NOMPI 
    MPI_Init(&argc, &argv);
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int group_rank = mpi_rank;
#endif

    int64_t totalHits = 0;
    ibis::horometer totTimer;
    if (endIdx == -1) {
        totTimer.start();
    	QueryProcessor queryProcessor(datafile, model, indexfile, 
		verboseness, NULL, logFile.str().c_str()); // the file handler
    	std::cout << "REPORT -- initiate the QueryProcessor object for file \"" << datafile << "\" ...\n";
	totalHits = queryProcessor.getNumHits(condstring, varPathStr, mpi_dim, mpi_len);
    } else {
#ifdef FQ_NOMPI
        totTimer.start();
        for (int i=startIdx; i<endIdx; i++) {
       	    unsigned int hits = 0;
            ibis::horometer timer;
            timer.start();
            std::ostringstream dataFileName("");;
            std::ostringstream indexFileName("");;
            dataFileName << datafile << i;
            if (indexfile != 0) {
                indexFileName << indexfile << i;
            }
    	    QueryProcessor queryProcessor(dataFileName.str(), model, indexFileName.str(), 
			verboseness, NULL, logFile.str().c_str()); // the file handler
    	    std::cout << "REPORT -- initiate the QueryProcessor object for file \"" << datafile << "\" ...\n";
            if (varPath!=0) {
		hits = queryProcessor.getNumHits(condstring, varPathStr);
            } else {
		hits = queryProcessor.getNumHits(condstring);
            }
            timer.stop();
	    totalHits += hits;
	    std::cout << "REPORT -- conditions \"" << condstring 
		<< "\" number of hits on file " << i << ":" << hits << std::endl;
        }
#else
	int64_t numHits = 0;
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
        if (logfile != 0) logFile << "/" << mpi_rank << ".log";
        totTimer.start();

        for (int fileIdx=groupId+startIdx; fileIdx<endIdx; fileIdx+=num_group) {
       	    int hits = 0;
            std::ostringstream dataFileName("");;
            std::ostringstream indexFileName("");;
            dataFileName << datafile << fileIdx;
            if (indexfile != 0) {
                indexFileName << indexfile << fileIdx;
            }
            ibis::horometer initTimer;
            ibis::horometer destructor;
            initTimer.start();
            QueryProcessor* queryProcessor;
            queryProcessor = new QueryProcessor(dataFileName.str(), model, indexFileName.str(),
                        verboseness, NULL, logFile.str().c_str(), mpi_comm); // the file handler
            MPI_Barrier(mpi_comm);
            initTimer.stop();

	    if (mpi_len > 0 ) hits = queryProcessor->getNumHits(condstring, varPathStr, mpi_dim, mpi_len);
	    else hits = queryProcessor->getNumHits(condstring, varPathStr);
            MPI_Barrier(mpi_comm);

            destructor.start();
            delete(queryProcessor);
            MPI_Barrier(mpi_comm);
            destructor.stop();
            if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC)) {
                LOGGER(true) << "Statistic\tinit\t"
                    << initTimer.CPUTime() << "\t" << initTimer.realTime();
                LOGGER(true) << "Statistic\tdestruct\t"
                    << destructor.CPUTime() << "\t" << destructor.realTime() << "\t" << hits;
            }
	    if (group_rank ==0) {
	    	std::cout << "REPORT -- conditions \"" << condstring 
		    << "\" number of hits on file " << fileIdx << ":" << hits << std::endl;
	    	numHits += hits;
	    }
        }
    	ibis::horometer idleTimer;
    	idleTimer.start();
    	MPI_Allreduce(&numHits, &totalHits, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    	idleTimer.stop();
    	MPI_Finalize();
    	if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC) && group_rank==0) {
            LOGGER(true) << "Statistic\tFinal_wait_time\t"
            	<< idleTimer.CPUTime() << "\t" << idleTimer.realTime();
    	}
    	totTimer.stop();
    	if (ibis::gParameters().isTrue(FQ_REPORT_STATISTIC) && group_rank==0) {
            LOGGER(true) << "Statistic\tTotal_time\t"
            	<< totTimer.CPUTime() << "\t" << totTimer.realTime() << "\t" << totalHits;
    	}
#endif
    }
    std::cout << "REPORT -- build indexes completed!\n";
} // main
