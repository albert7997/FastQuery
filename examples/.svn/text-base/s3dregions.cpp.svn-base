#include "BPArrayIODriver.h"
#include "queryProcessor.h"
#include "indexBuilder.h"
#include <iostream>

static const char *options="d:D:f:F:i:I:l:L:n:N:o:O:p:P:v:V:w:W:b:B:";
static std::string wherestring;
static std::string datafile;
static std::string indexfile;
static std::string logfile;
static std::string outf;
static std::string varPath;
static std::string varName;
static std::vector<double> bounds;

static FQ::FileFormat model = FQ::FQ_BP;
static int mpi_rank = 0;
static int mpi_size = 1;

static void  parseArgs(int argc, char **argv) {
    extern char *optarg;
    int c;
    double tmp;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'd':
	case 'D':
	case 'f':
	case 'F': datafile = optarg; break;
	case 'i':
	case 'I': indexfile = optarg; break;
	case 'l':
	case 'L': logfile = optarg; break;
	case 'n':
	case 'N': varName = optarg; break;
	case 'o':
	case 'O': outf = optarg; break;
	case 'p':
	case 'P': varPath = optarg; break;
	case 'w':
	case 'W': wherestring = optarg; break;
	case 'v':
	case 'V': ibis::gVerbose = atoi(optarg); break;
	case 'b':
	case 'B': {
	    tmp = strtod(optarg, 0);
	    bounds.push_back(tmp);
	    break;}
	default: break;
        } // switch
    } // while
} // parseArgs

static void usage(const char *name) {
    std::cout << "Usage:\n" << *name
	      << " -d data-file-name -i index-file-name -w where-clause"
	" -l log-file-name -n name-of-variable -p path-of-variable"
	" -b bound -v verboseness\n"
	"Identify a set of regions of interests defined by 'name > bound' and "
	"additional conditions specified by the option -w.\n"
	"Note that a data file, a variable name,"
	" and one or more bounds are required\n\n e.g:\n\t"
	      << name
	      << " -d s3d1234.bp -w 'px < 0.3' -n h2o2 -b 1e-8 -b 1e-7 -b 1e-6"
	      << std::endl;
} // usage

int  main(int argc, char **argv) {
    parseArgs(argc, argv);
    if (datafile.empty() || varName.empty() || bounds.empty()) {
	usage(*argv);
	return -1;
    }

#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
    if (ibis::gVerbose >= 0 && mpi_rank == 0) {
	std::cout << *argv << " will process " << varName
		  << " from data file \"" << datafile << "\"";
	if (! wherestring.empty())
	    std::cout << "  with where condition \"" << wherestring << '"';
	if (! indexfile.empty())
	    std::cout << "  with indexfile \"" << indexfile << "\"";
	if (! varPath.empty())
	    std::cout << "  variable path \"" << varPath << "\"";
	if (! logfile.empty())
	    std::cout << "  logfile = \"" << logfile << "\"";
	if (! outf.empty())
	    std::cout << "  output file = \"" << outf << "\"";
	std::cout << std::endl;
    }

    adios_init_noxml();
    adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, 50);
    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    QueryProcessor
	qp(datafile, model, indexfile, ibis::gVerbose, "", logfile.c_str());
    if (! qp.isValid()) {
	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << "Warning: failed to initiate the QueryProcessor on file \"" 
	    << datafile << "\"";

	return -1;
    }

    int ierr = qp.recordRegions(outf, varPath, wherestring, varName, bounds);
    LOGGER(ierr < 0 && ibis::gVerbose >= 0)
	<< "Warning -- " << *argv
	<< " failed to run QueryProcessor::recordRegions";

    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    return 0;
} // main
