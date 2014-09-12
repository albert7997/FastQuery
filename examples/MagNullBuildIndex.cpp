/**
   A simple program to test the capability of building and storing
   indexes to a file.

   Binning option is specified in the form of "<binning ... />".  NOTE that
   on most systems, the binning option needs to be quoted because it
   involved characters that have special meaning to commonly used shells.
*/

#include "BPArrayIODriver.h"
#include "indexBuilder.h"
#include <iostream>	

static const char *options="b:B:f:F:p:P:n:N:i:I:v:V:m:M:l:L:g:G:rR:t:T:";

static std::string datafile;
static std::string indexfile;
static std::ostringstream logfile;
static char *binning = 0;
static bool forcerebuild = false;
static char* varPath = 0;
static char* varName = 0;
static int mpi_len = 10000;
static int mpi_dim = 0;
static double threshold = 0.1;
static int mpi_size=0, mpi_rank=0;

static std::string varPathStr;
static std::string varNameStr;
static std::string grpName = "record";
static std::string new_fileName;
static std::string new_grpName = grpName;
static std::string new_varName = "/var/B";

static int firstTimestep = 0;
static int lastTimestep = 0;

static FQ::FileFormat model = FQ::FQ_BP;
static ADIOS_File* adiosDataFile = 0;
static IndexBuilder* indexBuilder = 0;

// forward declaration
static int read_data();
static int build_index();


static void 
parseArgs(int argc, char **argv) 
{
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	    //case 'm':
	    //case 'M': fileModel = optarg; break;
	case 'f':
	case 'F': datafile = optarg; break;
	case 'i':
	case 'I': indexfile = optarg; break;
	case 'g':
	case 'G': logfile << optarg; break;
	    //case 'p':
	    //case 'P': varPath = optarg; break;
	    //case 'n':
	    //case 'N': varName = optarg; break;
	case 'v':
	case 'V': ibis::gVerbose = atoi(optarg); break;
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
	    //<< " [-n variable-name] [-p variable-path]"
		  << " [-b '<binning nbins=1000 />' (default unbinned)]"
		  << " [-r (force-rebuild-index)] [-v verboseness]"
	    //<< " [-m fileModel [HDF5(default), H5PART, NETCDF]]"
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

    adios_init_noxml ();
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 1000);

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	<< "MagNullBuildIndex is running with mpi_size = " << mpi_size << "...";

    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    ibis::gParameters().add("fileManager.maxBytes", "6GB");
    if (forcerebuild) {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "true");
    } else {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "false");
    }

    if (! indexfile.empty()) {
	LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: using indexfile \"" << indexfile.c_str() << "\" ...";
    } else {
        std::ostringstream oss;
        oss << datafile;
        long pos = oss.tellp();
        oss.seekp(pos-3);
        oss << "_idx.bp";
        indexfile = oss.str();
	LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: using auto indexfile \"" << indexfile.c_str() << "\" ...";
    }

    if (binning !=0) {
      	LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: using binning option \"" << binning << "\" ...";
    }

    if (varPath !=0) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: specifying variable path \"" << varPath << "\" ...";
	varPathStr = varPath;
    }

    if (varName !=0) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: specifying variable name \"" << varName << "\" ...";
	varNameStr = varName;
    }

    if (logfile.str().empty() != true){
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: using logfile \"" << logfile.str().c_str() << "\" ...";
    	logfile << "-" << mpi_rank << ".log";
    }

    if (mpi_len != 10000 && mpi_size == 1) {
        LOGGER(ibis::gVerbose > 1 && mpi_rank == 0)
	    << "DEBUG: ignore mpi_len for a single process ...";
        mpi_len = 10000;
    }

    LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	<< "Threshold=" << threshold << ", Timestep = " << lastTimestep;

    std::ostringstream oss;
    oss << datafile;
    long pos = oss.tellp();
    oss.seekp(pos-3);
    oss << "_add.bp";
    
    new_fileName = oss.str();
    LOGGER(ibis::gVerbose >= 0)
	<< "file=" << new_fileName << ", group=" << new_grpName << ", variable=" << new_varName;

#if 0 // why it doesn't work properly with multiprocesses???
    if (mpi_rank==0 && ADIOS_File::exists(new_fileName)) {
        LOGGER(ibis::gVerbose >= 0)
	    << "File " << new_fileName << " already exists --- remove!";
        remove(new_fileName.c_str());
    }
    if (mpi_rank==0 && !indexfile.empty() && ADIOS_File::exists(indexfile)) {
        LOGGER(ibis::gVerbose >= 0)
	    << "File " << indexfile << " already exists --- remove!";
        remove(indexfile.c_str());
    }
#endif

    ibis::horometer timer;
    timer.start();
    ret = read_data();
    if (ret < 0) {
        LOGGER(ibis::gVerbose > 0)
	    << "DEBUG: read_data failed!";
    }
    else {
        ret = build_index();
        LOGGER(ret < 0 && ibis::gVerbose > 0)
	    << "DEBUG: build_index failed!";
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
    timeArray[IDX_READ_DATA_TIME] = BPArrayIODriver::getReadDataTimer().realTime();
    timeArray[IDX_WRITE_INDEX_TIME] = BPArrayIODriver::getWriteIndexTimer().realTime();
#ifndef FQ_NOMPI
    double timeArrayArray[idx*mpi_size];
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(timeArray, idx, MPI_DOUBLE, timeArrayArray, idx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

	printf("MPI size=%d, Mean Total=%.3f, CPU=%.3f, ReadData=%.3f, "
	       "WriteIndex=%.3f, Median Total=%.3f, CPU=%.3f, "
	       "ReadData=%.3f, WriteIndex=%.3f\n", mpi_size,
	       sumTimeArray[IDX_TOT_TIME]/mpi_size,
	       sumTimeArray[IDX_CPU_TIME]/mpi_size,
	       sumTimeArray[IDX_READ_DATA_TIME]/mpi_size,
	       sumTimeArray[IDX_WRITE_INDEX_TIME]/mpi_size,
	       util::compute_median(totTimeArray, mpi_size),
	       util::compute_median(cpuTimeArray, mpi_size),
	       util::compute_median(readDataTimeArray, mpi_size),
	       util::compute_median(writeIndexTimeArray, mpi_size));
    }

    if (adiosDataFile)
        delete adiosDataFile;
    if (indexBuilder)
        delete indexBuilder;

    adios_finalize(mpi_rank);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    return ret;
} // main


static int
read_data()
{
    adiosDataFile = new ADIOS_File(datafile);

    if (adiosDataFile == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to new ADIOS_File() \""
	    << datafile.c_str() << "\"";
	return -1;
    }
    if (0 != adiosDataFile->getHandle()) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to open the file \""
	    << datafile << "\"";
	return -2;
    }

    firstTimestep = 0;

    ADIOS_Var* v5 = adiosDataFile->getVariable("/var/v5");
    ADIOS_Var* v6 = adiosDataFile->getVariable("/var/v6");
    ADIOS_Var* v7 = adiosDataFile->getVariable("/var/v7");
    if (v5==0 || v6==0 || v7==0) {
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: No such variable v5,v6,v7 in record";
	return -3;
    }

    std::vector<uint64_t> dims(v5->getNumDimension());
    v5->getDimensionArray(&dims);
    FQ::DataType type = BPArrayIODriver::getType(v5->getType());

    // for collective writes, we choose smaller mpi_len than dims[1]/mpi_size
    // the last process takes some more if dims[1]%mpi_size != 0
    int new_mpi_len = dims[1]/mpi_size;

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
    start[1] = new_mpi_len*mpi_rank;
#if 0
    if (start[1] >= dims[1]) {
	LOGGER(ibis::gVerbose > 0)
	    << "Debug: the process is waivered for readData: " 
	    << "start[1]=" << start[1] << " and "
	    << "dims[1]=" << dims[1];
	return 1;
    }
    else if (start[1]+mpi_len <= dims[1])
	count[1] = mpi_len;
    else
	count[1] = dims[1]-start[1];

#else
    if (mpi_rank == mpi_size-1) {
	count[1] = dims[1]-start[1];
    }
    else {
	count[1] = new_mpi_len;
    }
    LOGGER(ibis::gVerbose > 2)
	<<"[" << mpi_rank << "] takes count[1]=" << count[1];
#endif

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
    double* new_data = new double[nlocalelem];

    //
    // define variable
    //

    const char* path = "";
    const char* option = "a";	// append 

    //int64_t adios_handle;
    uint64_t adios_groupsize=0;
    uint64_t adios_totalsize=0;
    //int adios_err;

    int64_t adios_file;
    int64_t adios_group;
    MPI_Comm mycomm = MPI_COMM_WORLD;

    adios_declare_group(&adios_group, new_grpName.c_str(), "iter",
			adios_flag_yes);
    adios_select_method(adios_group, "MPI", "", "");

    std::ostringstream oss, oss_dims, oss_globals, oss_offsets;
#if 0
    oss << new_varName << "_dim1";
    oss_globals << oss.str().c_str()+1;

    adios_define_var(adios_group, oss.str().c_str()+1, path,
		     adios_unsigned_integer, 0, 0, 0);
    adios_groupsize += sizeof(uint32_t);
    for (int i=2; i<dims.size(); i++) {
	std::ostringstream oss2;
	oss2 << new_varName << "_dim" << i;
	oss_globals << "," << oss2.str().c_str()+1;

	adios_define_var(adios_group, oss2.str().c_str()+1, path,
			 adios_unsigned_integer, 0, 0, 0);
	adios_groupsize += sizeof(uint32_t);
    }

    oss_dims << "iter," << oss_globals.str();
    oss_offsets << mpi_len*mpi_rank << ",0,0";
#else
    oss_dims << "iter," << count[1];
    oss_globals << dims[1];
    oss_offsets << new_mpi_len*mpi_rank;
    for (int i=2; i<dims.size(); i++) {
	oss_dims << "," << count[i];
	oss_globals << "," << count[i];
	oss_offsets << ",0";
    }
#endif
	
    LOGGER(ibis::gVerbose > 0)
	<< "adios_define_var is called for " << new_varName
	<< ", with dimensions=\"" << oss_dims.str()
	<< "\", globals=\"" << oss_globals.str() << "\", offsets=\""
	<< oss_offsets.str() << "\"";
    switch(type){
    case FQ::FQT_BYTE:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_byte, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(char);
	break;}
    case FQ::FQT_INT:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_integer, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(int32_t);
	break;}
    case FQ::FQT_UINT:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_unsigned_integer, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(uint32_t);
	break;}
    case FQ::FQT_LONG:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_long, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(int64_t);
	break;}
    case FQ::FQT_ULONG:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_unsigned_long, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(uint64_t);
	break;}
    case FQ::FQT_FLOAT:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_real, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(float);
	break;}
    case FQ::FQT_DOUBLE:{
	adios_define_var(adios_group, new_varName.c_str()+1, path,
			 adios_double, oss_dims.str().c_str(),
			 oss_globals.str().c_str(), oss_offsets.str().c_str());
	adios_groupsize += nelem*sizeof(double);
	break;}
    default:
	LOGGER(ibis::gVerbose > 2)
	    << "Debug -- read_data():"
	    << " unknown FQ data type " << type;
	return -4;
    }

    for (unsigned time=firstTimestep; time<=lastTimestep; time++) {
	start[0] = time-firstTimestep;

	int64_t nbytes = v5->readData(v5_data, -1, start, count);
	if (nbytes <= 0) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- Can't get data for v5";
	    return -5;
	}

	nbytes = v6->readData(v6_data, -1, start, count);
	if (nbytes <= 0) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- Can't get data for v6";
	    return -6;
	}

	nbytes = v7->readData(v7_data, -1, start, count);
	if (nbytes <= 0) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- Can't get data for v7";
	    return -7;
	}

	//int cnt=0;
	for (int k=0; k<nlocalelem; k++) {
	    new_data[k] = sqrt(v5_data[k]*v5_data[k] + v6_data[k]*v6_data[k]
			       + v7_data[k]*v7_data[k]);
	    //if (new_data[k] < threshold) ++cnt;
	}

	// Write index to the file!
	adios_open(&adios_file, new_grpName.c_str(), new_fileName.c_str(),
		   option, &mycomm);
	adios_group_size(adios_file, adios_groupsize, &adios_totalsize);

#if 0
	for (int i=1; i<dims.size(); i++) {
	    std::ostringstream oss2;
	    oss2 << new_varName << "_dim" << i;
	    ret = adios_write(adios_file, oss2.str().c_str(), (void*)&(dims[i]));
	    LOGGER(ibis::gVerbose > 0)
		<< "Writing variable " << oss2.str() << ": val = " << dims[i]
		<< ", ret = " << ret << "...";
	}
#endif

	(void) adios_write(adios_file, new_varName.c_str(), (void*)new_data);
	adios_close(adios_file);

	LOGGER(ibis::gVerbose > 0 && mpi_rank == 0)
	    << "read_data: timestep " << time << " is processed...";
    }

    delete[] v5_data;
    delete[] v6_data;
    delete[] v7_data;
    delete[] new_data;

    LOGGER(ibis::gVerbose > 5)
	<< "Debug -- " << "read_data():"
	<< " group size is " << adios_groupsize;

    return 0;
}

static int 
build_index()
{
    ibis::horometer indexTimer;
    indexTimer.start();

    bool berr = true;
    indexBuilder = new IndexBuilder
	(new_fileName, model, indexfile, ibis::gVerbose, NULL,
	 logfile.str().c_str()); // the file handler
    if (indexBuilder->isValid() == false) {
	berr = false;
	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to initiate the IndexBuilder object for file \"" 
	    << new_fileName.c_str() << "\"";
	return -1;
    } else {
	LOGGER(ibis::gVerbose > 1)
	    << "DEBUG: initiate the IndexBuilder object for file \"" 
	    << new_fileName.c_str() << "\"";
    }
    LOGGER(ibis::gVerbose > 0)
	<< (berr == false ? 
	    "Warning: Failed to complete building index" :
	    "REPORT: Successfully complete building index");

    for (unsigned time=firstTimestep; time<=lastTimestep; time++) {
	std::ostringstream oss;
	oss << "/time" << time << "/record/var/B";
	varNameStr = oss.str(); 
	std::vector<std::string> variables;
	int numVar = indexBuilder->getAllVariables
	    (variables, varPathStr, varNameStr);
	if (numVar <= 0) {
	    LOGGER(ibis::gVerbose >= 0)
		<< "[" << time << "] numVar <= 0";
	    return -2;
	}

	int builtVar = 0;
	if (mpi_len > 0) {
	    builtVar = indexBuilder->buildIndexes
		(binning, varPathStr, varNameStr, mpi_dim, mpi_len);
	} else {
	    if (varName!=0) {
		builtVar = indexBuilder->buildIndexes
		    (binning, varPathStr, varNameStr);
	    } else if (varPath!=0) {
		builtVar = indexBuilder->buildIndexes(binning, varPathStr);
	    } else {
		builtVar = indexBuilder->buildIndexes(binning);
	    }
	}
	LOGGER(builtVar != numVar && ibis::gVerbose > 0)
	    << "Warning: buildIndexes expected a return valeu of " << numVar
	    << ", but the actual return value is " << builtVar;

	LOGGER(ibis::gVerbose >= 0 && mpi_rank == 0)
	    << "build_index: timestep " << time << " completed";
    }

    return 0;
}
