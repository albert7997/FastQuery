/*
 * A utility program that can read values from dataset or attribute from file(s).
 */

#include "ADIOS_Wrapper.h"
#include "queryProcessor.h"

#include <iostream>
#include <stdlib.h>
#include <iomanip>

static const char *options="f:F:p:P:n:N:v:V:m:M:a:A:dDg:G:";
static bool xport = false;
static std::string datafile;
static std::ostringstream logfile;
static std::string varPath;
static std::string varName;
static std::string attrName;
static int mpi_size=0;
static int mpi_rank=0;

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
        case 'g':
        case 'G': logfile << optarg; break;
        case 'p':
        case 'P': varPath = optarg; break;
        case 'n':
        case 'N': varName = optarg; break;
        case 'a':
        case 'A': attrName = optarg; break;
        case 'v':
        case 'V': ibis::gVerbose = atoi(optarg); break;
	case 'x':
	case 'X': xport = true; break;
        } // switch
    } // while

    if (! xport)
	xport = (ibis::gVerbose > 3);
} // parseArgs

template <typename E> 
bool getData(QueryProcessor &queryProcessor, E *data,  uint64_t len, 
	     std::string &variable, std::string &filename)
{
    if ( ! queryProcessor.getData(variable, &data[0])) {
    	LOGGER(ibis::gVerbose > 0)
	    << "ERROR: Failed to get the data of variable \""
	    << variable.c_str() << "\" from file \""
	    << filename.c_str() << "\"" << std::endl;
        return false;
    }
    return true;
}

template <typename E> 
bool xportData(E *data, uint64_t len, std::vector<uint64_t> dims,
		 std::string &variable) {
    ibis::util::logger lg;
    lg() << "Data of variable \"" << variable.c_str() << "\"\n{";
    // determine the length to print per line.
    int count = 10; // default is 10
    if (dims.back() < count) count = dims.back();

    for(uint64_t idx=0; idx<len; idx++) {
	if (idx % count == 0 || idx % dims.back() == 0) {
	    std::vector<uint64_t> curr;
	    curr.resize(dims.size());
	    uint64_t tmp = idx;
	    for (int i=dims.size()-1; i>=0; i--) {
		curr[i] = tmp % dims[i];
		tmp = (tmp - curr[i])/dims[i];
	    }
	    lg() << "\n    (" << curr[0];
	    for (uint64_t i=1; i<dims.size(); i++) {
		lg() << ", " << curr[i];
	    }
	    lg() << ")";
        }
	lg() << "\t" << data[idx];
    }
    lg() << "\n}";
    return true;
}


template <typename E> 
bool getAttr(QueryProcessor &queryProcessor, E *data,  uint64_t len, 
	     std::string &attrName, std::string &varName,
	     std::string &varPath, std::string &filename)
{
    if (! queryProcessor.getAttribute(varName, attrName, &data[0],
				      varPath)) {
    	LOGGER(ibis::gVerbose > 0)
    	    << "Failed to get the information for attribute \""
	    << attrName.c_str() << " of variable \""
	    << varName.c_str() << "\" from file \""
	    << filename.c_str() << "\"" << std::endl;
        return false;
    }
    if (xport) {
	ibis::util::logger lg;
    	lg() << "Value of attribute \"" << attrName.c_str() 
	     << "\" of variable \"" << varName.c_str() 
	     << "\" from file \"" << filename.c_str() << "\"";
    	for(uint64_t i=0; i<len; i++) {
	    lg() << i << " " <<  data[i];
	}
    }
    return true;
}

/// The function to do the actual read operations.  Mostly to limite the
/// scope of queryProcessor so that the files are closed before
/// adios_finalize and MPI_Finalize are called.
bool doRead() {
    FQ::DataType type;
    std::string variable;
    std::vector<uint64_t> dims;
    FQ::FileFormat model = FQ::FQ_BP;
    bool berr = true;

    // the file handler
    QueryProcessor queryProcessor(datafile, model, "", ibis::gVerbose,
				  "", logfile.str().c_str());
    if (queryProcessor.isValid() == false) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- Failed to initiate query processor for file \"" 
	    << datafile.c_str() << "\"" << std::endl;
 	return false;
    }

    if (! attrName.empty()) {
	FQ::DataType attrType; 
	uint64_t len;
	berr = queryProcessor.getAttributeInfo
	    (varName, attrName, &len, &attrType, varPath);
	if (! berr) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning Failed to get the information for attribute "
		<< attrName.c_str() << " of variable "
		<< varName.c_str() << " from file " << datafile;
	}
	else {
	    switch (attrType) {
    	    case FQ::FQT_BYTE: {
	    	char data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_FLOAT: {
	    	float data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_DOUBLE: {
	    	double data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_INT: {
	    	int32_t data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_LONG: {
	    	int64_t data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_UINT: {
	    	uint32_t data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_ULONG: {
	    	uint64_t data[len];
		berr = getAttr(queryProcessor, data, len,
			       attrName, varName, varPath, datafile);
		break;}
	    default:
		LOGGER(ibis::gVerbose > 0)
		    << "Warning -- unsupported data type for attribute \""
		    << attrName.c_str() << "\" of variable \""
		    << varName.c_str() << "\" from file \""
		    << datafile.c_str() << "\"";
	    }
	}
    }

    berr = queryProcessor.getVariableInfo
	(varName, variable, dims, &type, varPath);
    if (! berr) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- Failed to get the information for variable \""
	    << variable.c_str() << "\" from file \""
	    << datafile.c_str() << "\"" << std::endl;
    }
    else {
	uint64_t len = dims[0];
	for (int i=1; i<dims.size(); i++) len *= dims[i];
	switch (type) {
	case FQ::FQT_BYTE: {
	    ibis::array_t<char> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	case FQ::FQT_FLOAT: {
	    ibis::array_t<float> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	case FQ::FQT_DOUBLE: {
	    ibis::array_t<double> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	case FQ::FQT_INT: {
	    ibis::array_t<int32_t> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	case FQ::FQT_LONG: {
	    ibis::array_t<int64_t> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	case FQ::FQT_UINT: {
	    ibis::array_t<uint32_t> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	case FQ::FQT_ULONG: {
	    ibis::array_t<uint64_t> data(len);
	    berr = getData(queryProcessor, data.begin(), len,
			   variable, datafile);
	    if (berr && xport && mpi_rank == 0) {
		xportData(data.begin(), len, dims, variable);
	    }
	    break;}
	default:
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- unsupported data type for variable \"" 
		<< varName.c_str() << "\" from file \""
		<< datafile << "\"" << std::endl;
	    berr = false;
	}
    }
    return berr;
} // doRead

int main(int argc, char **argv) {
    parseArgs(argc, argv);
    if (datafile.empty() || varName.empty()) {
	std::cout << "Usage:\n" << *argv
		  << " -f data-file-name"
		  << " -n variable-name"
		  << " [-p variable-path]"
		  << " [-a attribute-name]"
	    //<< " [-m file-format [HDF5(default), H5PART, NETCDF]"
		  << " [-v verboseness]"
		  << " [-x (xport)]\n\n"
		  << "\tFor More detailed usage description and examples,"
	    " please see file GUIDE"
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
	<< "BPread running with mpi_size = " << mpi_size << " ...";

    if (logfile.str().empty() != true){
        if (ibis::gVerbose > 1)
	    std::cout << "DEBUG: using logfile \""
		      << logfile.str().c_str() << "\" ...\n";
        logfile << "-" << mpi_rank << ".log";
    }

    bool berr = doRead();
    adios_finalize(mpi_rank);
    LOGGER(ibis::gVerbose >= 0)
	<< *argv << " invoked adios_finalize(" << mpi_rank << ")\n";
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    if (berr) {
    	LOGGER(ibis::gVerbose > 0)
	    << *argv << " successfully complete reading data";
    	return 0;	
    }
    else {
    	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- " << *argv << " failed to complete reading data";
    	return -1;	
    }
}
