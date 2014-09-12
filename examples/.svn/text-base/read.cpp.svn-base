/*
 * A utility program that can read values from dataset or attribute from file(s).
 */

#include "queryProcessor.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>

static const char *options="d:D:f:F:p:P:n:N:v:V:m:M:a:A:xXg:G:";
static int verboseness = 1;
static bool xport = false;
static char *fileModel = 0;
static std::string datafile;
static std::ostringstream logfile;
static std::string varPath;
static std::string varName;
static std::string attrName;
int mpi_size, mpi_rank;

void parseArgs(int argc, char **argv) {
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
        case 'm':
        case 'M': fileModel = optarg; break;
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
        case 'V': verboseness = atoi(optarg); break;
	case 'x':
	case 'X': xport = true; break;
        } // switch
    } // while
} // parseArgs

template <typename E> 
bool getData(QueryProcessor* queryProcessor, E *data,  uint64_t len, 
	const std::string &varName, const std::string varPath, std::string &filename)
{
    if ( ! queryProcessor->getData(varName, &data[0], varPath)) {
    	if (verboseness > 0) {
	    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
		<< "\" under path \"" << varPath.c_str() << "\" from file \"" << filename.c_str() << "\"" << std::endl;
    	}
        return false;
    }
    return true;
}

template <typename E> 
bool xportData(E *data, uint64_t len, std::vector<uint64_t> dims, std::string &variable) {
    std::cout << "Data of variable \"" << variable.c_str() << "\"\n{";
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
	    std::cout << "\n    (" << curr[0];
	    for (uint64_t i=1; i<dims.size(); i++) {
		std::cout << ", " << curr[i];
	    }
	    std::cout << ")";
        }
	std::cout << "\t" << data[idx];
    }
    std::cout << "\n}" << std::endl;
    return true;
}


template <typename E> 
bool getAttr(QueryProcessor* queryProcessor, E *data,  uint64_t len, 
	std::string &attrName, std::string &varName, std::string &varPath, std::string &filename)
{
    if (! queryProcessor->getAttribute(varName, attrName, &data[0], varPath)) {
    	if (verboseness > 0) {
    	    std::cout << "ERROR: Failed to get the information for attribute \"" << attrName.c_str() 
		<< " of variable \"" << varName.c_str() << "\" from file \"" << filename.c_str() << "\"" << std::endl;
    	}
        return false;
    }
    if (xport) {
    	std::cout << "Value of attribute \"" << attrName.c_str() 
	    << "\" of variable \"" << varName.c_str() 
	    << "\" from file \"" << filename.c_str() << "\"" << std::endl;
    	for(uint64_t i=0; i<len; i++) {
	    std::cout << i << " " <<  data[i] << std::endl;
	}
    }
    return true;
}

int main(int argc, char **argv){

    parseArgs(argc, argv);
    if (datafile.size() == 0 || varName.size() == 0) {
        std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
		  << " -n variable-name"
                  << " [-p variable-path]"
		  << " [-a attribute-name]"
		  << " [-m file-format [HDF5(default), H5PART, NETCDF, PNETCDF]"
		  << " [-v verboseness]"
		  << " [-x (xport)]\n\n"
		  << "\tFor More detailed usage description and examples, please see file GUIDE"
                  << std::endl;
        return -1;
    }
#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

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

    if (logfile.str().empty() != true){
#ifndef FQ_NOMPI
        logfile << "-" << mpi_rank << ".log";
#endif
        if (verboseness > 1) 
	    std::cout << "DEBUG: using logfile \"" << logfile.str().c_str() << "\""
		<< " with verboseness degree " << verboseness << "...\n";
    }

    bool berr = true;
    QueryProcessor* queryProcessor = new QueryProcessor(datafile, model, "", verboseness, "", logfile.str().c_str()); // the file handler
    if (queryProcessor->isValid() == false) {
	if (verboseness > 0) {
	    std::cout << "ERROR: Failed to initiate query processor for file \"" 
		<< datafile.c_str() << "\"" << std::endl;
	}
 	berr = false;
    }

    if ( attrName.empty()) {
	std::string variable;
	std::vector<uint64_t> dims;
  	FQ::DataType type;
	if (! queryProcessor->getVariableInfo(varName, variable, dims, &type, varPath)) {
	    if (verboseness > 0) {
		std::cout << "ERROR: Failed to get the information for variable \"" << variable.c_str() 
			<< "\" from file \"" << datafile.c_str() << "\"" << std::endl;
	    }
	    berr = false;
 	} else {
	    uint64_t len = dims[0];
	    for (int i=1; i<dims.size(); i++) len *= dims[i];
	    switch (type) {
    	    case FQ::FQT_BYTE: {
	    	char data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
    	    case FQ::FQT_FLOAT: {
	    	float data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
    	    case FQ::FQT_DOUBLE: {
	    	double data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
    	    case FQ::FQT_INT: {
	    	int32_t data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
    	    case FQ::FQT_LONG: {
	    	int64_t data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
    	    case FQ::FQT_UINT: {
	    	uint32_t data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
    	    case FQ::FQT_ULONG: {
	    	uint64_t data[len];
		berr = berr && getData(queryProcessor, data, len, varName, varPath, datafile);
	  	if (xport) {
#ifndef FQ_NOMPI
		if (mpi_rank == 0) {
#endif
		    xportData(data, len, dims, variable);
#ifndef FQ_NOMPI
		}
#endif
		}
		break;}
	    default:
		LOGGER(verboseness > 0)
		    << "ERROR: Unknown or unsupported data type for variable \"" 
		    << varName.c_str() << "\" from file \"" << datafile << "\"";
		berr = false;
	    }
	}
    } else {
	FQ::DataType attrType; 
	uint64_t len;
	if (! queryProcessor->getAttributeInfo(varName, attrName, &len, &attrType, varPath)) {
	    std::cout << "ERROR: Failed to get the information for attribute " << attrName.c_str() 
		<< " of variable " << varName.c_str() << " from file " << datafile << std::endl;
	    berr = false;
	} else {
	    switch (attrType) {
    	    case FQ::FQT_BYTE: {
	    	char data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_FLOAT: {
	    	float data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_DOUBLE: {
	    	double data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_INT: {
	    	int32_t data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_LONG: {
	    	int64_t data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_UINT: {
	    	uint32_t data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
    	    case FQ::FQT_ULONG: {
	    	uint64_t data[len];
		berr = berr && getAttr(queryProcessor, data, len, attrName, varName, varPath, datafile);
		break;}
	    default:
		if (verboseness > 0) {
		    std::cout << "ERROR: Unknown or unsupported data type for attribute \"" 
			<< attrName.c_str() << "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile.c_str() << "\"" << std::endl;
		}
		berr = false;
	    }
	}
    }
    delete(queryProcessor);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    if (berr) {
    	if (verboseness > 0) {
	    std::cout << "REPORT: Successfully completed reading data" << std::endl;
	}
    	return 0;	
    } else {
    	if (verboseness > 0) {
	    std::cout << "REPORT: Failed to complete reading data" << std::endl;
	}
    	return -1;	
    }
}
