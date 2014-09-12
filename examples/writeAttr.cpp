/*
 * A utility program that can write attribute values to a dataset in file(s).
 */

#include "indexBuilder.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>

static const char *options="f:F:n:N:p:P:a:A:v:V:m:M:k:K:i:I:t:T:l:L:cCg:G:";
static bool check = false;
static int verboseness = 1;

static char *fileModel = 0;
static std::string datafile;
static std::string inputfile;
static std::ostringstream logfile;
static std::string varName;
static std::string pathName;
static std::string attrName;
static std::string typeStr;
static int attrLen;

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
	case 'I': inputfile = optarg; break;
        case 'g':
        case 'G': logfile << optarg; break;
        case 'p':
        case 'P': pathName = optarg; break;
        case 'a':
        case 'A': attrName = optarg; break;
        case 'n':
        case 'N': varName = optarg; break;
        case 'v':
        case 'V': verboseness = atoi(optarg); break;
	case 't':
	case 'Y': typeStr = optarg; break;
	case 'l':
	case 'L': attrLen = atoi(optarg); break;
	case 'c':
	case 'C': check =true; break;
        } // switch
    } // while
} // parseArgs

template <typename E> 
bool verify(E *input, E *data, uint64_t len)
{
    for(uint64_t i=0; i<len; i++) {
	if (data[i] != input[i]) {
    	    if (verboseness > 0) {
	    	std::cout << "ERROR: write attribute does not match to the original input" << std::endl;
	    }
	    return false;
	}
    }
    if (verboseness > 1) {
	std::cout << "DEBUG: write attribute is verified successfully" << std::endl;
    }
    return true;
}
int main(int argc, char **argv){

    parseArgs(argc, argv);
    if (datafile.size() == 0 || inputfile.size() == 0 || varName.size() == 0 
	|| attrName.size() == 0 || typeStr.size() == 0 || attrLen <= 0) {
        std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
		  << " -i input-file-name"
		  << " -g log-file-name"
		  << " -n variable-name"
                  << " -p variable-path"
		  << " -a attribute-name"
		  << " -t attribute-type [double, float, int, long, uint, ulong, byte]"
		  << " -l attribute-length"
		  << " [-m file-format [HDF5(default), H5PART, NETCDF, PNETCDF]"
		  << " [-v verboseness]"
		  << " [-c (check-written-data)]\n\n"
		  << "\tFor More detailed usage description and examples, please see file GUIDE"
                  << std::endl;
        return -1;
    }

#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

    FQ::FileFormat model = FQ::FQ_HDF5;
    if (fileModel != 0) {
	std::string format = fileModel;
    	if (format.compare("HDF5") == 0) {
	    model = FQ::FQ_HDF5;
	    if (verboseness > 1) {
		std::cout << "DEBUG: use HDF5 data format" << std::endl;
	    }
    	} else if (format.compare("H5PART") == 0) {
	    model = FQ::FQ_H5Part;
	    if (verboseness > 1) {
		std::cout << "DEBUG: use H5PART data format" << std::endl;
	    }
    	} else if (format.compare("NETCDF") == 0) {
	    model = FQ::FQ_NetCDF;
	    if (verboseness > 1) {
		std::cout << "DEBUG: use NETCDF data format" << std::endl;
	    }
    	} else if (format.compare("PNETCDF") == 0) {
	    model = FQ::FQ_pnetCDF;
	    if (verboseness > 1) {
		std::cout << "DEBUG: use PNETCDF data format" << std::endl;
	    }
    	}
    }

    if (logfile.str().empty() != true){
        if (verboseness > 1) std::cout << "DEBUG: using logfile \"" << logfile.str().c_str() << "\" ...\n";
#ifndef FQ_NOMPI
        logfile << "-" << mpi_rank << ".log";
#endif
    }
    if (verboseness > 1) {
	std::cout << "DEBUG: initialize index builder using file " << datafile.c_str() << std::endl;
    }
    IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, "", verboseness,"",logfile.str().c_str()); // the file handler
    if (indexBuilder->isValid() == false ) {
	if (verboseness > 0) {
            std::cout << "ERROR: Failed to initiate index builder for file \""
            	<< datafile.c_str() << "\"" << std::endl;
	    std::cout << "REPORT: Failed to complete writing the attribute" << std::endl;
        }
	delete (indexBuilder);
#ifndef FQ_NOMPI
        MPI_Finalize();
#endif
	return -1;
    }

    FQ::DataType type = FQ::FQT_UNKNOWN;
    if (typeStr.compare("char") == 0) {
	type = FQ::FQT_BYTE;
    } else if (typeStr.compare("double") == 0) {
	type = FQ::FQT_DOUBLE;
    } else if (typeStr.compare("float") == 0) {
	type = FQ::FQT_FLOAT;
    } else if (typeStr.compare("int") == 0) {
	type = FQ::FQT_INT;
    } else if (typeStr.compare("long") == 0) {
	type = FQ::FQT_LONG;
    } else if (typeStr.compare("uint") == 0) {
	type = FQ::FQT_UINT;
    } else if (typeStr.compare("ulong") == 0) {
	type = FQ::FQT_ULONG;
    } else {
	if (verboseness > 0) {
	    std::cout << "ERROR: Invalid data type. Must be char, double, float, int, long, uint or ulong." << std::endl;
	    std::cout << "REPORT: Failed to complete writing the attribute" << std::endl;
	}
	delete (indexBuilder);
#ifndef FQ_NOMPI
        MPI_Finalize();
#endif
	return -1;
    }


    std::ifstream infile(inputfile.c_str());
    if (! infile) {
        if (verboseness > 0) {
            std::cout << "ERROR: Failed to open the input file \"" << inputfile.c_str() << "\"" << std::endl;
            std::cout << "REPORT: Failed to complete writing the attribute" << std::endl;
        }
	delete (indexBuilder);
#ifndef FQ_NOMPI
        MPI_Finalize();
#endif
        return -1;
    }

    bool berr = true;
    switch (type) {
    case FQ::FQT_BYTE: {
	char data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		char output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    case FQ::FQT_FLOAT: {
	float data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		float output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    case FQ::FQT_DOUBLE: {
	double data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		double output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    case FQ::FQT_INT: {
	int32_t data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		int32_t output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    case FQ::FQT_LONG: {
	int64_t data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		int64_t output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    case FQ::FQT_UINT: {
	uint32_t data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		uint32_t output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    case FQ::FQT_ULONG: {
	uint64_t data[attrLen];
	for (uint64_t i=0; i<attrLen; i++) infile >> data[i];
	berr = indexBuilder->setAttribute(varName, attrName, &data[0], attrLen, type, pathName);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the attribute \"" << attrName.c_str() 
		<< "\"of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		uint64_t output[attrLen];
	    	berr = indexBuilder->getAttribute(varName, attrName, &output[0], pathName);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the attribute \"" << attrName.c_str()
			<< "\" of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, attrLen);
		}
	    }
	}
	break;}
    default: {
	if (verboseness > 0) {
	    std::cout << "ERROR: Unknown or unsupported data type for variable \"" << varName.c_str()
                << "\" from file \"" << datafile << "\"" << std::endl;
        }
        berr = false;
	break;}
    }
    delete (indexBuilder);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    if (berr) {
        if (verboseness > 0) {
            std::cout << "REPORT: Successfully completed writing the attribute" << std::endl;
        }
        return 0;
    } else {
        if (verboseness > 0) {
            std::cout << "REPORT: Failed to complete writing the attribute" << std::endl;
        }
        return -1;
    }
}
