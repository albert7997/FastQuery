/*
 * A utility program that can write values to a dataset in file(s).
 */

#include "indexBuilder.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>

static const char *options="f:F:o:O:p:P:n:N:v:V:m:M:k:K:d:D:i:I:t:T:cCg:G:";
static bool check = false;
static int verboseness = 1;

static char *fileModel = 0;
static std::string datafile;
static std::ostringstream logfile;
static std::string inputfile;
static std::string varName;
static std::string varPath;
static std::string dimStr;
static std::string typeStr;

void parseArgs(int argc, char **argv) {
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
        case 'm':
        case 'M': fileModel = optarg; break;
        case 'O':
        case 'o':
        case 'f':
        case 'F': datafile = optarg; break;
	case 'i':
	case 'I': inputfile = optarg; break;
        case 'g':
        case 'G': logfile << optarg; break;
        case 'p':
        case 'P': varPath = optarg; break;
        case 'n':
        case 'N': varName = optarg; break;
        case 'v':
        case 'V': verboseness = atoi(optarg); break;
	case 'd':
	case 'D': dimStr = optarg; break;
	case 't':
	case 'Y': typeStr = optarg; break;
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
	    	std::cout << "ERROR: write data does not match to the original input" << std::endl;
	    }
	    return false;
	}
    }
    if (verboseness > 1) {
	std::cout << "DEBUG: write data is verified successfully" << std::endl;
    }
    return true;
}
int main(int argc, char **argv){

#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
    parseArgs(argc, argv);
    if (datafile.size() == 0 || inputfile.size() == 0 || varName.size() == 0
	|| varPath.size() == 0  || dimStr.size() == 0 || typeStr.size() == 0) {
        std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
		  << " -i input-data-file-name"
		  << " -n variable-name"
                  << " -p variable-path"
		  << " -d variable-dimension (e.g. 2:2)"
		  << " -t variable-type [double, float, int, long, uint, ulong, byte]"
		  << " [-m file-format [HDF5(default), H5PART, NETCDF, PNETCDF]"
		  << " [-v verboseness]"
		  << " [-c (check-written-data)]\n"
		  << "\tFor More detailed usage description and examples, please see file GUIDE"
                  << std::endl;
        return -1;
    }

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

    if (verboseness > 1) {
	std::cout << "DEBUG: initialize index builder using file " << datafile.c_str() << std::endl;
    }
    if (logfile.str().empty() != true){
        if (verboseness > 1) std::cout << "DEBUG: using logfile \"" << logfile.str().c_str() << "\" ...\n";
#ifndef FQ_NOMPI
    	logfile << "-" << mpi_rank << ".log";
#endif
    }

    IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, "", verboseness,"",logfile.str().c_str()); // the file handler
    if (indexBuilder->isValid() == false ) {
	if (verboseness > 0) {
            std::cout << "ERROR: Failed to initiate index builder for file \""
            	<< datafile.c_str() << "\"" << std::endl;
	    std::cout << "REPORT: Failed to complete writing the data" << std::endl;
        }
	delete (indexBuilder);
#ifndef FQ_NOMPI
    	MPI_Finalize();
#endif
	return -1;
    }

    std::vector<uint64_t> dims;
    int prePos = -1;
    do {
	prePos++;
    	int pos = dimStr.find(':', prePos);
	std::string str;
	if (pos != dimStr.npos) {
	    str = dimStr.substr(prePos, pos-prePos);
	} else {
	    str = dimStr.substr(prePos);
	}
	dims.push_back(atoi(str.c_str()));
	prePos = pos;
    } while(prePos != dimStr.npos); 
	
    uint64_t len = dims[0];
    std::cout << dims[0] << std::endl;
    for (int i=1; i<dims.size(); i++) {
	len *= dims[i];
	std::cout << dims[i] << std::endl;
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
	    std::cout << "REPORT: Failed to complete writing the data" << std::endl;
	}
	delete (indexBuilder);
#ifndef FQ_NOMPI
    	MPI_Finalize();
#endif
	return -1;
    }
    std::string variable = "";
    if (varPath.compare("/")==0) {
	variable = "/" + varName;
    } else {
    	variable = varPath+"/"+varName;
    }
    if ( ! indexBuilder->createNewVariable(variable, dims, type)) {
	if (verboseness > 0) {
	    std::cout << "ERROR: Failed to create new dataset for variable \"" << varName.c_str() << "\"" << std::endl;
	    std::cout << "REPORT: Failed to complete writing the data" << std::endl;
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
	    std::cout << "REPORT: Failed to complete writing the data" << std::endl;
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
	char data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERRPR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		char output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
		}
	    }
	}
	break;}
    case FQ::FQT_FLOAT: {
	float data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERROR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		float output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
		}
	    }
	}
	break;}
    case FQ::FQT_DOUBLE: {
	double data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERROR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		double output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
		}
	    }
	}
	break;}
    case FQ::FQT_INT: {
	int32_t data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	//for (uint64_t i=0; i<len; i++) infile >> data[i];
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERROR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		int32_t output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
		}
	    }
	}
	break;}
    case FQ::FQT_LONG: {
	int64_t data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERROR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		int64_t output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
		}
	    }
	}
	break;}
    case FQ::FQT_UINT: {
	uint32_t data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERROR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		uint32_t output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
		}
	    }
	}
	break;}
    case FQ::FQT_ULONG: {
	uint64_t data[len];
	for (uint64_t i=0; i<len; i++) {
	    infile >> data[i];
	    if (infile.eof()) data[i] = 0;
	}
	berr = indexBuilder->setData(varName, &data[0], varPath);
	if (! berr) {
	    std::cout << "ERROR: Failed to set the data of variable \"" << varName.c_str() 
		<< "\" to file \"" << datafile << "\"" << std::endl;
	} else {
	    if (check == true) {
		uint64_t output[len];
	    	berr = indexBuilder->getData(varName, &output[0], varPath);
		if (! berr) {
		    std::cout << "ERROR: Failed to get the data of variable \"" << varName.c_str() 
			<< "\" from file \"" << datafile << "\"" << std::endl;
		} else {
	    	    berr = verify(data, output, len);
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
            std::cout << "REPORT: Successfully completed writing the data" << std::endl;
        }
        return 0;
    } else {
        if (verboseness > 0) {
            std::cout << "REPORT: Failed to complete writing the data" << std::endl;
        }
        return -1;
    }
}
