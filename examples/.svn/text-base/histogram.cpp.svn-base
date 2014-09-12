/**
   Print out the histogram for a specified conditions.
   If more than one time step is present in the data file, the same set of
   conditions will be applied to each time step separatedly.

   1D Histograms only in this version  **********

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
   std::set.

*/

#include "queryProcessor.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#define DEBUG 0

static const char *options="f:F:n:N:m:M:q:Q:p:P:i:I:v:V:l:L:g:G:x:X:y:Y:e:E:s:S:cCk:K:t:T:z:Z:";
static int verboseness = 1;
static char *condstring = 0;
static std::string datafile;
static std::string indexfile;
static std::ostringstream logfile;
static char *fileModel = 0;
static char *varPath = 0;
static char *varName1 = 0;
static char *varName2 = 0;
static char *varName3 = 0;
static char *hist_path = 0;
static int mpi_len = 1000;
static int mpi_dim = 0;

static int dimension = 0;	// Option x
static double begin1 = 0;		// option y
static double end1 = 0;		// option e
static double stride1 = 0;		// option s

static double begin2 = 0;		// option y
static double end2 = 0;		// option e
static double stride2 = 0;		// option s

static double begin3 = 0;       // option y
static double end3 = 0;     // option e
static double stride3 = 0;      // option s

static char *varName, *begin, *end, *stride;

static bool easyToShow = true;
static bool verification = true;

void parseArgs(int argc, char **argv) 
{
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
	switch (c) {
	case 'f':
	case 'F': datafile = optarg; break;
	case 'i':
	case 'I': indexfile = optarg; break;
	case 'g':
	case 'G': logfile << optarg; break;
	case 'n':
	case 'N': varName = optarg; break;
	case 'p':
	case 'P': varPath = optarg; break;
	case 'q':
	case 'Q': condstring = optarg; break;
	case 'v':
	case 'V': verboseness = atoi(optarg); break;
	case 'm':
	case 'M': fileModel = optarg; break;
	case 'l':
	case 'L': mpi_len = atoi(optarg);break; 
	case 'x':
	case 'X': dimension = atoi(optarg);break; 
	case 'y':
	    //case 'Y': begin1 = begin2 = begin3 = atoi(optarg);break; 
	case 'Y': begin = optarg; break; 
	case 'e':
	    //case 'E': end1 = end2 = end3 = atoi(optarg);break; 
	case 'E': end = optarg; break; 
	case 's':
	    //case 'S': stride1 = stride2 = stride3 = atof(optarg);break;
	case 'S': stride = optarg; break;
	case 'c': 
	case 'C': verification = true; break;
	case 'k':
	case 'K': hist_path = optarg; break;
	default: break;
	} // switch
    } // while
} // parseArgs

int main(int argc, char **argv) 
{
    std::string varPathStr;
    std::string varNameStr1;
    std::string varNameStr2;
    std::string varNameStr3;
    parseArgs(argc, argv);
    std::vector<double> beginList;
    std::vector<double> endList;
    std::vector<double> strideList;
    beginList.resize(dimension);
    endList.resize(dimension);
    strideList.resize(dimension);
    if (datafile.empty() || condstring == 0 || dimension == 0 ) {
	std::cerr << "Usage:\n" << *argv 
		  << " -f data-file-name" 
		  << " -q query-conditions-in-a-single-string"
		  << " -x histogram-dimension"
		  << " -y begin"
		  << " -e end"
		  << " -s stride"
		  << " [-i index-file-name]"
		  << " [-g log-file-name]"
		  << " [-n name-of-variable]"
		  << " [-p path-of-variable]" 
		  << " [-m file model [HDF5(default), H5PART, NETCDF, PNETCDF]"
		  << " [-b use-boundingbox-data-selection]"
		  << " [-v verboseness]"
		  << " [-l mpi-subarray-length]"
	    //<< "\n e.g:   ./histogram -f h5uc-data-index.h5 -q 'px < 0.3' -n y -p TimeStep2 -x 1\n"
		  <<   "\n e.g:   ./histogram -f h5uc-data.h5 -i indexfile -q 'px<0.3 && py>0' -x 2 -n py,pz;"
		  << " -y '0,-0.5;' -s '0.1,0.02;' -e '1,0;' -p TimeStep2\n\n"
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

    ibis::gParameters().add(FQ_REPORT_STATISTIC, "true");
	
    ibis::horometer totTimer;
    totTimer.start();

    FQ::FileFormat model = FQ::FQ_HDF5;
    if (fileModel != 0) {
	std::string format = fileModel;
	if (format.compare("HDF5") == 0) {
	    model = FQ::FQ_HDF5;
	} 
	else if (format.compare("H5PART") == 0) {
	    model = FQ::FQ_H5Part;
	} 
	else if (format.compare("NETCDF") == 0) {
	    model = FQ::FQ_NetCDF;
	}
	else if (format.compare("PNETCDF") == 0) {
	    model = FQ::FQ_pnetCDF;
	}
    }

    if (! indexfile.empty()) {
	if (verboseness > 1) 
	    std::cout << "DEBUG: using indexfile \"" << indexfile.c_str() << "\" ... \n";
    }

    if (varPath != 0) {
	if (verboseness > 1)  std::cout << "Debug: use variable path \"" << varPath << "\"\n";
	varPathStr = varPath;
    }
	
    //	std::cout << "varName:" << varName << " begin:" << begin << " end:" << end << " stride:" << stride << std::endl;
    varName1 = strtok(varName, ",;");
    if (dimension>1) varName2 = strtok(NULL, ",;");
    if (dimension>2) varName3 = strtok(NULL, ",;");
    begin1 = atof(strtok(begin, ",;"));
    if (dimension>1) begin2 = atof(strtok(NULL, ",;"));
    if (dimension>2) begin3 = atof(strtok(NULL, ",;"));
    end1 = atof(strtok(end, ",;"));
    if (dimension>1) end2 = atof(strtok(NULL, ",;"));
    if (dimension>2) end3 = atof(strtok(NULL, ",;"));
    stride1 = atof(strtok(stride, ",;"));
    if (dimension>1) stride2 = atof(strtok(NULL, ",;"));
    if (dimension>2) stride3 = atof(strtok(NULL, ",;"));	

    //std::cout << "varName is " << varName1 << ", " << varName2 << ", " << varName3 << std::endl;
    //std::cout << "begin   is " << begin1 << ", " << begin2 << ", " << begin3 << std::endl;
    //std::cout << "end     is " << end1 << ", " << end2 << ", " << end3 << std::endl;
    //std::cout << "stride  is " << stride1 << ", " << stride2 << ", " << stride3 << std::endl;

    if (varName1!=0) varNameStr1 = varName1;
    if (dimension>1 && varName2!=0) varNameStr2 = varName2;
    if (dimension>2 && varName3!=0) varNameStr3 = varName3;
    /*	
	if (mpi_rank==0) {
	unsigned int dims1 = static_cast<uint32_t>(1+floor((end1-begin1)/stride1));
    	unsigned int dims2 = static_cast<uint32_t>(1+floor((end2-begin2)/stride2));
    	unsigned int dims3 = static_cast<uint32_t>(1+floor((end3-begin3)/stride3));
	std::cout << "dims1 * dims2 * dims3 = " << dims1 << " * " << dims2 << " * " << dims3 << std::endl;
	}
    */


    if (logfile.str().empty() != true) {
#ifndef FQ_NOMPI
	logfile << mpi_rank << ".log";
#endif
	if (verboseness > 1) std::cout << "Debug: using logfile \"" << logfile.str().c_str() << "\"\n";
    }
    if (verboseness >1) {
	std::cout << "open the file handler" << std::endl;
    }
    // open the named file
    QueryProcessor* queryProcessor = new QueryProcessor(datafile, model, indexfile, verboseness, "", logfile.str().c_str()); // the file handler

    if (queryProcessor->isValid() == false) {
	if (verboseness > 0) {
	    std::cout << "ERROR: failed to initiate the QueryProcessor object for file \"" 
		      << datafile.c_str() << "\" ...\n";
	    std::cout << "REPORT: failed to complete processing query" << std::endl;
	}
	delete(queryProcessor);
#ifndef FQ_NOMPI
	MPI_Finalize();
#endif
	return -1;
    }

    uint64_t hits = 0;
    // getNumHits
    ibis::horometer timer;
    timer.start();
    hits = queryProcessor->getNumHits(condstring, varPathStr, mpi_dim, mpi_len);
    timer.stop();
    if (verboseness > 1)
	std::cout << "Debug: conditions \"" << condstring 
		  << "\" number of hits " << hits << std::endl;;

    if (hits == 0) {
	if (verboseness > 1) {
	    std::cout << "Warning -- No element is seleteced ==>"
		      << " the rest of the test is skipped!" << std::endl;
	}
	if (verboseness > 0) {
#ifndef FQ_NOMPI        
	    if (mpi_rank==0) {
#endif
		std::cout << "REPORT: successfully completed processing query with " 
			  << hits << " hits" << std::endl;
#ifndef FQ_NOMPI        
	    }
#endif
	}
	delete(queryProcessor);
#ifndef FQ_NOMPI
	MPI_Finalize();
#endif
	return hits;
    }

    // executeQuery
    std::vector<uint64_t> coords;
    std::vector<uint32_t> counts;
    bool herr = true;
    //	if (mpi_rank==0) std::cout<<"histogram starting..."<<std::endl;	
    if (varPath != 0) {
	//coords.reserve(hits*dims.size());
	// hits1 = queryProcessor->executeQuery((char*)condstring, coords, varPathStr, FQ::POINTS_SELECTION, mpi_dim, mpi_len);
	if (dimension==1) { 
	    //counts.assign(static_cast<uint32_t>(1+floor(end1-begin1)/stride1), 0);
	    herr = queryProcessor->get1DHistogram
		((char*) condstring, varNameStr1, varPathStr, begin1, end1, stride1, 
		 counts, mpi_dim, mpi_len);
	} else if (dimension==2) {
	    //			if (mpi_rank==0) std::cout << "in 2Dhistogram" << std::endl;
	    //counts.assign(static_cast<uint32_t>(1+floor(end1-begin1)/stride1)*
	    //          static_cast<uint32_t>(1+floor(end2-begin2)/stride2), 0);
	    herr = queryProcessor->get2DHistogram
		((char*) condstring, varPathStr, varNameStr1, begin1, end1, stride1, 
		 varNameStr2, begin2, end2, stride2, 
		 counts, mpi_dim, mpi_len);       
	    //			if (mpi_rank==0) std::cout << "out 2Dhistogram" << std::endl;
	} else if (dimension==3) {
	    herr = queryProcessor->get3DHistogram
		((char*) condstring, varPathStr, 
		 varNameStr1, begin1, end1, stride1, 
		 varNameStr2, begin2, end2, stride2, 
		 varNameStr3, begin3, end3, stride3,
		 counts, mpi_dim, mpi_len);
	}
	if (! herr) {
	    LOGGER(ibis::gVerbose >= 0)
		<< *argv << " failed to compute the histogram";
	    return -2;
	}

	/************************/
	/*  verify part         */
	/************************/

	if (verification) {

	    //verify the Histogram
	    //if (mpi_rank==0) std::cout << "starting verify the histogram..." << std::endl;
	    uint64_t len = 1;
	    if (len) {
		//std::cout<<"Warning: May use too large memory. Can only check sum.\n";				
	    } else {
		/*				double data[len];
		//
		std::vector<uint32_t> temp_counts;
		temp_counts.assign(static_cast<uint32_t>(1+floor(end-begin)/stride), 0);
		bool verr = true;
		//std::cout << "starting getData...." << std::endl;
		verr = queryProcessor->getData(varNameStr, &data[0], varPathStr);

		#ifndef FQ_NOMPI  
		if (mpi_rank==0) {
		#endif
		if (len<=1000000) {
		//std::cout << "getData success" << std::endl;		
		std::cout << "temp Histogram" << std::endl;
		    
		// copy from fasbit parth.cpp get1DHistogram
		if (len != 0) {
		for (uint32_t i = 0; i < len; ++ i) {
		++ temp_counts[static_cast<uint32_t>((data[i] - begin) / stride)];
		}
		}
		//
		std::cout << "temp Histogram" << std::endl;
		std::cout << "temp_counts.size is "<< temp_counts.size() << std::endl;
		for (int i=0; i<temp_counts.size(); i++) {
		std::cout << "[" << begin+i*stride << ", " << begin+(i+1)*stride << "]:\t" << temp_counts[i] << std::endl;
		}       
		            
		std::cout << "test Histogram" << std::endl;
		std::cout << "counts.size is "<< counts.size() << std::endl;
		for (int i=0; i<counts.size(); i++) {
		std::cout << "[" << begin+i*stride << ", " << begin+(i+1)*stride << "]:\t" << counts[i] << std::endl;
		}
		// verify two histogram vectors	
		if (counts!=temp_counts) {
		std::cout << "ERROR:Vector is not match.Histogram fail." << std::endl;
		} else {
		std::cout << "histogram success" << std::endl;
		}
		}
		#ifndef FQ_NOMPI     
		} 
		#endif
		*/			}
	    //unsigned int hits = 0 ;
	    //hits = queryProcessor->getNumHits(condstring, varPathStr, mpi_dim, mpi_len);
	    uint64_t hits1 = 0;
	    for (int i=0; i<counts.size(); i++) {
		hits1 += counts[i];
	    }
	    if (hits1 != hits) {
		std::cout<<"Error:\tcheck sum failed. Num of Hit is " << hits << ",and histogram number is " << hits1<<std::endl; 
	    } else 
		std::cout<<"verification result is correct.\n";
	}
		
	std::fstream histogramFile;
#ifndef FQ_NOMPI        
	if (mpi_rank==0) {
#endif
	    if (dimension==1) {
		std::fstream file;
		//char fileName[100]="";
		//char path[]="/global/homes/v/vidcina/fq/example/";
		//fileName<<dimension<<"D"<<"histogram["<<begin1<<":"<<stride1<<":"<<end1<<"].out";
		//sprintf(fileName, "%s%d%s%d%s%d%s%d%s", path, dimension, "Dhistogram[", begin1, ":", stride1, ":", end1, "].out");
		//std::string temp="";
		//temp.push_back(fileName.str());
		std::ostringstream fileName;
		fileName << hist_path << "_"<< dimension << "D" << "histogram["
			 << begin1 << ":" << stride1 << ":" <<end1 << "].out";
		std::string str =  fileName.str();
		const char* chr = str.c_str();
		file.open(chr, std::ios::out);
		if ( file.fail() ) {
		    std::cout << str  << std::endl;
		    std::cout << "openFile fail" << std::endl;
		} else {
		    for (int i=0; i<counts.size(); i++) {
			file << begin1+i*stride1 << "\t" << begin1+(i+1)*stride1 << "\t" << counts[i] << std::endl;
		    }
		}
		//histogramFile.close();
	    } else if (dimension==2) {
		std::cout << "2DHistogram "
			  << "Variable1 "<< varName1 << " begin " << begin1
			  << " to " << end1 <<" stride is " << stride1
			  << "Variable2 "<< varName2 << " begin " << begin2
			  << " to " << end2 <<" stride is " << stride2
			  << std::endl ;
		std::cout << "counts.size is "<< counts.size() << std::endl;
		unsigned int imax = static_cast<uint32_t>(1+floor((end1-begin1)/stride1));
		unsigned int jmax = static_cast<uint32_t>(1+floor((end2-begin2)/stride2));
				
		for (unsigned int i=0; i<imax; i++) {
		    for (unsigned int j=0; j<jmax; j++) {
			std::cout << "[" << begin1+i*stride1 << ", "
				  << begin1+(i+1)*stride1 << "), ["
				  << begin2+j*stride2 << ", "
				  << begin2+(j+1)*stride2 << "):\t"
				  << counts[i*jmax+j] << std::endl;
		    }
		}
	    } else if (dimension==3) {
		std::cout << "3DHistogram " 
			  << "Variable1 "<< varName1 << " begin " << begin1
			  << " to " << end1 << " stride is " << stride1 
			  << "Variable2 "<< varName2 << " begin " << begin2
			  << " to " << end2 << " stride is " << stride2 
			  << "Variable3 "<< varName3 << " begin " << begin3
			  << " to " << end2 << " stride is " << stride3
			  << std::endl ;

		std::cout << "counts.size is "<< counts.size() << std::endl;
		unsigned int imax = static_cast<uint32_t>
		    (1+floor((end1-begin1)/stride1));
		unsigned int jmax = static_cast<uint32_t>
		    (1+floor((end2-begin2)/stride2));
		unsigned int kmax = static_cast<uint32_t>
		    (1+floor((end3-begin3)/stride3));
#ifndef FQ_NOMPI
		if (mpi_rank==0 && imax*jmax*kmax!=counts.size()) {
		    std::cout<<"ERROR: counts.size not match."<<std::endl;
		    delete(queryProcessor);
		    MPI_Finalize();
		    return 0;
		}
#endif	
		for (unsigned int i=0; i<imax; i++) {
		    for (unsigned int j=0; j<jmax; j++) {
			for (unsigned int k=0; k<kmax; k++) {
			    if (easyToShow) {
				if (counts[i*jmax*kmax + j*kmax + k]!=0) {
				    std::cout << "[" << begin1+i*stride1 << ", "
					      << begin1+(i+1)*stride1 << "), [" 
					      << begin2+j*stride2 << ", "
					      << begin2+(j+1)*stride2 << "), [" 
					      << begin3+k*stride3 << ", "
					      << begin3+(k+1)*stride3 << "):\t" 
					      << counts[i*jmax*kmax + j*kmax + k]
					      << std::endl;
				}
			    } else {
				std::cout << "[" << begin1+i*stride1 << ", "
					  << begin1+(i+1)*stride1 << "), [" 
					  << begin2+j*stride2 << ", "
					  << begin2+(j+1)*stride2 << "), [" 
					  << begin3+k*stride3 << ", "
					  << begin3+(k+1)*stride3 << "):\t" 
					  << counts[i*jmax*kmax + j*kmax + k]
					  << std::endl;
			    }	
			}
		    }
		}
		std::cout << "successfuly printed histogram" << std::endl;	

	    }

#ifndef FQ_NOMPI        
	} 
#endif
	//		}//end else

    }//end if(!varPath)

 
    //	MPI_Barrier(MPI_COMM_WORLD);
    /*
      if (hits != hits1) 
      {
      std::cout << "Error -- number of hits does not match!" << std::endl;
      std::cout << "REPORT: failed to complete processing query" << std::endl;
      delete(queryProcessor);
      #ifndef FQ_NOMPI
      MPI_Finalize();
      #endif
      return -1;
      }
    */
   
    if (verboseness > 0) {
#ifndef FQ_NOMPI        
	if (mpi_rank==0) {
#endif
	    std::cout << "REPORT: successfully completed get1DHistogram with " 
		      << counts.size() << " histogram size" << std::endl;
#ifndef FQ_NOMPI        
	}
#endif
    }
    delete(queryProcessor);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    totTimer.stop();
    LOGGER(FastQuery::reportTiming())
	<< "Statistic\thistogram::totTimer\t"
        << totTimer.CPUTime() << "\t" << totTimer.realTime()
        << "\t";
    return hits;
} // main
