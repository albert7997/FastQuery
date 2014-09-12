#include "fq.h"

#ifdef FQ_HAVE_HDF5
#include "hdf5file.h"
#endif

#ifdef FQ_HAVE_NETCDF
#include "netCDFfile.h"
#endif

#ifdef FQ_HAVE_PNETCDF
#include "pnetCDFfile.h"
#endif

#ifdef FQ_HAVE_BP
#include "BPArrayIODriver.h"
#endif

#include <cmath>        // std::ceil

// The default value of report_timing is false.
bool FastQuery::report_timing = false;

/*******************************************
 * Constructor and Destructor
 ********************************************/
#ifdef FQ_NOMPI
FastQuery::FastQuery(const std::string& dataFileName,
                     const FQ::FileFormat ffmt,
                     const std::string& indexFileName,
                     const int v, const char *rcfile,
                     const char *logfile,
                     bool readOnly,
                     void *extra)
#else
FastQuery::FastQuery(const std::string& dataFileName,
                     const FQ::FileFormat ffmt,
                     const std::string& indexFileName,
                     const int v, const char *rcfile,
                     const char *logfile,
                     bool readOnly,
                     const MPI_Comm comm,
                     void *extra)
#endif
{
    ibis::util::setVerboseLevel(v);
    ibis::init(rcfile, logfile);

    dataFile = 0;
    indexFile = 0;
    metadataMgr = 0;
    // true if indexFile needs to be initiated.
    bool indexing = (indexFileName.compare("") != 0 &&
                     indexFileName.compare(dataFileName) != 0);

#ifndef FQ_NOMPI
    mpi_comm = comm;
    MPI_Group mpi_group;
    MPI_Comm_group(mpi_comm, &mpi_group);
    MPI_Group_size(mpi_group, &mpi_size);
    MPI_Group_rank(mpi_group, &mpi_rank);
#endif

    LOGGER(ibis::gVerbose > 0)
        << "FastQuery constructor invoked with datafileName=" << dataFileName
        << ", fileFormat=" << ffmt << ", readOnly=" << readOnly;

    // open the file
    std::string indexPath = "";
    switch (ffmt) {

#ifdef FQ_HAVE_HDF5
    case FQ::FQ_H5Part: indexPath = "/__H5PartIndex__";
    case FQ::FQ_HDF5: {
#ifdef FQ_NOMPI
        if (readOnly == false && indexing == true) {
            dataFile = new HDF5(dataFileName, true, indexPath);
        }
        else {
            dataFile = new HDF5(dataFileName, readOnly, indexPath);
        }
#else
        if (readOnly == false && indexing == true) {
            dataFile = new HDF5(dataFileName, mpi_comm, true, indexPath);
        }
        else {
            dataFile = new HDF5(dataFileName, mpi_comm, readOnly,
                                indexPath);
        }
#endif
        if (dataFile->isValid() == false) {
            delete (dataFile);
            dataFile = 0;
        }
        if (indexing) {
#ifdef FQ_NOMPI
            indexFile = new HDF5(indexFileName, readOnly, indexPath);
#else
            indexFile = new HDF5(indexFileName, mpi_comm, readOnly,
                                 indexPath);
#endif
            if (indexFile->isValid() == false) {
                delete (indexFile);
                indexFile = 0;
            }
        }
        else {
            indexFile = dataFile;
        }
        break;}
#endif


#ifdef FQ_HAVE_NETCDF
    case FQ::FQ_NetCDF: {
#ifndef FQ_NOMPI
        // netcdf is not supported in MPI mode yet
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::FastQuery:"
            << " MPI is not supported for NetCDF yet";
#else
        if (readOnly == false && indexing == true) {
            dataFile = new NETCDF(dataFileName, true, "");
        }
        else {
            dataFile = new NETCDF(dataFileName, readOnly, "");
        }
        if (dataFile->isValid() == false) {
            delete (dataFile);
            dataFile = 0;
        }
        if (indexing) {
            indexFile = new NETCDF(indexFileName, readOnly, "");
            if (indexFile->isValid() == false) {
                delete (dataFile);
                dataFile = 0;
                delete (indexFile);
                indexFile = 0;
            }
        }
        else {
            indexFile = dataFile;
        }
#endif
        break;}
#endif


#ifdef FQ_HAVE_PNETCDF
    case FQ::FQ_pnetCDF: {
        if (readOnly == false && indexing == true) {
            dataFile = new PNETCDF(dataFileName, true, "", mpi_comm);
        } else {
            dataFile = new PNETCDF(dataFileName, readOnly, "", mpi_comm);
        }
        if (dataFile->isValid() == false) {
            delete (dataFile);
            dataFile = 0;
        }
        if (indexing) {
            indexFile = new PNETCDF(indexFileName, readOnly, "", mpi_comm);
            if (indexFile->isValid() == false) {
                delete (dataFile);
                dataFile = 0;
                delete (indexFile);
                indexFile = 0;
            }
        } else {
            indexFile = dataFile;
        }
        break;}
#endif


#ifdef FQ_HAVE_BP
    case FQ::FQ_BP: {
        MPI_Comm comm = MPI_COMM_WORLD;
#ifndef FQ_NOMPI
        comm = mpi_comm;
#endif
        bool streaming = true;
        float timeout = 0.0;
        enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_BP;
        if (extra != 0) {
            const BPExtras &bpx = *static_cast<BPExtras*>(extra);
            streaming = bpx.streaming;
            read_method = bpx.read_method;
            if (bpx.timeout != 0.0)
                timeout = bpx.timeout;
            else if (read_method == ADIOS_READ_METHOD_DATASPACES ||
                     read_method == ADIOS_READ_METHOD_DIMES ||
                     read_method == ADIOS_READ_METHOD_FLEXIO)
                timeout = FQ_ADIOS_STREAM_TIMEOUT;
        }

        if (readOnly == false && indexing == true) {
            dataFile = new BPArrayIODriver
                (dataFileName, "", comm, read_method, timeout, streaming);
        }
        else {
            dataFile = new BPArrayIODriver
                (dataFileName, "", comm, read_method, timeout, streaming);
        }
        if (dataFile->isValid() == false) {
            delete (dataFile);
            dataFile = 0;
        }
        if (indexing) {
            indexFile = new BPArrayIODriver
                (indexFileName, "", comm, read_method, timeout, streaming);
            if (indexFile->isValid() == false) {
                delete (dataFile);
                dataFile = 0;
                delete (indexFile);
                indexFile = 0;
            }
        }
        else {
            indexFile = dataFile;
        }
        break;}
#endif
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::FastQuery: unsupport file model";
        break;}
    }

    if (dataFile == 0 || indexFile == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery:"
            << " failed to initialize the FastQuery object";
        return;
    }
#if defined(DEBUG) && DEBUG+0 > 0
    report_timing = true;
#else
#ifdef FQ_REPORT_STATISTIC
    if (! report_timing)
        report_timing = ibis::gParameters().isTrue(FQ_REPORT_STATISTIC);
#endif
    if (ibis::gVerbose < 0) {
        report_timing = false;
    }
    else if (! report_timing) {
        if (ibis::gVerbose > 3) {
            report_timing = true;
        }
        else {
            report_timing =
                ibis::gParameters().isTrue("FastQuery.reportTiming");
        }
    }
#endif

#ifndef FQ_NOMPI
    if (report_timing)
        report_timing = (ibis::gVerbose > 5 || mpi_rank == 0);
#endif

    ibis::horometer timer;
    if (report_timing)
        timer.start();
    // initialize information manager
    metadataMgr = new MetadataMgr(*dataFile);
    if (report_timing) {
        timer.stop();
        LOGGER(true) << "Statistic\tFQ::init\t"
                     << timer.CPUTime() << "\t" << timer.realTime();
    }
    LOGGER(ibis::gVerbose > 2)
        << "FastQuery: successfully initialized the FastQuery object";
} // FastQuery::FastQuery

FastQuery::~FastQuery(){
#ifdef DEBUG
    LOGGER(ibis::gVerbose > 2)
        << "FastQuery::~FastQuery: destructor is called";
#endif
    // close the file
    if (indexFile != dataFile) {
        if ( indexFile) {
            delete indexFile;
#ifdef DEBUG
            LOGGER(ibis::gVerbose > 2)
                << "~FastQuery: freed index file handler";
#endif
        }
        indexFile = 0;
    }
    if ( dataFile) {
        delete dataFile;
#ifdef DEBUG
        LOGGER(ibis::gVerbose > 2)
            << "~FastQuery: freed data file handler";
#endif
    }
    dataFile = 0;
    delete(metadataMgr);
    metadataMgr = 0;
} // FastQuery::~FastQuery

unsigned int FastQuery::getAllVariables
(std::vector<std::string> &variables,
 const std::string &varPathStr, const std::string &varNameStr)
{
    if (! isValid("FastQuery::getAllVariables")) return 0;

    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    for(unsigned int i=0; i<varInfoList.size(); i++) {
        variables.push_back(varInfoList[i].getPath());
    }
    LOGGER(ibis::gVerbose > 6)
        << "FastQuery::getAllVariable("
        << varPathStr.c_str() << ", " << varNameStr.c_str()
        << ") found " << numVar << " variable" << (numVar > 1 ? "s" : "")
        << " matching the specified name";
    return numVar;
} // FastQuery::getAllVariables

bool FastQuery::getVariableInfo
(const std::string &varNameStr, std::string &variable,
 std::vector <uint64_t> &dims, FQ::DataType *type,
 const std::string &varPathStr)
{
    if (! isValid("FastQuery::getVariableInfo")) return false;

    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList,varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getVariableInfo("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") failed to find a named variable";
        return false;
    }
    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList[i].getPath().size()) {
                idx = i;
                len = varInfoList[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- FastQuery::getVariableInfo("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") found multiple variables with the given name, use \""
            << varInfoList[idx].getPath() << "\"";
    }
    *type = varInfoList[idx].getType();
    variable = varInfoList[idx].getPath().c_str();
    const std::vector<uint64_t> &cnts(varSpaceList[idx].getCounts());
    dims.resize(cnts.size());
    std::copy(cnts.begin(), cnts.end(), dims.begin());
    LOGGER(ibis::gVerbose > 6)
        << "FastQuery::getVariableInfo("
        << varPathStr.c_str() << ", " << varNameStr.c_str()
        << ") successfully got variable information (" << variable
        << ", " << dims.size() << ')';
    return true;
} // FastQuery::getVariableInfo

bool FastQuery::checkForVariable
(const std::string &varNameStr, const std::string &varPathStr)
{
    if (! isValid("FastQuery::checkForVariable")) return false;

    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    return (numVar == 0);
} // FastQuery::checkForVariable

bool FastQuery::getData(const std::string &varNameStr, void *data,
                        const std::string &varPathStr,
                        const bool collective, const bool bcast)
{
    if (! isValid("FastQuery::getData")) return false;

    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") failed to find the named variable";
        return false;
    }
    int idx=0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList[i].getPath().size()) {
                idx = i;
                len = varInfoList[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- FastQuery::getData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") found multiple variables with the name, use \""
            << varInfoList[idx].getPath() << "\"";
    }

    VarInfo varInfo = varInfoList[idx];
    VarSpace varSpace = varSpaceList[idx];
    bool berr;
#ifndef FQ_NOMPI
    if (collective == false) {
#endif
    if (varInfo.getSize() == varSpace.getSize()) {
        berr = dataFile->getData(varInfo.getPath(), data);
    }
    else {
        berr = dataFile->getArrayData(varInfo.getPath(),
                                      varSpace.getOffsets(),
                                      varSpace.getCounts(),
                                      varSpace.getStrides(), data);
    }
#ifndef FQ_NOMPI
    } else {
    int mpi_dim = 0; // split data by the first dimension
    std::vector<uint64_t> offsets = varSpace.getOffsets();
    std::vector<uint64_t> counts = varSpace.getCounts();
    std::vector<uint64_t> strides = varSpace.getStrides();
    int nElements = varSpace.getSize();
    unsigned int totalCount = counts[mpi_dim];
    if (totalCount < mpi_size) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- IndexBuilder::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << "): number of MPI tasks cannot be less than"
            << " the dimension length " << totalCount;
        return false;
    }
    int mpi_len = std::ceil((double)totalCount/(double)mpi_size);
    if (mpi_rank >= totalCount%mpi_len && totalCount%mpi_len != 0) {
        counts[mpi_dim] = mpi_len - 1;
        offsets[mpi_dim] += (mpi_len*mpi_rank -
                             (mpi_rank-totalCount%mpi_len)) *
            strides[mpi_dim];
    }
    else {
        counts[mpi_dim] = mpi_len;
        offsets[mpi_dim] += (mpi_len*mpi_rank)*strides[mpi_dim];
    }
    bool mpi_berr = false;
    FQ::DataType fqType = varInfo.getType();
    void* mpi_data;
    int data_len = counts[mpi_dim]*(nElements/totalCount);
    switch(fqType){
    case FQ::FQT_FLOAT: {
        mpi_data = new float[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_DOUBLE: {
        mpi_data = new double[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_BYTE: {
        mpi_data = new signed char[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_UBYTE: {
        mpi_data = new unsigned char[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_SHORT: {
        mpi_data = new int16_t[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_USHORT: {
        mpi_data = new uint16_t[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_INT: {
        mpi_data = new int32_t[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_UINT: {
        mpi_data = new uint32_t[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_LONG: {
        mpi_data = new int64_t[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    case FQ::FQT_ULONG: {
        mpi_data = new uint64_t[data_len];
        mpi_berr = dataFile->getArrayData(varInfo.getPath(), offsets,
                                          counts, strides, mpi_data);
        break;}
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") can not handle FQ data type " << fqType;
        mpi_berr = false;
    }

    MPI_Allreduce(&mpi_berr, &berr, 1, MPI_BYTE, MPI_BAND, mpi_comm);

    if (berr) {
        int recvcounts[mpi_size];
        int displs[mpi_size];
        recvcounts[0] = mpi_len*(nElements/totalCount);
        displs[0] = 0;
        for(int i=1; i<mpi_size; i++) {
            if (i >= totalCount%mpi_len && totalCount%mpi_len != 0)
                recvcounts[i] = (mpi_len - 1)*(nElements/totalCount);
            else
                recvcounts[i] = mpi_len*(nElements/totalCount);
            displs[i] = displs[i-1] + recvcounts[i-1];
        }

        switch(fqType){
        case FQ::FQT_FLOAT: {
            MPI_Allgatherv(mpi_data, data_len, MPI_FLOAT, data, recvcounts,
                           displs, MPI_FLOAT, mpi_comm);
            break;}
        case FQ::FQT_DOUBLE: {
            MPI_Allgatherv(mpi_data, data_len, MPI_DOUBLE, data,
                           recvcounts, displs, MPI_DOUBLE, mpi_comm);
            break;}
        case FQ::FQT_BYTE: {
            MPI_Allgatherv(mpi_data, data_len, MPI_CHAR, data, recvcounts,
                           displs, MPI_CHAR, mpi_comm);
            break;}
        case FQ::FQT_UBYTE: {
            MPI_Allgatherv(mpi_data, data_len, MPI_UNSIGNED_CHAR, data,
                           recvcounts, displs, MPI_UNSIGNED_CHAR, mpi_comm);
            break;}
        case FQ::FQT_SHORT: {
            MPI_Allgatherv(mpi_data, data_len, MPI_SHORT, data,
                           recvcounts, displs, MPI_SHORT, mpi_comm);
            break;}
        case FQ::FQT_USHORT: {
            MPI_Allgatherv(mpi_data, data_len, MPI_UNSIGNED_SHORT,
                           data, recvcounts, displs,
                           MPI_UNSIGNED_SHORT, mpi_comm);
            break;}
        case FQ::FQT_INT: {
            MPI_Allgatherv(mpi_data, data_len, MPI_INT, data, recvcounts,
                           displs, MPI_INT, mpi_comm);
            break;}
        case FQ::FQT_UINT: {
            MPI_Allgatherv(mpi_data, data_len, MPI_UNSIGNED, data,
                           recvcounts, displs, MPI_UNSIGNED, mpi_comm);
            break;}
        case FQ::FQT_LONG: {
            MPI_Allgatherv(mpi_data, data_len, MPI_LONG_LONG, data,
                           recvcounts, displs, MPI_LONG_LONG, mpi_comm);
            break;}
        case FQ::FQT_ULONG: {
            MPI_Allgatherv(mpi_data, data_len, MPI_UNSIGNED_LONG_LONG,
                           data, recvcounts, displs,
                           MPI_UNSIGNED_LONG_LONG, mpi_comm);
            break;}
        default:
            berr = false;
            break;
        }
    }
    }
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") failed to get data from dataset \""
            << varInfo.getPath().c_str() << "\"";
        return false;
    }
    else {
        LOGGER(ibis::gVerbose > 6)
            << "FastQuery::getData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ") successfully retrieved data from dataset \""
            << varInfo.getPath().c_str() << "\"";
        return true;
    }
} // FastQuery::getData

bool FastQuery::getAttribute(const std::string &varNameStr,
                             const std::string &attrNameStr,
                             void *values,
                             const std::string &varPathStr)
{
    if (! isValid("FastQuery::getAttribute")) return false;
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getAttribute("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") found no variable with the given name";
        return false;
    }
    int idx=0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList[i].getPath().size()) {
                idx = i;
                len = varInfoList[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- FastQuery::getAttribute("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") found multiple variables with the name, use \""
            << varInfoList[idx].getPath() << "\"";
    }
    bool berr = dataFile->getAttribute
        (varInfoList[idx].getPath(), attrNameStr, values);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getAttribute("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") failed to get attribute from the dataset";
        return false;
    }
    else {
        LOGGER(ibis::gVerbose > 2)
            << "FastQuery::getAttribute("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") successfully got attribute value";
        return true;
    }
} // FastQuery::getAttribute

bool FastQuery::getAttributeInfo(const std::string &varNameStr,
                                 const std::string &attrNameStr,
                                 uint64_t *length,
                                 FQ::DataType *type,
                                 const std::string &varPathStr)
{
    if (! isValid("FastQuery::getAttribute")) return false;
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getAttributeInfo("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") failed to find a variable with the given name";
        return false;
    }
    int idx=0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList[i].getPath().size()) {
                idx = i;
                len = varInfoList[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getAttributeInfo("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << " found multiple variables with the name, use \""
            << varInfoList[idx].getPath() << "\"";
    }
    if (! dataFile->getAttributeInfo(varInfoList[idx].getPath(),
                                     attrNameStr, length, type)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FastQuery::getAttributeInfo("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") failed to get attribute from the dataset";
        return false;
    }
    else {
        LOGGER(ibis::gVerbose > 6)
            << "FastQuery::getAttributeInfo("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << ", " << attrNameStr.c_str()
            << ") successfully got attribute information";
        return true;
    }
} // FastQuery::getAttributeInfo

bool FastQuery::isValid() {
    return isValid("FastQuery");
} // FastQuery::isValid

bool FastQuery::isValid(const std::string &func) const {
    if (!dataFile) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << func
            << ": data file did not open successfully";
        return false;
    }
    if (!indexFile) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << func
            << ": index file did not open successfully";
        return false;
    }
    if (!metadataMgr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << func
            << ": info manager is not initialized yet";
        return false;
    }
    return true;
} // FastQuery::isValid

static int compare_double(const void* lhs, const void* rhs)
{
    double* l = (double*) lhs;
    double* r = (double*) rhs;
    return l > r ? 1 : l < r ? -1 : 0;
}

double
util::compute_median(double* arr, int num, bool allow_modify)
{
    if (num <= 0) return 0;
    if (num == 1) return arr[0];

    double* a = arr;
    if (! allow_modify) {
        a = new double[num];
        for (int i=0; i<num; i++) {
            a[i] = arr[i];
        }
    }
    qsort(a, num, sizeof(double), compare_double);

    double ret=0;
    if ((num%2)==0) {
        ret = (a[num/2-1]+a[num/2])/2;
    }
    else {
        ret = a[num/2];
    }
    if (a != arr) delete[] a;
    return ret;
}


