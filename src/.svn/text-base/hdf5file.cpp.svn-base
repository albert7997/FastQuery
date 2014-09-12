#include "hdf5file.h"
#include "fq.h"         // FastQuery::reportTiming()
#include <stack>

/*******************************************
 * Constructor & De-Constructor
 ********************************************/
#ifdef FQ_NOMPI
HDF5::HDF5(const std::string fileName, const bool readOnly,
           const std::string indexPath)
#else
HDF5::HDF5(const std::string fileName, MPI_Comm comm,
           const bool readOnly, const std::string indexPath)
#endif
{
    _fileName = fileName;
    _indexPath = indexPath;
    _fileId = -1;

#ifndef FQ_NOMPI
    mpi_comm = comm;
    MPI_Group mpi_group;
    MPI_Comm_group(comm, &mpi_group);
    MPI_Group_size(mpi_group, &mpi_size);
    MPI_Group_rank(mpi_group, &mpi_rank);
#endif

    if(! __openFile(readOnly)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5: failed to open file \""
            << _fileName.c_str() << "\"";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5: successfully opened the file \""
            << _fileName.c_str() << "\"";
    }
}

HDF5::~HDF5()
{
    if(! __closeFile()) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::~HDF5: failed to close the file \""
            << _fileName.c_str() << "\"";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::~HDF5: successfully closed the file \""
            << _fileName.c_str() << "\"";
    }
}

/*******************************************
 * Public Functions
 ********************************************/
bool HDF5::getAttribute(const std::string &variable,
                        const std::string &attrName, void *values)
{
    return __getAttribute(variable, attrName, values);
} // HDF5::getAttribute

bool HDF5::setAttribute(const std::string &variable,
                        const std::string &attrName,
                        const void *values, const uint64_t len,
                        const FQ::DataType fqType)
{
    return __setAttribute(variable, attrName, values, (hsize_t)len, fqType);
} // HDF5::setAttribute

bool HDF5::getAttributeInfo(const std::string &variable,
                            const std::string &attrName,
                            uint64_t *length, FQ::DataType *fqType)
{
    if ( ! __getAttributeInfo(variable, attrName, length, fqType)) {
        return false;
    }
    switch(*fqType){
    case FQ::FQT_FLOAT:
        *length /= sizeof(float);
        break;
    case FQ::FQT_DOUBLE:
        *length /= sizeof(double);
        break;
    case FQ::FQT_BYTE:
    case FQ::FQT_UBYTE:
        *length /= sizeof(char);
        break;
    case FQ::FQT_SHORT:
    case FQ::FQT_USHORT:
        *length /= sizeof(int16_t);
        break;
    case FQ::FQT_INT:
    case FQ::FQT_UINT:
        *length /= sizeof(int32_t);
        break;
    case FQ::FQT_LONG:
    case FQ::FQT_ULONG:
        *length /= sizeof(int64_t);
        break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getAttributeInfo("
            << variable.c_str() << "):"
            << " unknown FQ data type " << fqType;
        *length =0;
        return false;
    }
    return true;
} // HDF5::getAttributeInfo

bool HDF5::createDataset(const std::string& datasetName,
                         const std::vector<uint64_t > dims,
                         const FQ::DataType fqType)
{
    if (dims.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createDataset(" << datasetName
            << "): the number of dataset dimensions is 0";
        return false;
    }
    if (! __createDataset(datasetName, dims, fqType, false)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createDataset(" << datasetName
            << ") failed to create the dataset";
        return false;
    }
    return true;
} // HDF5::createDataset

bool HDF5::setData(const std::string &variable, const void *data)
{
    if (! __writeData(variable, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setData(" << variable
            << ") failed to set data";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::setData(" << variable
            << ") successfully set data";
        return true;
    }
} // HDF5::setData

bool HDF5::getData(const std::string &variable, void *data)
{
    if (! __readData(variable, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getData:"
            << " failed to get data from dataset \""
            << variable.c_str() << "\"";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::getData:"
            << " successfully got data from dataset \""
            << variable.c_str() << "\"";
        return true;
    }
} // HDF5::getData

bool HDF5::setArrayData(const std::string &variable,
                        const std::vector<uint64_t> &offsets,
                        const std::vector<uint64_t> &counts,
                        const std::vector<uint64_t> &strides,
                        const void *data)
{
    if (! __writeArrayData(variable, offsets, counts, strides, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setArrayData("
            << variable.c_str() << "):"
            << " failed to set data";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 3)
            << "HDF5::setArrayData("
            << variable.c_str() << "):"
            << " successfully set data";
        return true;
    }
} // HDF5::setDataSection

bool HDF5::getArrayData(const std::string &variable,
                        const std::vector<uint64_t> &offsets,
                        const std::vector<uint64_t> &counts,
                        const std::vector<uint64_t> &strides, void *data)
{
    if (! __readArrayData(variable, offsets, counts, strides, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getArrayData("
            << variable.c_str() << "):"
            << " failed to get data";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 3)
            << "HDF5::getArrayData("
            << variable.c_str() << "):"
            << " successfully got data";
        return true;
    }
} // HDF5::getDataSection

bool HDF5::getPointData(const std::string &variable,
                        const std::vector<uint64_t> &coords, void *data)
{
    if (! __readPointData(variable, coords, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getPointData("
            << variable.c_str() << "):"
            << " failed to get data";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::getPointData("
            << variable.c_str() << "):"
            << " successfully got data";
        return true;
    }
} // HDF5::getPointData

bool HDF5::getBitmapKeys(const std::string &variable, void *keys,
                         const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    berr =  __readData(datasetName, keys);
#else
    uint64_t start, end;
    if (! __getOffsets(variable, BitmapKeyColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getBitmapKeys("
            << variable.c_str() << "):"
            << " failed to get bitmap keys length";
        return false;
    }

    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapKeys";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.resize(1);
    counts.resize(1);
    strides.resize(1);
    offsets[0] = start;
    counts[0] = end-start;
    strides[0] = 1;
    berr = __readArrayData(datasetName, offsets, counts, strides, keys);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getBitmapKeys("
            << variable.c_str() << "):"
            << " failed to get bitmap keys";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getBitmapKeys("
        << variable.c_str() << "):"
        << " successfully got bitmap keys";
    return true;
} // HDF5::getBitmapKeys

bool HDF5::setBitmapKeys(const std::string &variable, const void *keys,
                         const uint64_t nkeys, const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    berr =  __writeData(datasetName, keys);
#else
    uint64_t start, end;
    if (! __getOffsets(variable, BitmapKeyColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapKeys("
            << variable.c_str() << "):"
            << " failed to get bitmap keys length";
        return false;
    }

    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapKeys";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.resize(1);
    counts.resize(1);
    strides.resize(1);
    offsets[0] = start;
    counts[0] = end-start;
    strides[0] = 1;
    berr = __writeArrayData(datasetName, offsets, counts, strides, keys,
                            HDF5_COLLECTIVE_WRITE);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapKeys("
            << variable.c_str() << "):"
            << " failed to set bitmap keys";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::setBitmapKeys("
        << variable.c_str() << "):"
        << " successfully set bitmap keys";
    return true;
} // HDF5::setBitmapKeys

// must be called collectively in MPI mode
bool HDF5::createBitmapKeys(const std::string &variable, const uint64_t nkeys,
                            const FQ::DataType fqType, const uint64_t mpi_iter,
                            const uint64_t mpi_idx)
{
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";

    std::vector<uint64_t> dims;
    dims.resize(1);
    dims[0] = nkeys;
    if (! __createDataset(datasetName, dims, fqType, false)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmapKeys("
            << variable.c_str() << "):"
            << " failed to create bitmap keys dataset";
        return false;
    }
#else
    // get bitmapKey length
    uint64_t curLen = nkeys;
    if (mpi_rank == 0 && mpi_iter > 0) {
        uint64_t start, end;
        __getOffsets(variable, BitmapKeyColIdx, &start, &end,
                     mpi_iter*mpi_size-1);
        curLen = end;
    }
    MPI_Bcast(&curLen, 1, MPI_UNSIGNED_LONG, 0, mpi_comm);
    if (mpi_iter > 0 && curLen == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmapKeys("
            << variable.c_str() << "):"
            << " failed to get bitmap keys length";
        return false;
    }
    uint64_t extraLen = 0;
    uint64_t val = nkeys;
    MPI_Allreduce(&val, &extraLen, 1, MPI_UNSIGNED_LONG, MPI_SUM, mpi_comm);

    // create or extend the bitmapKey dataset collectively
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapKeys";
    std::vector<uint64_t> dims;
    dims.push_back(curLen);
    uint64_t totLen;
    // first create a dataset with the chunk size = nkeys from the first
    // processor
    if (mpi_iter == 0) {
        totLen = extraLen;
        if (! __createDataset(datasetName, dims, fqType, true)) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- HDF5::setBitmapKeys("
                << variable.c_str() << "):"
                << " failed to create the bitmap keys dataset";
            return false;
        }
    } else {
        totLen = curLen+extraLen;
    }
    // then extend the dataset to the total length
    if (! __extendDataset(datasetName, dims, totLen)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapKeys("
            << variable.c_str() << "):"
            << " failed to extend the bitmap keys dataset";
        return false;
    }
    if (FastQuery::reportTiming()) {
        LOGGER(true) << "Statistic\tBitmapKeySize\t" << extraLen;
    }
#endif
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::createBitmapKeys("
        << variable.c_str() << "):"
        << " successfully create bitmap keys dataset";
    return true;
} // HDF5::createBitmapKeys

bool HDF5::getBitmapKeyLength(const std::string &variable, uint64_t *nkeys,
                              const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    std::vector<uint64_t> dims;
    berr = __getDatasetDimension(datasetName, dims);
    if (berr) *nkeys = dims[0];
#else
    berr = __getOffsetLength(variable, BitmapKeyColIdx, nkeys, mpi_idx);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- HDF5::getBitmapKeyLength("
            << variable.c_str() << "):"
            << " failed to get bitmap keys length";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getBitmapKeyLength("
        << variable.c_str() << "):"
        << " successfully got bitmap keys length " << (uint64_t)(*nkeys);
    return true;
} // HDF5::getBitmapKeysLength

#ifndef FQ_NOMPI
//must be called collectively in MPI mode
bool HDF5::setBitmapKeyLength(const std::string &variable,
                              const uint64_t nkeys, const uint64_t mpi_iter)
{
    bool berr = true;
    uint64_t* nkeysArray;
    if (mpi_rank == 0) nkeysArray = new uint64_t[mpi_size];
    uint64_t val = nkeys;
    MPI_Gather(&val, 1, MPI_UNSIGNED_LONG, nkeysArray, 1,
               MPI_UNSIGNED_LONG, 0, mpi_comm);

    if (mpi_rank == 0) {
        // first mpi task sets the bitmapKey length for all other tasks
        uint64_t curLen = 0;
        if (mpi_iter > 0) {
            uint64_t start, end;
            berr = __getOffsets(variable, BitmapKeyColIdx, &start, &end,
                                mpi_iter*mpi_size-1);
            if (! berr) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::setBitmapKeyLength("
                    << variable.c_str() << "):"
                    << " failed to get bitmap keys length";
            } else {
                curLen = end;
            }
        }
        if (berr) {
            nkeysArray[0] += curLen;
            for(unsigned int i=1; i<mpi_size; i++) {
                nkeysArray[i] += nkeysArray[i-1];
            }
            std::vector<uint64_t> offsets;
            std::vector<uint64_t> counts;
            std::vector<uint64_t> strides;
            offsets.resize(2);
            counts.resize(2);
            strides.resize(2);
            offsets[0] = mpi_iter*mpi_size+1;
            counts[0] = mpi_size;
            strides[0] = 1;
            offsets[1] = BitmapKeyColIdx;
            counts[1] = 1;
            strides[1] = 1;
            std::string datasetName = _indexPath;
            datasetName += variable;
            datasetName += ".MPIoffsetTable";
            berr = __writeArrayData(datasetName, offsets, counts, strides,
                                    nkeysArray);
        }
    }
    if (mpi_rank == 0) delete[] nkeysArray;
    MPI_Bcast(&berr, 1, MPI_INT, 0, mpi_comm);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapKeyLength("
            << variable.c_str() << "):"
            << " failed to set bitmap keys length";
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::setBitmapKeyLength("
        << variable.c_str() << "):"
        << " successfully set bitmap keys length";
    return true;
}
#endif

bool HDF5::getBitmapOffsets(const std::string &variable, void *bitmapOffsets,
                            const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";
    berr =  __readData(datasetName, bitmapOffsets);
#else
    uint64_t start, end;
    if (! __getOffsets(variable, BitmapOffsetColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to get bitmap offsets length";
        return false;
    }

    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapOffsets";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.resize(1);
    counts.resize(1);
    strides.resize(1);
    offsets[0] = start;
    counts[0] = end-start;
    strides[0] = 1;
    berr = __readArrayData(datasetName, offsets, counts, strides,
                           bitmapOffsets);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to get bitmap offsets";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getBitmapOffsets("
        << variable.c_str() << "):"
        << " successfully got bitmap offsets";
    return true;
} // HDF5::getBitmapOffsets

bool HDF5::setBitmapOffsets(const std::string &variable,
                            const void *bitmapOffsets,
                            const uint64_t noffsets, const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";
    berr =  __writeData(datasetName, bitmapOffsets);
#else
    uint64_t start, end;
    if (! __getOffsets(variable, BitmapOffsetColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to get bitmap offsets length";
        return false;
    }

    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapOffsets";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.resize(1);
    counts.resize(1);
    strides.resize(1);
    offsets[0] = start;
    counts[0] = end-start;
    strides[0] = 1;
    berr = __writeArrayData(datasetName, offsets, counts, strides,
                            bitmapOffsets, HDF5_COLLECTIVE_WRITE);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to set bitmap offsets";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::setBitmapOffsets("
        << variable.c_str() << "):"
        << " successfully set bitmap offsets with size " << noffsets;
    return true;
} // HDF5::setBitmapOffsets

// must be called collectively in MPI mode
bool HDF5::createBitmapOffsets(const std::string &variable,
                               const uint64_t noffsets,
                               const FQ::DataType fqType,
                               const uint64_t mpi_iter,
                               const uint64_t mpi_idx)
{
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";

    std::vector<uint64_t> dims;
    dims.resize(1);
    dims[0] = noffsets;
    if (! __createDataset(datasetName, dims, fqType, false)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to create bitmap offsets dataset";
        return false;
    }
#else
    // get bitmapOffset length
    uint64_t curLen = noffsets;
    if (mpi_rank == 0 && mpi_iter > 0) {
        uint64_t start, end;
        __getOffsets(variable, BitmapOffsetColIdx, &start, &end,
                     mpi_iter*mpi_size-1);
        curLen = end;
    }
    MPI_Bcast(&curLen, 1, MPI_UNSIGNED_LONG, 0, mpi_comm);
    if (mpi_iter > 0 && curLen == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to get bitmap offsets length";
        return false;
    }
    uint64_t extraLen = 0;
    uint64_t val = noffsets;
    MPI_Allreduce(&val, &extraLen, 1, MPI_UNSIGNED_LONG, MPI_SUM, mpi_comm);

    // create or extend the bitmapOffset dataset collectively
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapOffsets";
    std::vector<uint64_t> dims;
    dims.push_back(curLen);
    uint64_t totLen;
    if (mpi_iter == 0) {
        totLen = extraLen;
        if (! __createDataset(datasetName, dims, fqType, true)) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- HDF5::createBitmapOffsets("
                << variable.c_str() << "):"
                << " failed to create the bitmap offsets dataset";
            return false;
        }
    } else {
        totLen = curLen+extraLen;
    }
    if (! __extendDataset(datasetName, dims, totLen)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to extend the bitmap offsets dataset";
        return false;
    }
    if (FastQuery::reportTiming()) {
        LOGGER(true) << "Statistic\tBitmapOffsetSize\t" << extraLen;
    }
#endif
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::createBitmapOffsets("
        << variable.c_str() << "):"
        << " successfully created the bitmap offsets dataset";
    return true;
} // HDF5::createBitmapOffsets

bool HDF5::getBitmapOffsetLength(const std::string &variable,
                                 uint64_t *noffsets, const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";
    //berr = getDataLength(datasetName, noffsets);
    std::vector<uint64_t> dims;
    berr = __getDatasetDimension(datasetName, dims);
    if (berr) *noffsets = dims[0];
#else
    berr = __getOffsetLength(variable, BitmapOffsetColIdx, noffsets, mpi_idx);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- HDF5::getBitmapOffsetsLength("
            << variable.c_str() << "):"
            << " failed to get bitmap offsets length";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getBitmapOffsetsLength("
        << variable.c_str() << "):"
        << " successfully got bitmap offsets length "
        << (uint64_t)(*noffsets);
    return true;
} // HDF5::getBitmapOffsetsLength

#ifndef FQ_NOMPI
//must be called collectively in MPI mode
bool HDF5::setBitmapOffsetLength(const std::string &variable,
                                 const uint64_t noffsets,
                                 const uint64_t mpi_iter)
{
    bool berr = true;
    uint64_t* noffsetsArray;
    if (mpi_rank == 0) noffsetsArray = new uint64_t[mpi_size];
    uint64_t val = noffsets;
    MPI_Gather(&val, 1, MPI_UNSIGNED_LONG, noffsetsArray, 1,
               MPI_UNSIGNED_LONG, 0, mpi_comm);

    if (mpi_rank == 0) {
        // first mpi task sets the bitmapOffset length for all other tasks
        uint64_t curLen = 0;
        if (mpi_iter > 0) {
            uint64_t start, end;
            berr = __getOffsets(variable, BitmapOffsetColIdx, &start, &end,
                                mpi_iter*mpi_size-1);
            if (! berr){
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::setBitmapOffsetLength("
                    << variable.c_str() << "):"
                    << " failed to get bitmap offsets length";
            } else {
                curLen = end;
            }
        }
        if (berr) {
            noffsetsArray[0] += curLen;
            for(unsigned int i=1; i<mpi_size; i++) {
                noffsetsArray[i] += noffsetsArray[i-1];
            }
            std::vector<uint64_t> offsets;
            std::vector<uint64_t> counts;
            std::vector<uint64_t> strides;
            offsets.resize(2);
            counts.resize(2);
            strides.resize(2);
            offsets[0] = mpi_iter*mpi_size+1;
            counts[0] = mpi_size;
            strides[0] = 1;
            offsets[1] = BitmapOffsetColIdx;
            counts[1] = 1;
            strides[1] = 1;
            std::string datasetName = _indexPath;
            datasetName += variable;
            datasetName += ".MPIoffsetTable";
            berr = __writeArrayData(datasetName, offsets, counts,
                                    strides, noffsetsArray);
        }
    }
    if (mpi_rank == 0) delete[] noffsetsArray;
    MPI_Bcast(&berr, 1, MPI_INT, 0, mpi_comm);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapOffsetLength("
            << variable.c_str() << "):"
            << " failed to set bitmap offsets length";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::setBitmapOffsetLength("
        << variable.c_str() << "):"
        << " successfully set bitmap offsets length";
    return true;
}
#endif

FQ::DataType HDF5::getBitmapOffsetType(const std::string& variable)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
#ifdef FQ_NOMPI
    datasetName += ".bitmapOffsets";
#else
    datasetName += ".MPIbitmapOffsets";
#endif

    FQ::DataType fqType;
    if (!  __getDatasetType(datasetName, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getBitmapOffsetType("
            << variable.c_str() << "):"
            << " failed to get bitmap offset type";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::getBitmapOffsetType("
            << variable.c_str() << "):"
            << " successfully got bitmap offset FQ type " << fqType;
    }
    return fqType;
} // HDF5::getBitmapOffsetType

bool HDF5::readBitmap(const std::string& variable, const uint64_t startoffset,
                      const uint64_t endoffset, uint32_t *data,
                      const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
#ifdef FQ_NOMPI
    datasetName += ".bitmap";
#else
    datasetName += ".MPIbitmap";
#endif

    uint64_t count = endoffset - startoffset;
    if (count == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::readBitmap("
            << variable.c_str() << "):"
            << " offset size is 0";
        return false;
    }

    std::vector<uint64_t> offsets;
    offsets.resize(1);
    offsets[0] = startoffset;

    std::vector<uint64_t> counts;
    counts.resize(1);
    counts[0] = count;

    std::vector<uint64_t> strides;
    strides.resize(1);
    strides[0] = 1;

#ifndef FQ_NOMPI
    uint64_t start, end;
    if (! __getOffsets(variable, BitmapColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::readBitmap("
            << variable.c_str() << "):"
            << " failed to get bitmap length";
        return false;
    }
    offsets[0] += start;
#endif

    if (! __readArrayData(datasetName, offsets, counts, strides, (void*)data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::readBitmap("
            << variable.c_str() << "):"
            << " failed to read bitmap with size " << count
            << "[" << startoffset << " - " << endoffset << "]";
        return false;
    }
    LOGGER(ibis::gVerbose > 3)
        << "HDF5::readBitmap("
        << variable.c_str() << "):"
        << " successfully read bitmap with size " << count
        << "[" << startoffset << " - " << endoffset << "]";
    return true;
} // HDF5::readBitmap

bool HDF5::writeBitmap(const std::string& variable,
                       const uint64_t startoffset,
                       const uint64_t endoffset,
                       const uint32_t *data, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
#ifdef FQ_NOMPI
    datasetName += ".bitmap";
#else
    datasetName += ".MPIbitmap";
#endif

    uint64_t count = endoffset - startoffset;
    if (count == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::writeBitmap("
            << variable.c_str() << "):"
            << " offset size is 0";
        return false;
    }

    std::vector<uint64_t> offsets;
    offsets.resize(1);
    offsets[0] = startoffset;

    std::vector<uint64_t> counts;
    counts.resize(1);
    counts[0] = count;

    std::vector<uint64_t> strides;
    strides.resize(1);
    strides[0] = 1;

#ifndef FQ_NOMPI
    uint64_t start, end;
    if (! __getOffsets(variable, BitmapColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::writeBitmap("
            << variable.c_str() << "):"
            << " failed to get bitmap length";
        return false;
    }
    offsets[0] += start;
#endif

    if (! __writeArrayData(datasetName, offsets, counts, strides,
                           (void*)data, HDF5_COLLECTIVE_WRITE)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::writeBitmap("
            << variable.c_str() << "):"
            << " failed to write bitmap with size " << (endoffset-startoffset)
            << "[" << startoffset << " - " << endoffset << "]";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::writeBitmap("
        << variable.c_str() << "):"
        << " successfully wrote bitmap with size " << (endoffset-startoffset)
        << "[" << startoffset << " - " << endoffset << "]";
    return true;
} // HDF5::writeBitmap

// must be called collectively in MPI mode
bool HDF5::createBitmap(const std::string &variable,
                        const uint64_t nElements,
                        const uint64_t mpi_iter, const uint64_t mpi_idx)
{
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmap";

    std::vector<uint64_t> dims;
    dims.resize(1);
    dims[0] = nElements;
    if (! __createDataset(datasetName, dims, FQ::FQT_UINT, false)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmap("
            << variable.c_str() << "):"
            << " failed to create bitmap dataset";
        return false;
    }
#else
    // get bitmap length
    uint64_t curLen = nElements;
    if (mpi_rank == 0 && mpi_iter > 0) {
        uint64_t start, end;
        __getOffsets(variable, BitmapColIdx, &start, &end,
                     mpi_iter*mpi_size-1);
        curLen = end;
    }
    MPI_Bcast(&curLen, 1, MPI_UNSIGNED_LONG, 0, mpi_comm);
    if (mpi_iter > 0 && curLen == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmap("
            << variable.c_str() << "):"
            << " failed to get bitmap length";
        return false;
    }
    uint64_t extraLen = 0;
    uint64_t val = nElements;
    MPI_Allreduce(&val, &extraLen, 1, MPI_UNSIGNED_LONG, MPI_SUM, mpi_comm);

    // create or extend the bitmap dataset collectively
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmap";
    std::vector<uint64_t> dims;
    dims.push_back(curLen);
    uint64_t totLen;
    if (mpi_iter == 0) {
        totLen = extraLen;
        if (! __createDataset(datasetName, dims, FQ::FQT_UINT, true)) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- HDF5::createBitmap("
                << variable.c_str() << "):"
                << " failed to create bitmap dataset";
            return false;
        }
    } else {
        totLen = curLen+extraLen;
    }
    if (! __extendDataset(datasetName, dims, totLen)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createBitmap("
            << variable.c_str() << "):"
            << " failed to extend the bitmap dataset";
        return false;
    }
    if (FastQuery::reportTiming()) {
        LOGGER(true) << "Statistic\tBitmapSize\t" << extraLen;
    }
#endif
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::createBitmap("
        << variable.c_str() << "):"
        << " successfully created bitmap dataset";
    return true;
} // HDF5::createBitmap

bool HDF5::getBitmapLength(const std::string &variable,
                           uint64_t *nElements, const uint64_t mpi_idx)
{
    bool berr = true;
#ifdef FQ_NOMPI
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmap";
    //berr = getDataLength(datasetName, nElements);
    std::vector<uint64_t> dims;
    berr = __getDatasetDimension(datasetName, dims);
    if (berr) *nElements = dims[0];
#else
    berr = __getOffsetLength(variable, BitmapColIdx, nElements, mpi_idx);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- HDF5::getBitmapLength("
            << variable.c_str() << "):"
            << " failed to get bitmap length";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getBitmapLength("
        << variable.c_str() << "):"
        << " successfully got bitmap length " << *nElements;
    return true;
} // HDF5::getBitmapLength

#ifndef FQ_NOMPI
//must be called collectively in MPI mode
bool HDF5::setBitmapLength(const std::string &variable,
                           const uint64_t nElements, const uint64_t mpi_iter)
{
    bool berr = true;
    uint64_t* nElementsArray;
    if (mpi_rank == 0) nElementsArray = new uint64_t[mpi_size];
    uint64_t val = nElements;
    MPI_Gather(&val, 1, MPI_UNSIGNED_LONG, nElementsArray, 1,
               MPI_UNSIGNED_LONG, 0, mpi_comm);

    if (mpi_rank == 0) {
        // first mpi task sets the bitmap length for all other tasks
        uint64_t curLen = 0;
        if (mpi_iter > 0) {
            uint64_t start, end;
            berr = __getOffsets(variable, BitmapColIdx, &start, &end,
                                mpi_iter*mpi_size-1);
            if (! berr){
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::setBitmapLength("
                    << variable.c_str() << "):"
                    << " failed to set bitmap length";
            } else {
                curLen = end;
            }
        }
        if (berr) {
            nElementsArray[0] += curLen;
            for(unsigned int i=1; i<mpi_size; i++) {
                nElementsArray[i] += nElementsArray[i-1];
            }
            std::vector<uint64_t> offsets;
            std::vector<uint64_t> counts;
            std::vector<uint64_t> strides;
            offsets.resize(2);
            counts.resize(2);
            strides.resize(2);
            offsets[0] = mpi_iter*mpi_size+1;
            counts[0] = mpi_size;
            strides[0] = 1;
            offsets[1] = BitmapColIdx;
            counts[1] = 1;
            strides[1] = 1;
            std::string datasetName = _indexPath;
            datasetName += variable;
            datasetName += ".MPIoffsetTable";
            berr = __writeArrayData(datasetName, offsets, counts, strides,
                                    nElementsArray);
        }
    }
    if (mpi_rank == 0) delete[] nElementsArray;
    MPI_Bcast(&berr, 1, MPI_INT, 0, mpi_comm);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::setBitmapLength("
            << variable.c_str() << "):"
            << " failed to set bitmap length";
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::setBitmapLength("
        << variable.c_str() << "):"
        << " successfully set bitmap length";
    return true;
}

// must be called collectively in MPIO mode
bool HDF5::createOffsetTable(const std::string &variable,
                             const uint64_t mpi_max_iter)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIoffsetTable";
    uint64_t nElements = mpi_max_iter * (uint64_t)mpi_size + 1;
    std::vector<uint64_t> dims;
    dims.push_back(nElements); // the first element is always 0
    dims.push_back(3); // 3 columns to store the offset of bitmap, bitmapoffsets and bitmapkeys

    if (! __createDataset(datasetName, dims, FQ::FQT_ULONG, false)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::createOffsetTable("
            << variable.c_str() << "):"
            << " failed to create offset table dataset";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::createOffsetTable("
        << variable.c_str() << "):"
        << " successfully created offset table dataset with length "
        << nElements;
    return true;
}
#endif

bool HDF5::getActualRange(const std::string &variable, void *range)
{
#ifndef FQ_NOMPI
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- " << "HDF5::getActualRange("
        << variable.c_str() << "):"
        << " it is not support in MPI mode yet";
    return false;
#else
    std::string attrName = "actualRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getAttribute(datasetName, attrName, range)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::getActualRange("
            << variable.c_str() << "):"
            << " failed to get actual range";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::getActualRange("
            << variable.c_str() << "):"
            << " successfully got actual range";
        return true;
    }
#endif
} // HDF5::getActualRange

bool HDF5::setActualRange(const std::string &variable, const void *range,
                          const FQ::DataType fqType)
{
#ifndef FQ_NOMPI
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- " << "HDF5::setActualRange("
        << variable.c_str() << "):"
        << " it is not support in MPI mode yet";
    return false;
#else
    std::string attrName = "actualRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __setAttribute(datasetName, attrName, range, 2, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::setActualRange("
            << variable.c_str() << "):"
            << " failed to set actual range";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::setActualRange("
            << variable.c_str() << "):"
            << " successfully set actual range";
        return true;
    }
#endif
} // HDF5::setActualRange

bool HDF5::getExpectedRange(const std::string &variable, void *range)
{
#ifndef FQ_NOMPI
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- " << "HDF5::getExpectedRange("
        << variable.c_str() << "):"
        << " it is not support in MPI mode yet";
    return false;
#else
    std::string attrName = "expectedRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getAttribute(datasetName, attrName, range)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::getExpectedRange("
            << variable.c_str() << "):"
            << " failed to get expected range";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::getExpectedRange("
            << variable.c_str() << "):"
            << " successfully got expected range";
        return true;
    }
#endif
} // HDF5::getExpectedRange

bool HDF5::setExpectedRange(const std::string &variable, const void *range,
                            const FQ::DataType fqType)
{
#ifndef FQ_NOMPI
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- " << "HDF5::setExpectedRange("
        << variable.c_str() << "):"
        << " it is not support in MPI mode yet";
    return false;
#else
    std::string attrName = "expectedRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __setAttribute(datasetName, attrName, range, 2, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::setExpectedRange("
            << variable.c_str() << "):"
            << " failed to set expected range";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::setExpectedRange("
            << variable.c_str() << "):"
            << " successfully set expected range";
        return true;
    }
#endif
} // HDF5::setExpectedRange

bool HDF5::getVariableInfo(const std::string &variable,
                           std::vector <uint64_t> &dims, FQ::DataType *fqType)
{
    if(! __getDatasetType(variable, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::getVariableInfo(" << variable
            << "): failed to get variable type";
        return false;
    }
    if(! __getDatasetDimension(variable, dims)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::getVariableInfo(" << variable
            << "): failed to get variable dimension";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getVariableInfo(" << variable
        << "): successfully get variable information dimension "
        << dims.size() << " and FQ type " << *fqType;
    return true;
} // HDF5::getVariableInfo

bool HDF5::getAllVariables(const std::string &path,
                           std::vector<std::string> &variables)
{
    if(! __getAllVariables(path, variables)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::getAllVariables(" << path
            << "): failed to get variables";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::getAllVariables(" << path.c_str()
            << "): successfully found number of variables: "
            << variables.size();
        return true;
    }
}

std::string HDF5::getSortedFieldName()
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::getSortedFieldName"
            << " file is not opened " << _fileName.c_str();
        return "";
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::getSortedFieldName"
        << " it is not implemented yet";
    return "";
} // HDF5::getSortedFieldName

/*******************************************
 * Private Functions
 ********************************************/
bool HDF5::__openFile(const bool readOnly)
{
    ibis::horometer timer;
    timer.start();
    //make sure the file is a valid HDF5 file and exists...
    htri_t  temp;
    H5E_BEGIN_TRY{
        temp = H5Fis_hdf5(_fileName.c_str());
    }H5E_END_TRY;

    hid_t acc_tpl = H5P_DEFAULT;
#ifndef FQ_NOMPI
    acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HDF5_MPI_POSIX
    //H5Pset_fapl_mpiposix(acc_tpl, mpi_comm, MPI_INFO_NULL);
    H5Pset_fapl_mpiposix(acc_tpl, mpi_comm, 0);
#else
    //H5Pset_fapl_mpio(acc_tpl, mpi_comm, MPI_INFO_NULL);
    H5Pset_fapl_mpio(acc_tpl, mpi_comm, 0);
#endif
#ifdef HDF5_DATA_ALIGNMENT
    H5Pset_alignment(acc_tpl, 0, FS_STRIPE_SIZE);
#endif
#endif

#ifdef HDF5_DISABLE_METADATA_FLUSH
    // disable metadata flush for performance
#ifdef FQ_NOMPI
    acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
#endif
    H5AC_cache_config_t mdc_config;
    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(acc_tpl, &mdc_config);
    mdc_config.evictions_enabled = false;
    mdc_config.flash_incr_mode = H5C_flash_incr__off;
    mdc_config.incr_mode = H5C_incr__off;
    mdc_config.decr_mode = H5C_decr__off;
    H5Pset_mdc_config(acc_tpl, &mdc_config);
#endif

    if(temp > 0){
        if ( readOnly ) {
            LOGGER(ibis::gVerbose > 2)
                << "HDF5::__openFile:"
                << " attampt to open file \"" << _fileName.c_str()
                << "\" with read permission";
            _fileId = H5Fopen(_fileName.c_str(), H5F_ACC_RDONLY, acc_tpl);
        } else {
            LOGGER(ibis::gVerbose > 2)
                << "HDF5::__openFile:"
                << " attampt to open file \"" << _fileName.c_str()
                << "\" with write permission";
            _fileId = H5Fopen(_fileName.c_str(), H5F_ACC_RDWR, acc_tpl);
        }
    } else {
        if (readOnly) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "HDF5::__openFile:"
                << " cannot open non-exist or invalid HDF5 file \""
                << _fileName.c_str() << "\" with read permission";
            H5Pclose(acc_tpl);
            return false;
        } else {
            LOGGER(ibis::gVerbose > 2)
                << "HDF5::__openFile:"
                << " attampt to create file \"" << _fileName.c_str()
                << "\" with write permission";
            _fileId = H5Fcreate(_fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                                acc_tpl);
        }
    }
    H5Pclose(acc_tpl);
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__openFile:"
            << " fail to open/create file \""
            << _fileName.c_str() << "\"";
        return false;
    }
    timer.stop();
    if (FastQuery::reportTiming()) {
        LOGGER(true) << "Statistic\tHDF5::openFile\t"
                     << timer.CPUTime() << "\t" << timer.realTime();
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__openFile:"
        << " successfully open file \"" << _fileName.c_str() << "\""
        << " used " << timer.CPUTime() << " sec CPU time and "
        << timer.realTime() << " sec elapsed time";
    return true;
} // HDF5::__openFile

bool HDF5::__closeFile(){
    ibis::horometer timer;
    timer.start();
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 2)
            << "HDF5::__closeFile: "
            << "file \"" << _fileName.c_str() << "\" was not opened";
        return true;
    }
#if (! defined(FQ_NOMPI)) && defined(HDF5_DISABLE_METADATA_FLUSH)
    H5Fflush(_fileId, H5F_SCOPE_LOCAL);
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__closeFile: "
        << "flush the metadata to file \"" << _fileName.c_str() << "\"";
#endif
    if (H5Fclose(_fileId) < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__closeFile: "
            << "failed to close file \"" << _fileName.c_str() << "\"";
        return false;
    }
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tHDF5::closeFile\t"
        << timer.CPUTime() << "\t" << timer.realTime();
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__closeFile:"
        << " successfully closed file \"" << _fileName.c_str() << "\""
        << " used " << timer.CPUTime() << " sec CPU time and "
        << timer.realTime() << " sec elapsed time";
    return true;
} // HDF5::__closeFile

hid_t HDF5::__getHDF5DataType(const FQ::DataType fqType)
{
    hid_t type = -1;
    switch(fqType) {
    case FQ::FQT_FLOAT: {
        type =  H5T_NATIVE_FLOAT;
        break;}
    case FQ::FQT_DOUBLE: {
        type = H5T_NATIVE_DOUBLE;
        break;}
    case FQ::FQT_BYTE: {
        type = H5T_NATIVE_SCHAR;
        break;}
    case FQ::FQT_UBYTE: {
        type = H5T_NATIVE_UCHAR;
        break;}
    case FQ::FQT_SHORT: {
        type = H5T_NATIVE_INT16;
        break;}
    case FQ::FQT_USHORT: {
        type = H5T_NATIVE_UINT16;
        break;}
    case FQ::FQT_INT: {
        type = H5T_NATIVE_INT32;
        break;}
    case FQ::FQT_UINT: {
        type = H5T_NATIVE_UINT32;
        break;}
    case FQ::FQT_LONG: {
        type = H5T_NATIVE_INT64;
        break;}
    case FQ::FQT_ULONG: {
        type = H5T_NATIVE_UINT64;
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getHDF5DataType("
            << fqType << "):"
            << " HDF5 does not yet support the FQ data type";
    }
    }
    return type;
} // HDF5::__getHDF5DataType

/// Convert an atomic HDF5 data type to FQ::DataType.
FQ::DataType HDF5::__getFQDataType(const hid_t type) {
    FQ::DataType fqType = FQ::FQT_UNKNOWN;
    hid_t native = H5Tget_native_type(type, H5T_DIR_ASCEND);
    if (H5Tequal(native, H5T_NATIVE_FLOAT) > 0) {
        fqType = FQ::FQT_FLOAT;
    }
    else if (H5Tequal(native, H5T_NATIVE_DOUBLE) > 0) {
        fqType = FQ::FQT_DOUBLE;
    }
    else if (H5Tequal(native, H5T_NATIVE_SCHAR) > 0) {
        fqType = FQ::FQT_BYTE;
    }
    else if (H5Tequal(native, H5T_NATIVE_UCHAR) > 0) {
        fqType = FQ::FQT_UBYTE;
    }
    else if (H5Tequal(native, H5T_NATIVE_INT16) > 0) {
        fqType = FQ::FQT_SHORT;
    }
    else if (H5Tequal(native, H5T_NATIVE_UINT16) > 0) {
        fqType = FQ::FQT_USHORT;
    }
    else if (H5Tequal(native, H5T_NATIVE_INT32) > 0) {
        fqType = FQ::FQT_INT;
    }
    else if (H5Tequal(native, H5T_NATIVE_UINT32) > 0) {
        fqType = FQ::FQT_UINT;
    }
    else if (H5Tequal(native, H5T_NATIVE_INT64) > 0) {
        fqType = FQ::FQT_LONG;
    }
    else if (H5Tequal(native, H5T_NATIVE_UINT64) > 0) {
        fqType = FQ::FQT_ULONG;
    }
    else {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getFQDataType does not yet support HDF5 type "
            << (unsigned long) type;
    }
    H5Tclose(native);
    return fqType;
} // HDF5::__getFQDataType

H5O_type_t HDF5::__getHDF5ObjType(const std::string &path)
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getHDF5ObjType("
            << path.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return H5O_TYPE_UNKNOWN;
    }
    if (path.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getHDF5ObjType("
            << path.c_str() << "):"
            << " path is an empty string";
        return H5O_TYPE_UNKNOWN;
    }
    if (path[0] != '/') {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getHDF5ObjType("
            << path.c_str() << "):"
            << " path must be given as an absolute path starting with '/'";
        return H5O_TYPE_UNKNOWN;
    }
    size_t pos = 0;
    while(1) {
        pos = path.find('/', pos+1);
        if (pos == path.npos) break;
        std::string tmp = path.substr(0,pos);
        if (! H5Lexists(_fileId, tmp.c_str(), H5P_DEFAULT)){
            return H5O_TYPE_UNKNOWN;
        }
    }
    if (! H5Lexists(_fileId, path.c_str(), H5P_DEFAULT)) {
        return H5O_TYPE_UNKNOWN;
    }

    H5O_info_t infobuf;
    if (H5Oget_info_by_name(_fileId, path.c_str(), &infobuf, H5P_DEFAULT)
        < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getHDF5ObjType("
            << path.c_str() << "):"
            << " failed to get object information";
        return H5O_TYPE_UNKNOWN;
    }
    return infobuf.type;
}

hid_t HDF5::__getDatasetId(const std::string &datasetName)
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getDatasetId("
            << datasetName.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return -1;
    }
    if ( __getHDF5ObjType(datasetName) != H5O_TYPE_DATASET ) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- " << "HDF5::__getDatasetId("
            << datasetName.c_str() << "):"
            << " dataset does not exist";
        return -1;
    }

    hid_t _curDatasetId = H5Dopen(_fileId, datasetName.c_str(), H5P_DEFAULT);
    if (_curDatasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getDatasetId("
            << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return -1;
    } else {
        return _curDatasetId;
    }
} // HDF5::__getDatasetId

bool HDF5::__getDatasetType(const std::string &datasetName,
                            FQ::DataType *fqType)
{
    *fqType = FQ::FQT_UNKNOWN;

    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getDatasetType("
            << datasetName.c_str() << "):"
            << " cannot open the dataset";
        return false;
    }

    hid_t type = H5Dget_type(datasetId);
    H5Dclose(datasetId);
    *fqType = __getFQDataType(type);
    if (fqType < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getDatasetType("
            << datasetName.c_str() << "):"
            << " failed to get the HDF5 data type";
        return false;
    }
    return true;
} // HDF5::__gtDatasetType

bool HDF5::__getDatasetSize(const std::string &datasetName, uint64_t *size)
{
    *size = 0;
    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 1) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- HDF5::__getDatasetSize"
            << "(" << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }
    *size = (uint64_t)H5Dget_storage_size(datasetId);
    if (*size == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- HDF5::__getDatasetSize"
            << "(" << datasetName.c_str() << "):"
            << " dataset size is 0";
        return false;
    }
    return true;
} // HDF5::__getDatasetSize

bool HDF5::__getDatasetDimension(const std::string &datasetName,
                                 std::vector <uint64_t> &dims)
{
    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- " << "HDF5::__getDatasetDimension"
            << "(" << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }

    hid_t spaceId = H5Dget_space(datasetId);
    H5Dclose(datasetId);
    if (spaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getDatasetDimension"
            << "(" << datasetName.c_str() << "):"
            << " failed to open the dataset space";
        return false;
    }

    int numDim = H5Sget_simple_extent_ndims(spaceId);
    if (numDim < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getDatasetDimension"
            << "(" << datasetName.c_str() << "):"
            << " failed to get the number of dimensions";
        H5Sclose(spaceId);
        return false;
    }

    hsize_t* stddims = new hsize_t[numDim];
    hsize_t* maxdims = new hsize_t[numDim];
    int ierr = H5Sget_simple_extent_dims(spaceId, stddims, maxdims);
    H5Sclose(spaceId);
    if (ierr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getDatasetDimension"
            << "(" << datasetName.c_str() << "):"
            << " failed to get the dimension size";
        delete[] stddims;
        delete[] maxdims;
        return false;
    }

    //now read the attributes for this data set to get the number of dim
    for(int i=0; i<numDim; i++){
        dims.push_back(stddims[i]);
    }
    delete[] stddims;
    delete[] maxdims;
    return true;
} // HDF5::__getDatasetDimension

bool HDF5::__getAttribute(const std::string &path,
                          const std::string &attrName, void *values)
{
    H5O_type_t objType = __getHDF5ObjType(path);
    if(objType != H5O_TYPE_DATASET && objType != H5O_TYPE_GROUP) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttribute"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " not a valid path to a dataset or group";
        return false;
    }

    hid_t objectId = H5Oopen(_fileId, path.c_str(), H5P_DEFAULT);
    if (objectId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttribute"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " failed to open the object";
        return false;
    }

    if ( H5Aexists(objectId, attrName.c_str()) <= 0) {
        // attribute does not exist
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttribute"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " attribute does not exist";
        H5Oclose(objectId);
        return false;
    }

    hid_t attrId = H5Aopen_name(objectId, attrName.c_str());
    if (attrId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttribute"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " failed to open the attribute";
        H5Oclose(objectId);
        return false;
    }
    hid_t attrType = H5Aget_type(attrId);
    if ( H5Aread(attrId, attrType, values) < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttribute"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " failed to read the attribute";
        H5Aclose(attrId);
        H5Oclose(objectId);
        return false;
    }

    H5Aclose(attrId);
    H5Oclose(objectId);
    return true;
} // HDF5::__getAttribute

bool HDF5::__setAttribute(const std::string &path,
                          const std::string &attrName,
                          const void *values, const hsize_t len,
                          const FQ::DataType fqType)
{
    int berr = true;
#ifndef FQ_NOMPI
    if(mpi_rank==0) {
#endif
        if (len == 0) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- HDF5::__setAttribute"
                << "(" << path.c_str() << "," << attrName.c_str() << "):"
                << " cannot create an attribute with size 0";
            berr = false;
        }

        if (berr) {
            H5O_type_t objType = __getHDF5ObjType(path);
            if(objType != H5O_TYPE_DATASET && objType != H5O_TYPE_GROUP) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::__setAttribute"
                    << "(" << path.c_str() << "," << attrName.c_str() << "):"
                    << " not a valid path to a dataset or group";
                berr = false;
            }
        }

        // delete the existing attribute
        hid_t objId;
        if (berr) {
            objId = H5Oopen(_fileId, path.c_str(), H5P_DEFAULT);
            if ( H5Aexists(objId, attrName.c_str()) ) {
                if (H5Adelete(objId, attrName.c_str()) < 0) {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- HDF5::__setAttribute"
                        << "(" << path.c_str() << "," << attrName.c_str() << "):"
                        << " failed to delete the existing attribute";
                    berr = false;
                }
            }
        }

        hid_t dataspaceId;
        if (berr) {
            dataspaceId = H5Screate_simple(1, &len, NULL);
            if (dataspaceId < 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::__setAttribute"
                    << "(" << path.c_str() << "," << attrName.c_str() << "):"
                    << " failed to create the attribute dataspace";
                berr = false;
            }
        }

        hid_t attrType;
        if (berr) {
            attrType = __getHDF5DataType(fqType);
            if (attrType < 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::__setAttribute"
                    << "(" << path.c_str() << "," << attrName.c_str() << "):"
                    << " not yet supported FQ data type " << fqType;
                berr = false;
            }
        }

        hid_t attrId;
        if (berr) {
            attrId = H5Acreate(objId, attrName.c_str(), attrType, dataspaceId,
                               H5P_DEFAULT, H5P_DEFAULT);
            if (attrId < 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::__setAttribute"
                    << "(" << path.c_str() << "," << attrName.c_str() << "):"
                    << " failed to create the attribute";
                berr = false;
            }
        }

        if (H5Awrite(attrId, attrType, values) < 0) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- HDF5::__setAttribute"
                << "(" << path.c_str() << "," << attrName.c_str() << "):"
                << " failed to write the attribute";
            berr = false;
        }
        H5Oclose(objId);
        H5Sclose(dataspaceId);
        H5Aclose(attrId);
#ifndef FQ_NOMPI
    }
    MPI_Bcast(&berr, 1, MPI_BYTE, 0, mpi_comm);
#endif
#if (! defined(FQ_NOMPI)) && defined(HDF5_DISABLE_METADATA_FLUSH)
    H5Fflush(_fileId, H5F_SCOPE_LOCAL);
#endif
    return berr;
} // HDF5::__setAttribute

bool HDF5::__getAttributeInfo(const std::string &path,
                              const std::string &attrName,
                              uint64_t *size, FQ::DataType *fqType)
{
    H5O_type_t objType = __getHDF5ObjType(path);
    if(objType != H5O_TYPE_DATASET && objType != H5O_TYPE_GROUP) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttributeInfo"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " not a valid path to a dataset or group";
        return false;
    }

    hid_t objectId = H5Oopen(_fileId, path.c_str(), H5P_DEFAULT);
    if (objectId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttributeInfo"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " failed to open the object";
        return false;
    }

    if ( H5Aexists(objectId, attrName.c_str()) <= 0) {
        // attribute does not exist
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttributeInfo"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " attribute does not exist";
        H5Oclose(objectId);
        return false;
    }

    hid_t attrId = H5Aopen_name(objectId, attrName.c_str());
    if (attrId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttributeInfo"
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " failed to open the attribute";
        H5Aclose(attrId);
        H5Oclose(objectId);
        return false;
    }
    hid_t type = H5Aget_type(attrId);
    *fqType = __getFQDataType(type);
    if (*fqType == FQ::FQT_UNKNOWN) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAttributeInfo("
            << "(" << path.c_str() << "," << attrName.c_str() << "):"
            << " failed to get the attribute type";
        return false;
    }
    *size = H5Aget_storage_size(attrId);

    H5Aclose(attrId);
    H5Oclose(objectId);
    return true;
} // HDF5::__getAttributeInfo

bool HDF5::__readData(const std::string &datasetName, void *data)
{
    ibis::horometer openTimer;
    openTimer.start();

    if (data == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- HDF5::__readData("
            << datasetName.c_str() << "):"
            << " no data needs to be read";
        return true;
    }
    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readData("
            << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }
    hid_t type = H5Dget_type(datasetId);
    if (type < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readData("
            << datasetName.c_str() << "):"
            << " unknown HDF5 data type";
        return false;
    }
    openTimer.stop();

    LOGGER(FastQuery::reportTiming())
        << "Statistic\thdf5File::__readData:H5Dopen\t"
        << openTimer.CPUTime() << "\t" << openTimer.realTime() << "\t";

    ibis::horometer readTimer;
    readTimer.start();

    /*#ifndef FQ_NOMPI
      hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      herr_t herr = H5Dread(datasetId, type, H5S_ALL, H5S_ALL, xfer_plist, data);
      #else*/
    herr_t herr = H5Dread(datasetId, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          data);
    //#endif

    readTimer.stop();

    LOGGER(FastQuery::reportTiming())
        << "Statistic\thdf5File::__readData:H5Dread\t"
        << readTimer.CPUTime() << "\t" << readTimer.realTime() << "\t";

    ibis::horometer closeTimer;
    closeTimer.start();


    H5Dclose(datasetId);

    closeTimer.stop();

    LOGGER(FastQuery::reportTiming())
        << "Statistic\thdf5File::__readData:H5Dclose\t"
        << closeTimer.CPUTime() << "\t" << closeTimer.realTime() << "\t";

    if (herr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readData("
            << datasetName.c_str() << "):"
            << " failed to read data";
        return false;
    }
    return true;
} // HDF5::__readData

bool HDF5::__writeData(const std::string &datasetName, const void *data)
{
    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__writeData("
            << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }
    hid_t type = H5Dget_type(datasetId);
    if (type < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__writeData("
            << datasetName.c_str() << "):"
            << " unknown HDF5 data type";
        return false;
    }

    /*#ifndef FQ_NOMPI
      hid_t xfer_plist = H5P_DEFAULT;
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      #endif*/

    herr_t herr = H5Dwrite(datasetId, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    /*#ifndef FQ_NOMPI
      H5Pclose(xfer_plist);
      #endif*/
    H5Dclose(datasetId);
    if (herr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__writeData("
            << datasetName.c_str() << "):"
            << " failed to write data";
        return false;
    }
    return true;
} // HDF5::__writeData

bool HDF5::__createDataset(const std::string& datasetName,
                           const std::vector<uint64_t> dims,
                           const FQ::DataType fqType, const bool UNLIMITED)
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__createDatset("
            << datasetName.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false;
    }

    if (dims.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__createDataset("
            << datasetName.c_str() << "):"
            << " number of dimension is 0";
        return false;
    }

    if (datasetName.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::createDataset("
            << datasetName.c_str() << "):"
            << " dataset name is not given";
        return false;
    }

    if (datasetName[0] != '/') {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::createDataset("
            << datasetName.c_str() << "):"
            << " dataset name  must be given as an absolute path starting with '/'";
        return false;
    }

    std::string pathName = "";
    std::string varName = "";
    size_t pos = datasetName.find_last_of('/');
    if (pos != 0) {
        pathName = datasetName.substr(0, pos);
        varName = datasetName.substr(pos+1);
    } else {
        pathName = "";
        varName = datasetName;
    }

    //check if dataset already exist
    H5O_type_t objType = __getHDF5ObjType(datasetName);
    if (objType == H5O_TYPE_DATASET) {
        // remove the existing dataset
        __removeDataset(datasetName);
    } else {
        // create the intermediate groups to the dataset
        if (pathName.compare("") != 0) {
            if ( __getHDF5ObjType(pathName.c_str()) != H5O_TYPE_GROUP) {
                hid_t grp_crt_plist = H5Pcreate(H5P_LINK_CREATE);
                /* Set flag for intermediate group creation */
                H5Pset_create_intermediate_group(grp_crt_plist, true);
                hid_t groupId = H5Gcreate2(_fileId, pathName.c_str(), grp_crt_plist, H5P_DEFAULT, H5P_DEFAULT);
                if (groupId < 0) {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- HDF5::__createDataset("
                        << datasetName.c_str() << "):"
                        << " fail to create intermediate group for the dataset";
                    H5Pclose(grp_crt_plist);
                    return false;
                }
                H5Gclose(groupId);
                H5Pclose(grp_crt_plist);
            }
        }
    }

    hid_t type = __getHDF5DataType(fqType);
    if (type < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__createDataset("
            << datasetName.c_str() << "):"
            << " not yet supported FQ data type " << fqType;
        return false;
    }

    uint64_t nDims = dims.size();
    hsize_t *size = new hsize_t[nDims];
    uint64_t totSize = 1;
    uint64_t totLen = dims[0];
    for(uint64_t i=0; i<nDims; i++) {
        size[i] = dims[i];
        totSize *= dims[i];
    }

    bool CHUNKED = false;
    hsize_t *maxdims;

    if (UNLIMITED) {
        CHUNKED = true;
    }
    if (totSize > FQ_CHUNK_SIZE) {
        CHUNKED = true;
	size[0] = floor(FQ_CHUNK_SIZE/(totSize/dims[0]));
    }

    if (CHUNKED) {
        maxdims = new hsize_t[nDims];
        for(uint64_t i=0; i<nDims; i++) maxdims[i] = H5S_UNLIMITED;
    } else {
        maxdims = NULL;
    }

    hid_t spaceId = H5Screate_simple(nDims, size, maxdims);
    delete[] maxdims;
    if (spaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__createDataset("
            << datasetName.c_str() << "):"
            << " failed to create memory space with dimension " << nDims
            << " and size " << totSize;
        delete[] size;
        return false;
    }

    hid_t plist = H5P_DEFAULT;
    if (CHUNKED) {
        plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist, dims.size(), size);
    }

    hid_t datasetId;
    if (pathName.compare("") == 0) {
        // create dataset under root with dataset name only
        datasetId = H5Dcreate(_fileId, varName.c_str(), type, spaceId,
                              H5P_DEFAULT, plist, H5P_DEFAULT);
    } else {
        datasetId = H5Dcreate(_fileId, datasetName.c_str(), type, spaceId,
                              H5P_DEFAULT, plist, H5P_DEFAULT);
    }
    if (CHUNKED) {
        H5Pclose(plist);
    }
    H5Sclose(spaceId);
    H5Dclose(datasetId);

    if (size[0] != totLen) {
        std::vector<uint64_t> curDims = dims;
        curDims[0] = (uint64_t)size[0];
        if(! __extendDataset(datasetName.c_str(), curDims, totLen)) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- HDF5::__createDataset("
                << varName.c_str() << "):"
                << " failed to extend the dataset";
            delete[] size;
            return false;
        }
    }
    delete[] size;

    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__createDataset("
            << datasetName.c_str() << "):"
            << " failed to create the dataset with size " << totSize;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__createDataset("
        << datasetName.c_str() << "):"
        << " successfully created the dataset with size " << totSize;
#if (! defined(FQ_NOMPI)) && defined(HDF5_DISABLE_METADATA_FLUSH)
    H5Fflush(_fileId, H5F_SCOPE_LOCAL);
#endif
    return true;
} // HDF5::__createDataset

bool HDF5::__removeDataset(const std::string &datasetName)
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__removeDatset("
            << datasetName.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false;
    }

    H5O_type_t objType = __getHDF5ObjType(datasetName);
    if (objType != H5O_TYPE_DATASET) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__removeDataset("
            << datasetName.c_str() << "):"
            << " dataset does not exist";
        return false;
    }
    if ( H5Ldelete(_fileId, datasetName.c_str(), H5P_DEFAULT) < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__removeDataset("
            << datasetName.c_str() << "):"
            << " failed to remove the dataset link";
        return false;
    }
#if (! defined(FQ_NOMPI)) && defined(HDF5_DISABLE_METADATA_FLUSH)
    H5Fflush(_fileId, H5F_SCOPE_LOCAL);
#endif
    return true;
} // HDF5::__removeDataset

bool HDF5::__readPointData(const std::string &datasetName, const std::vector<uint64_t> &coords, void *data)
{
    if (coords.size() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " no element is selected";
        return true;
    }


    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }

    hid_t spaceId = H5Dget_space(datasetId);
    if (spaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " failed to open the dataspace";
        H5Dclose(datasetId);
        return false;
    }

    int numDim = H5Sget_simple_extent_ndims(spaceId);
    if (numDim < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " failed to get dataspace dimension";
        H5Dclose(datasetId);
        H5Sclose(spaceId);
        return false;
    }

    hsize_t nElements = coords.size() / numDim;
    hsize_t* indices = new hsize_t[coords.size()];
    for (unsigned int i=0; i<coords.size(); i++) {
        indices[i] = coords[i];
    }
    herr_t herr = H5Sselect_elements(spaceId, H5S_SELECT_SET, nElements, indices);
    if (herr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " failed to select data";
        delete[] indices;
        H5Dclose(datasetId);
        H5Sclose(spaceId);
        return false;
    }
    delete[] indices;

    herr = H5Sselect_valid(spaceId);
    if (herr == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " not a valid selection out of dataspace range";
        H5Sclose(spaceId);
        H5Dclose(datasetId);
        return false;
    }
    hid_t memspaceId = H5Screate_simple(1, &nElements, NULL);
    if (memspaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " failed to create memory space";
        H5Sclose(spaceId);
        H5Dclose(datasetId);
        return false;
    }
    herr = H5Sselect_all(memspaceId);
    hid_t type = H5Dget_type(datasetId);

    /*#ifndef FQ_NOMPI
      hid_t xfer_plist = H5P_DEFAULT;
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
      #endif*/

    herr = H5Dread(datasetId, type, memspaceId, spaceId, H5P_DEFAULT, data);

    /*#ifndef FQ_NOMPI
      H5Pclose(xfer_plist);
      #endif*/

    H5Dclose(datasetId);
    H5Sclose(spaceId);
    H5Sclose(memspaceId);

    if (herr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readPointData("
            << datasetName.c_str() << "):"
            << " failed to read data";
        return false;
    }
    return true;
} // HDF5::__readPointData

bool HDF5::__readArrayData(const std::string& datasetName, const std::vector<uint64_t> &offsets,
                           const std::vector<uint64_t> &counts, const std::vector<uint64_t> &strides, void *data)
{
    if (counts.size() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " data count is empty";
        return true;
    }

    uint64_t totalElements = counts[0];
    for (unsigned int i=1; i<counts.size(); i++) {
        totalElements*=counts[i];
    }
    if (totalElements == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " size of total data count is 0";
        return true;
    }

    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " cannot open the dataset";
        return false;
    }

    hid_t spaceId = H5Dget_space(datasetId);
    if (spaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " cannot open the dataspace";
        H5Dclose(datasetId);
        return false;
    }

    if( H5Sselect_hyperslab(spaceId, H5S_SELECT_SET, (hsize_t*)(&offsets[0]), (hsize_t*)(&strides[0]), (hsize_t*)(&counts[0]), NULL) < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " cannot select the data using hyperslab";
        H5Dclose(datasetId);
        H5Sclose(spaceId);
        return false;
    }
    hid_t memspaceId = H5Screate_simple(counts.size(), (hsize_t*)(&counts[0]), NULL);
    if (memspaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " cannot create memory space with size " << counts.size();
        H5Dclose(datasetId);
        H5Sclose(spaceId);
        return false;
    }
    H5Sselect_all(memspaceId);
    hid_t type = H5Dget_type(datasetId);

    hid_t xfer_plist = H5P_DEFAULT;

#ifndef FQ_NOMPI
    if (1) {
        xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    }
#endif

    hid_t herr = H5Dread(datasetId, type, memspaceId, spaceId, H5P_DEFAULT, data);

#ifndef FQ_NOMPI
    H5Pclose(xfer_plist);
#endif

    H5Dclose(datasetId);
    H5Sclose(spaceId);
    H5Sclose(memspaceId);

    if ( herr < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__readArrayData "
            << "(" << datasetName.c_str() << "):"
            << " failed to read data";
        return false;
    }
    return true;
} // HDF5::__readArrayData

bool HDF5::__writeArrayData(const std::string& datasetName, const std::vector<uint64_t> &offsets,
                            const std::vector<uint64_t> &counts, const std::vector<uint64_t> &strides, const void *data, const bool collective)
{
    if (counts.size() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " data count is empty";
        return true;
    }

    uint64_t totalElements = counts[0];
    for (unsigned int i=1; i<counts.size(); i++) {
        totalElements*=counts[i];
    }
    if (totalElements == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " size of total data count is 0";
        return true;
    }

    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }

    hid_t spaceId = H5Dget_space(datasetId);
    if (spaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " failed to open the dataspace";
        H5Dclose(datasetId);
        return false;
    }

    hsize_t* offsetsPtr = (hsize_t*)&(offsets[0]);
    hsize_t* stridesPtr = (hsize_t*)&(strides[0]);
    hsize_t* countsPtr = (hsize_t*)&(counts[0]);

    if( H5Sselect_hyperslab(spaceId, H5S_SELECT_SET, (hsize_t*)(&offsets[0]), (hsize_t*)(&strides[0]), (hsize_t*)(&counts[0]), NULL) < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " invalid dataset space selection";
        H5Dclose(datasetId);
        H5Sclose(spaceId);
        return false;
    }
    hid_t memspaceId = H5Screate_simple(counts.size(), (hsize_t*)(&counts[0]), NULL);
    if (memspaceId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " failed to create the memory space with size " << counts.size();
        H5Dclose(datasetId);
        H5Sclose(memspaceId);
        return false;
    }
    if ( H5Sselect_all(memspaceId) < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " invalid memory space selection";
        H5Dclose(datasetId);
        H5Sclose(spaceId);
        H5Sclose(memspaceId);
        return false;
    }

    hid_t type = H5Dget_type(datasetId);

    hid_t xfer_plist = H5P_DEFAULT;
#ifndef FQ_NOMPI
    if (collective) {
        xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    }
#endif

    hid_t herr = H5Dwrite(datasetId, type, memspaceId, spaceId, H5P_DEFAULT, data);

#ifndef FQ_NOMPI
    H5Pclose(xfer_plist);
#endif

    H5Dclose(datasetId);
    H5Sclose(spaceId);
    H5Sclose(memspaceId);
    if ( herr < 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__writeArrayData "
            << "(" << datasetName.c_str() << "):"
            << " failed to write array data with size " << counts.size();
        return false;
    }
    return true;
} // HDF5::__writeArrayData

bool HDF5::__extendDataset(const std::string& datasetName, const std::vector<uint64_t> &curDims, const uint64_t totLen)
{
    hid_t datasetId = __getDatasetId(datasetName);
    if (datasetId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__extendDataset "
            << "(" << datasetName.c_str() << "):"
            << " failed to open the dataset";
        return false;
    }
    if (curDims.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__extendDataset "
            << "(" << datasetName.c_str() << "):"
            << " number of dimension is 0";
        H5Dclose(datasetId);
        return false;
    }

    // extend the dataset by maximum of 4GB at a time
    unsigned int step = FQ_CHUNK_SIZE;
    hsize_t *size = new hsize_t[curDims.size()];
    for (int i=0; i<curDims.size(); i++) size[i] = (hsize_t)curDims[i];
    size[0] = curDims[0];
    while(size[0] < totLen) {
        if (size[0]+step > totLen) {
            size[0] += (totLen - size[0]);
        } else {
            size[0] += step;
        }
        if (H5Dset_extent(datasetId, size) < 0) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "HDF5::__extendDataset "
                << "(" << datasetName.c_str() << "):"
                << " failed to extend the dataset";
            delete[] size;
            H5Dclose(datasetId);
            return false;
        }
    }
    delete[] size;
    H5Dclose(datasetId);

    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__extendDataset "
        << "(" << datasetName.c_str() << "):"
        << " successfully extend the dataset";
#if (! defined(FQ_NOMPI)) && defined(HDF5_DISABLE_METADATA_FLUSH)
    //    H5Fflush(_fileId, H5F_SCOPE_LOCAL);
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__extendDataset "
        << "(" << datasetName.c_str() << "):"
        << " flush metadata";
#endif
    return true;
} // HDF5::__extendDataset

#ifndef FQ_NOMPI
bool HDF5::__getOffsets(const std::string &variable, const int column,
                        uint64_t *start, uint64_t *end, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIoffsetTable";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.push_back(mpi_idx);
    counts.push_back(2);  //read the start and end offset
    strides.push_back(1);
    offsets.push_back(column);
    counts.push_back(1);
    strides.push_back(1);
    uint64_t vals[2];
    if (! __readArrayData(datasetName, offsets, counts, strides, vals)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getOffsets("
            << variable.c_str() << "):"
            << " failed to get offset table";
        return false;
    }
    if (vals[1] < vals[0]) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getOffsets("
            << variable.c_str() << "):"
            << " invalid offset";
        return false;
    }
    *start = vals[0];
    *end = vals[1];
    LOGGER(ibis::gVerbose > 2)
        << "HDF5::__getOffsets("
        << variable.c_str() << "):"
        << " successfully got offset at position "
        << offsets[0] << " with values "
        << "("<<*start<<"-"<<*end<<"):";
    return true;
} // HDF5::__getOffsets

bool HDF5::__getOffsetLength(const std::string &variable, const int column,
                             uint64_t* len, const uint64_t mpi_idx)
{
    uint64_t start;
    uint64_t end;
    if (! __getOffsets(variable, column, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getOffsetLength("
            << variable.c_str() << "):"
            << " failed to get offset table";
        return false;
    }
    *len = (end-start);
    return true;
} // HDF5::__getOffsetLength
#endif

bool HDF5::__getAllVariables(const std::string &path, std::vector<std::string> &variables)
{
    variables.clear();

    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- HDF5::__getAllVariables("
            << path.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false;
    }
    // if the path is not valid file location, return false
    if (path.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getAllVariables("
            << path.c_str() << "):"
            << " path is an empty string";
        return false;
    }
    if (path[0] != '/') {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getAllVariables("
            << path.c_str() << "):"
            << " not an absolute path";
        return false;
    }

    H5O_type_t objType = __getHDF5ObjType(path);
    if (objType == H5O_TYPE_DATASET) {
        // if the path refers to a dataset, add the variable.
        variables.push_back(path);
        return true;
    } else if (objType == H5O_TYPE_GROUP || path.compare("/") == 0) {
        // if the path refers to a group, add all variables under the group.
        Node* rootNode = __getVariableTree(path);
        if (rootNode != NULL) {
            // traverse the object tree and collect datasets into the variables vector
            std::stack<Node*> s;
            s.push(rootNode);
            while(! s.empty()) {
                Node* node = s.top();
                s.pop();
                if (node->isDataset() == true) {
                    std::string varName;
                    Node* tmp = node;
                    while(1) {
                        varName.insert(0, tmp->getName());
                        tmp = tmp->getParent();
                        if (tmp == NULL) break;
                        varName.insert(0, "/");
                        if (tmp->getName().compare("/") == 0) break;
                    }
                    variables.push_back(varName);
                } else {
                    for (unsigned int i=0; i<node->getNumChildren(); i++) {
                        Node* childNode = node->getChild(i);
                        s.push(childNode);
                    }
                }
            }
            // free memory space
            delete(rootNode);
        }
        return true;
    } else {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "HDF5::__getAllVariables("
            << path.c_str() << "):"
            << " path does not refer to a dataset or group";
        return false;
    }
} // HDF5::__getAllVariables

Node* HDF5::__getVariableTree(const std::string &path)
{
    H5O_info_t infobuf;

    hid_t root = H5Gopen(_fileId, path.c_str(), H5P_DEFAULT);
    if (root < 0) return NULL;
    herr_t status = H5Oget_info (root, &infobuf);
    Node* node = new Node(NULL, infobuf.addr, path.c_str(), false);
    status = H5Literate (root, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func, (void *) node);
    status = H5Gclose (root);
    return node;
} // __getVariableTree


/************************************************************

  Operator function.  This function prints the name and type
  of the object passed to it.  If the object is a group, it
  is first checked against other groups in its path using
  the group_check function, then if it is not a duplicate,
  H5Literate is called for that group.  This guarantees that
  the program will not enter infinite recursion due to a
  circular path in the file.

************************************************************/
herr_t op_func (hid_t loc_id, const char *name, const H5L_info_t *info,
                void *operator_data)
{
    herr_t return_val = 0;
    H5O_info_t infobuf;
    Node* node = (Node*) operator_data;
    H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);
    switch (infobuf.type) {
    case H5O_TYPE_GROUP:
        /*
         * Check group address against linked list of operator
         * data structures.  We will always run the check, as the
         * reference count cannot be relied upon if there are
         * symbolic links, and H5Oget_info_by_name always follows
         * symbolic links.  Alternatively we could use H5Lget_info
         * and never recurse on groups discovered by symbolic
         * links, however it could still fail if an object's
         * reference count was manually manipulated with
         * H5Odecr_refcount.
         */
        if ( group_check (node, infobuf.addr) ) {
            printf ("Warning: Loop detected!\n");
        } else {
            /*
             * Initialize new operator data structure and begin
             * recursive iteration on the discovered group
             */
            Node* childNode = new Node(node, infobuf.addr, name, false);
            node->addChild(childNode);
            return_val = H5Literate_by_name (loc_id, name, H5_INDEX_NAME,
                                             H5_ITER_NATIVE, NULL, op_func, (void *) childNode,
                                             H5P_DEFAULT);
        }
        break;
    case H5O_TYPE_DATASET: {
        std::string str = name;
        if (str.find(".bitmap") == str.npos &&
            str.find(".bitmapKeys") == str.npos &&
            str.find(".bitmapOffsets") == str.npos ) {
            Node* childNode = new Node(node, infobuf.addr, name, true);
            node->addChild(childNode);
        }
        break;}
    case H5O_TYPE_NAMED_DATATYPE:
        break;
    default:
        break;
    }
    return return_val;
} // op_func


/************************************************************

  This function recursively searches the linked list of
  opdata structures for one whose address matches
  target_addr.  Returns 1 if a match is found, and 0
  otherwise.

************************************************************/
int group_check (Node *node, haddr_t target_addr)
{
    if (node->getAddr() == target_addr)
        return 1;       /* Addresses match */
    else if (node->getParent() == NULL)
        return 0;       /* Root group reached with no matches */
    else
        return group_check (node->getParent(), target_addr);
    /* Recursively examine the next node */
} // group_check
