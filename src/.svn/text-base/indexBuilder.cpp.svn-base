#include "indexBuilder.h"
#include <cmath>

#ifdef FQ_NOMPI
IndexBuilder::IndexBuilder(const std::string &dataFileName,
                           const FQ::FileFormat dataModel,
                           const std::string &indexFileName,
                           const int v,
                           const char *rcfile,
                           const char *logfile,
                           void *extra)
    : FastQuery(dataFileName, dataModel, indexFileName, v, rcfile, logfile,
                false, extra)
#else
IndexBuilder::IndexBuilder(const std::string &dataFileName,
                           const FQ::FileFormat dataModel,
                           const std::string &indexFileName,
                           const int v,
                           const char *rcfile,
                           const char *logfile,
                           MPI_Comm comm,
                           void *extra)
      : FastQuery(dataFileName, dataModel, indexFileName, v, rcfile, logfile,
                  false, comm, extra)
#endif
{
} // IndexBuilder::IndexBuilder

int IndexBuilder::buildIndexes(const char *binning,
                               const std::string &varPathStr,
                               const std::string &varNameStr,
                               uint64_t mpi_dim,
                               uint64_t mpi_len)
{
    int ret = 0;
    if (! isValid("IndexBuilder::buildIndexes")) return ret;
    ibis::horometer timer;
    timer.start();
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);

    if (numVar < 1) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::buildIndexes "
            " failed to find any variable named \"" << varNameStr
            << "\" under path \"" << varPathStr.c_str() << "\" in file \""
            << metadataMgr->getFileName() << '"';
        return ret;
    }
    // build one variable at a time
#ifdef FQ_NOMPI
    for(unsigned int i=0; i<varInfoList.size(); i++) {
        std::vector<VarInfo> newVarInfoList;
        newVarInfoList.push_back(varInfoList[i]);
        std::vector<VarSpace> newVarSpaceList;
        newVarSpaceList.push_back(varSpaceList[i]);
        FQ_Part part(newVarInfoList, newVarSpaceList, *dataFile, *indexFile);
        ret += part.buildIndexes(*indexFile, binning);
    }
#else
    int val = 0;
    if (mpi_len == 0)
        mpi_len = FQ_DEFAULT_MPI_LEN;
    do { // while (nextStep())
        LOGGER(ibis::gVerbose > 1)
            << "indexBulder::buildIndexes starting to process step "
            << dataFile->currentStep() << " from " << dataFile->getFileName();

        for(unsigned int i=0; i<varInfoList.size(); i++) {
            ibis::horometer varTimer;
            varTimer.start();
            std::vector<VarInfo> newVarInfoList;
            newVarInfoList.push_back(varInfoList[i]);
            std::vector<VarSpace> newVarSpaceList;
            bool flag = true;

            std::vector<uint64_t> offsets = varSpaceList[i].getOffsets();
            std::vector<uint64_t> counts = varSpaceList[i].getCounts();
            std::vector<uint64_t> strides = varSpaceList[i].getStrides();

            uint64_t totalCount = varSpaceList[i].getSize();
            // if (totalCount <= (mpi_size-1)*mpi_len) {
            //     LOGGER(ibis::gVerbose > 0)
            //         << "Warning -- IndexBuilder::buildIndexes:"
            //         " (mpi_size-1)*mpi_len must be larger than"
            //         " the dimension length " << totalCount;
            //     return -1;
            // }
            uint64_t mpi_max_iter =
                std::ceil((double)totalCount/((double)mpi_len*mpi_size));
            std::string mpi_text = varSpaceList[i].getText();
            LOGGER(ibis::gVerbose > 2)
                << "Debug --- indexBuilder::buildIndexes: totalCount="
                << totalCount << ", mpi_max_iter=" << mpi_max_iter
                << ", mpi_size=" << mpi_size
                << ", mpi_len=" << mpi_len
                << ", mpi_dim=" << mpi_dim;

            for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
                uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
                std::vector<uint64_t> splitOffsets = offsets;
                std::vector<uint64_t> splitCounts = counts;
                splitOffsets[mpi_dim] = mpi_idx*mpi_len*strides[mpi_dim] +
                    offsets[mpi_dim];
                if (mpi_idx*mpi_len >= totalCount) {
                    // out of the range
                    splitCounts[mpi_dim] = 0;
                }
                else if ((mpi_idx+1)*mpi_len>totalCount &&
                         (totalCount%mpi_len)!= 0) {
                    // the last array may not have the full length
                    splitCounts[mpi_dim] = totalCount%mpi_len;
                } else {
                    splitCounts[mpi_dim] = mpi_len;
                }
                VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
                splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                         mpi_iter, mpi_max_iter, mpi_comm);
                newVarSpaceList.clear();
                newVarSpaceList.push_back(splitVarSpace);
                FQ_Part part(newVarInfoList, newVarSpaceList,
                             *dataFile, *indexFile);
                flag = (part.buildIndexes(*indexFile, binning) == 1);
            }
            varTimer.stop();
            if (flag) val++;
            LOGGER(FastQuery::reportTiming())
                << "Statistic\tIndexBuilder::buildIndex\t"
                << varTimer.CPUTime() << "\t" << varTimer.realTime()
                << "\t" << i;
        }
    } while (dataFile->nextStep());
    MPI_Allreduce(&val, &ret, 1, MPI_INT, MPI_MIN, mpi_comm);
#endif
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tIndexBuilder::buildIndexes\t"
        << timer.CPUTime() << "\t" << timer.realTime();
    LOGGER(ibis::gVerbose > 2)
        << "IndexBuilder::buildIndexes: successfully build indixes for "
        << ret << " variables used " << timer.CPUTime() << " sec CPU time and "
        << timer.realTime() << " sec elapsed time";
    return ret;
} // IndexBuilder::buildIndexes


bool IndexBuilder::createNewVariable
(const std::string &variable,
 std::vector<uint64_t> dims, FQ::DataType type)
{
    if (! isValid("IndexBuilder::createNewVariable")) return false;
    if (! dataFile->createDataset(variable, dims, type)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::createNewVariable:"
            << " failed to create the dataset for the new variable \""
            << variable.c_str() << "\"";
        return false;
    }
    if (! metadataMgr->addVariable(variable)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::createNewVariable:"
            << " failed to register the new variable \""
            << variable.c_str() << "\"";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "IndexBuilder::createNewVariable:"
            << " successfully create new variable \""
            << variable.c_str() << "\"";
        return true;
    }
} // IndexBuilder

bool IndexBuilder::setData(const std::string &varNameStr, const void* data,
                           const std::string &varPathStr, const bool collective)
{
    if (! isValid("IndexBuilder::setData")) return false;
    bool berr = true;
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found from the file";
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
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- " << "IndexBuilder::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList[idx].getPath() << "\"";
    }
    VarInfo varInfo = varInfoList[idx];
    VarSpace varSpace = varSpaceList[idx];
#ifndef FQ_NOMPI
    if (collective==false) {
#endif
    if (varInfo.getSize() == varSpace.getSize()) {
        berr = dataFile->setData(varInfo.getPath(), data);
    } else {
        berr = dataFile->setArrayData
            (varInfo.getPath(), varSpace.getOffsets(),
             varSpace.getCounts(), varSpace.getStrides(), data);
    }
#ifndef FQ_NOMPI
    } else {
    int mpi_dim = 0; // split data by the first dimension
    std::vector<uint64_t> offsets = varSpace.getOffsets();
    std::vector<uint64_t> counts = varSpace.getCounts();
    std::vector<uint64_t> strides = varSpace.getStrides();
    uint64_t totalCount = counts[mpi_dim];
    if (totalCount < mpi_size) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " number of MPI tasks cannot be less than"
            << " the dimension length " << totalCount;
        return false;
    }
    int mpi_len = ceil((double)totalCount/(double)mpi_size);
    uint64_t data_pos;
    if (mpi_rank >= totalCount%mpi_len && totalCount%mpi_len != 0) {
        counts[mpi_dim] = mpi_len - 1;
        offsets[mpi_dim] += (mpi_len*mpi_rank -
                             (mpi_rank-totalCount%mpi_len))*strides[mpi_dim];
        data_pos = (mpi_len*mpi_rank - (mpi_rank-totalCount%mpi_len))
            *(varSpace.getSize() / totalCount);
    } else {
        counts[mpi_dim] = mpi_len;
        offsets[mpi_dim] += (mpi_len*mpi_rank)*strides[mpi_dim];
        data_pos = (mpi_len*mpi_rank)*(varSpace.getSize()/totalCount);
    }
    bool mpi_berr = false;
    FQ::DataType fqType = varInfo.getType();
    switch(fqType){
    case FQ::FQT_FLOAT: {
        float* mpi_data = (float*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_DOUBLE: {
        double* mpi_data = (double*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_BYTE: {
        signed char* mpi_data = (signed char*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_UBYTE: {
        unsigned char* mpi_data = (unsigned char*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_SHORT: {
        int16_t* mpi_data = (int16_t*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_USHORT: {
        uint16_t* mpi_data = (uint16_t*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_INT: {
        int32_t* mpi_data = (int32_t*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_UINT: {
        uint32_t* mpi_data = (uint32_t*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_LONG: {
        int64_t* mpi_data = (int64_t*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    case FQ::FQT_ULONG: {
        uint64_t* mpi_data = (uint64_t*)data;
        mpi_berr = dataFile->setArrayData
            (varInfo.getPath(), offsets, counts, strides, &mpi_data[data_pos]);
        break;}
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "FastQuery::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << "): unknown FQ data type " << fqType;
        mpi_berr = false;
    }
    MPI_Allreduce(&mpi_berr, &berr, 1, MPI_BYTE, MPI_BAND, mpi_comm);
    }
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "FastQuery::setData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << "): failed to set data to dataset \""
            << varInfo.getPath().c_str() << "\"";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "FastQuery::setData(" << varPathStr.c_str() << ", "
            << varNameStr.c_str() << "): successfully set data to dataset \""
            << varInfo.getPath().c_str() << "\"";
        return true;
    }
} // IndexBuilder::setData

bool IndexBuilder::setAttribute
(const std::string &varNameStr, const std::string &attrNameStr,
 const void *values, const uint64_t len, const FQ::DataType fqType,
 const std::string &varPathStr)
{
    if (! isValid("IndexBuilder::setAttribute")) return false;
    bool berr = true;
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::setAttribute:"
            << " no valid variable is found from the file"
            << " according to the input variable information "
            << " varNameStr \"" << varNameStr.c_str() << "\""
            << " varPathStr \"" << varPathStr.c_str() << "\"";
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
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- " << "IndexBuilder::setAttribute("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << ", "
            << attrNameStr.c_str() << "):"
            << " more than one valid variables are found from the file"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList[idx].getPath() << "\"";
    }
#ifdef FQ_NOMPI
    berr = dataFile->setAttribute(varInfoList[idx].getPath(), attrNameStr,
                                  values, len, fqType);
#else
    // only the first processor needs to set the attribute
    berr = dataFile->setAttribute(varInfoList[idx].getPath(), attrNameStr,
                                  values, len, fqType);
    MPI_Bcast(&berr, 1, MPI_BYTE, 0, mpi_comm);
#endif
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "IndexBuilder::setAttribute:"
            << " fail to write attribute \"" << attrNameStr.c_str()
            << "\" to the dataset \"" << varInfoList[idx].getPath() << "\"";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "IndexBuilder::setAttribute: successfully wrote attribute \""
            << attrNameStr.c_str() << "\" to the dataset \""
            << varInfoList[idx].getPath() << "\"";
        return true;
    }
} // IndexBuilder::setAttribute
