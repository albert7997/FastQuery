#include "queryProcessor.h"
#include <fstream>
#include <cmath>

/*!
  \brief Open a file and provide query and selected data retrieving
  functions.
  The file must be already exist, and it will be opened with read
  permission only.

  \param dataFileName IN: Name of a file containing data.

  \param filefmt IN: Enumerate data format model for the file.

  \param indexFileName IN: Name of a file containing indexes.

  \param v: verbosenes for print out message.

  \param rcfile: Name of the runtime control file for FastQuery and FastBit.

  \param logfile: Name of the log file for FastQuery and FastBit.

  \param comm IN: MPI communication channel for the file.
  It is only used for MPI parallel mode.
*/
#ifdef FQ_NOMPI
QueryProcessor::QueryProcessor(const std::string& dataFileName,
                               const FQ::FileFormat filefmt,
                               const std::string& indexFileName,
                               const int v, const char *rcfile,
                               const char *logfile, void *extra)
    : FastQuery(dataFileName, filefmt, indexFileName, v,
                rcfile, logfile, true, extra)
#else
QueryProcessor::QueryProcessor(const std::string& dataFileName,
                               const FQ::FileFormat filefmt,
                               const std::string& indexFileName,
                               const int v, const char *rcfile,
                               const char *logfile, MPI_Comm comm,
                               void *extra)
      : FastQuery(dataFileName, filefmt, indexFileName, v,
                  rcfile, logfile, true, comm, extra)
#endif
{
} // QueryProcessor::QueryProcessor

/*!
  \brief Get the values of data selected from a variable.

  \param varNameStr IN: A variable name string in the file structure.

  \param coords IN: Vector containing the coordinates of the points/boxes
  for the data values to be retrieved.
  All coordinates have the same dimensionality (rank)
  as the dataspace they are located within.
  The list of boxes is formatted as follows:
  <"start" coordinate>, immediately followed by
  <"opposite" corner coordinate>, followed by
  the next "start" and "opposite" coordinates, etc.
  until all of the blocks have been listed.

  \param data OUT: Data values. The memory size of the data must be reserved.

  \param varPathStr IN: A variable path string in the file
  structure (optional). If varPathStr is not given, it means any
  variable path.

  \param selectForm: The selection form of data selection.
  - If selectForm is POINTS_SELECTION, the data is retrieved by a vector
  of points.
  - If selectForm is LINES_SELECTION, the data is retrieved by a vector
  of bounding lines.
  - If selectForm is BOXES_SELECTION, the data is retrieved by a vector
  of bounding boxes.
  If selectForm is not specified, the default is POINTS_SELECTION.

  \return True if success, otherwise return False.
*/
bool QueryProcessor::getSelectedData
(const std::string &varNameStr,
 const std::vector<uint64_t> &coords, void *data,
 const std::string &varPathStr, const FQ::SelectForm selectForm)
{
    if (! isValid("QueryProcessor::getSelectedData")) return false;

    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList, varPathStr, varNameStr);
    if ( numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::getSelectedData("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
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
            << "Warning -- QueryProcessor::getSelectedData("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList[idx].getPath() << "\"";
    }
    VarInfo varInfo = varInfoList[idx];
    VarSpace varSpace = varSpaceList[idx];

    if (selectForm == FQ::POINTS_SELECTION) {
        unsigned int nDims = varInfo.getNDims();
        unsigned int nElements = coords.size()/nDims;
        if (varInfo.getSize() == varSpace.getSize()) {
            if (! dataFile->getPointData(varInfo.getPath(), coords, data)) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::getSelectedData("
                    << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
                    << " failed to get data";
                return false;
            }
        } else {
            // convert to absolute coords based on space info.
            std::vector<uint64_t> absCoords;
            absCoords.reserve(coords.size());
            std::vector<uint64_t> spaceOffsets = varSpace.getOffsets();
            std::vector<uint64_t> spaceCounts = varSpace.getCounts();
            std::vector<uint64_t> spaceStrides = varSpace.getStrides();
            unsigned int idx = 0;
            for (unsigned int i=0; i<nElements; i++) {
                for (unsigned int j=0; j<nDims; j++) {
                    if (coords[idx] >= spaceCounts[j]) {
                        LOGGER(ibis::gVerbose > 0)
                            << "Warning -- HDF5::getSelectedData("
                            << varPathStr.c_str() << ", " << varNameStr.c_str()
                            << "): data coordinate is out of range: "
                            << spaceCounts[j];
                        return false;
                    }
                    absCoords[idx] = coords[idx]* spaceStrides[j] +
                        spaceOffsets[j];
                    idx++;
                }
            }
            if (! dataFile->getPointData(varInfo.getPath(), absCoords, data)) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::getSelectedData("
                    << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
                    << " failed to get data";
                return false;
            }
        }
        LOGGER(ibis::gVerbose > 2)
            << "QueryProcessor::getSelectedData("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " get data in number of points: " << nElements;
    } else if (selectForm == FQ::BOXES_SELECTION) {
        unsigned int nDims = varInfo.getNDims();
        unsigned int nBlocks = coords.size()/nDims/2;
        std::vector<uint64_t> spaceOffsets = varSpace.getOffsets();
        std::vector<uint64_t> spaceCounts = varSpace.getCounts();
        std::vector<uint64_t> spaceStrides = varSpace.getStrides();
        unsigned int nElements = 0;
        for(unsigned int i=0; i<nBlocks; i++) {
            std::vector<uint64_t> offsets;
            std::vector<uint64_t> counts;
            std::vector<uint64_t> strides;
            offsets.resize(nDims);
            counts.resize(nDims);
            strides.resize(nDims);
            uint64_t idx = 1;
            for (unsigned int dim=0; dim<nDims; dim++) {
                if (coords[i*nDims*2+dim] > spaceCounts[dim] ||
                    coords[i*nDims*2+nDims+dim] > spaceCounts[dim]) {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- HDF5::getSelectedData("
                        << varPathStr.c_str() << ", " << varNameStr.c_str()
                        << "): data coordinate is out of range: "
                        << spaceCounts[dim];
                    return false;
                }
                offsets[dim] = coords[i*nDims*2+dim];
                counts[dim] = coords[i*nDims*2+nDims+1] - coords[i*nDims*2+dim];
                offsets[dim] += spaceOffsets[dim];
                strides[dim] = spaceStrides[dim];
                idx *= counts[dim];
            }

#ifdef DEBUG
            for (int i=0; i<offsets.size(); i++)
                std::cout << "offset: " << offsets[i] << std::endl;
            for (int i=0; i<counts.size(); i++)
                std::cout << "count: " << counts[i] << std::endl;
            for (int i=0; i<strides.size(); i++)
                std::cout << "stride: " << strides[i] << std::endl;
#endif
            bool berr = false;
            switch (varInfo.getType()) {
            case FQ::FQT_DOUBLE:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((double*)data)[nElements]));
                break;
            case FQ::FQT_FLOAT:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((float*)data)[nElements]));
                break;
            case FQ::FQT_BYTE:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((signed char*)data)[nElements]));
                break;
            case FQ::FQT_UBYTE:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((unsigned char*)data)[nElements]));
                break;
            case FQ::FQT_SHORT:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((int16_t*)data)[nElements]));
                break;
            case FQ::FQT_USHORT:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((uint16_t*)data)[nElements]));
                break;
            case FQ::FQT_INT:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((int32_t*)data)[nElements]));
                break;
            case FQ::FQT_UINT:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((uint32_t*)data)[nElements]));
                break;
            case FQ::FQT_LONG:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((int64_t*)data)[nElements]));
                break;
            case FQ::FQT_ULONG:
                berr = dataFile->getArrayData
                    (varInfo.getPath(), offsets, counts, strides,
                     &(((uint64_t*)data)[nElements]));
                break;
            default:
                berr = false;
            }
            if (! berr) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- HDF5::getSelectedData("
                    << varPathStr.c_str() << ", " << varNameStr.c_str()
                    << "): failed to get data";
                return false;
            }
            nElements += idx;
        }
        LOGGER(ibis::gVerbose > 2)
            << "QueryProcessor::getSelectedData("
            << varPathStr.c_str() << ", " << varNameStr.c_str()
            << "): get data in number of blocks: " << nBlocks;
    }
    return true;
} // QueryProcessor::getSelectedData

uint64_t QueryProcessor::getNumHits
(const std::string &query, const std::string &varPathStr,
 uint64_t mpi_dim, uint64_t mpi_len)
{
    if (! isValid("QueryProcessor::getNumHits")) return 0;

    ibis::horometer timer;
    timer.start();
    std::string queryFB;
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList, varSpaceList, mpi_dim) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::getNumHits("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }
#ifdef FQ_NOMPI
    FQ_Part part(varInfoList, varSpaceList, *dataFile, *indexFile);
    const char* tok = part.createQuery(queryFB);
    int64_t totalHits = part.submitQuery(tok);
#else
    int64_t numHits = 0;
    int64_t totalHits = 0;
    //std::vector<uint64_t> dims = varSpaceList[0].getCounts();
    uint64_t totalCount = varSpaceList[0].getSize();
    uint64_t mpi_max_iter = varSpaceList[0].getMpiMaxIter();
    if (mpi_len == 0)
	mpi_len = FQ_DEFAULT_MPI_LEN;
    if (mpi_len * mpi_size * mpi_max_iter < totalCount)
	mpi_max_iter = std::ceil((double)totalCount / (mpi_len*mpi_size));
    for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
        ibis::horometer splitTimer;
        splitTimer.start();
        uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
        std::vector<VarSpace> varSplitSpaceList;
        for(unsigned int i=0; i<varSpaceList.size(); i++) {
            std::string mpi_text = varSpaceList[i].getText();
            std::vector<uint64_t> splitOffsets = varSpaceList[i].getOffsets();
            std::vector<uint64_t> splitCounts = varSpaceList[i].getCounts();
            std::vector<uint64_t> strides = varSpaceList[i].getStrides();

            splitOffsets[mpi_dim] += mpi_idx*mpi_len*strides[mpi_dim];
            if (mpi_idx*mpi_len >= totalCount) {
                // out of the range
                splitCounts[mpi_dim] = 0;
            } else if ((mpi_idx+1)*mpi_len>totalCount &&
                       (totalCount%mpi_len)!= 0) {
                // the last array may not have the full length
                splitCounts[mpi_dim] = totalCount%mpi_len;
            } else {
                splitCounts[mpi_dim] = mpi_len;
            }
            VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
            splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                     mpi_iter, mpi_max_iter, mpi_comm);
            varSplitSpaceList.push_back(splitVarSpace);
        }
        splitTimer.stop();
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::getNumHits:init\t"
            << splitTimer.CPUTime() << "\t" << splitTimer.realTime()
            << "\t" << mpi_rank << "\t" << mpi_size;

        ibis::horometer initTimer, submitTimer;
        initTimer.start();
        FQ_Part part(varInfoList, varSplitSpaceList, *dataFile, *indexFile,
                     mpi_rank);
        const char* tok = part.createQuery(queryFB);
        initTimer.stop();
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::getNumHits:createQuery\t"
            << initTimer.CPUTime() << "\t" << initTimer.realTime()
            << "\t" << mpi_rank;

        submitTimer.start();
        int64_t hits = part.submitQuery(tok);
        if (hits < 0) hits = 0;
        numHits += hits;
        LOGGER(ibis::gVerbose > 2)
            << "QueryProcessor::getNumHits("
            << query << ", " << varPathStr.c_str() << "):"
            << " hit count[" << mpi_iter << ", " << mpi_idx <<"]: " << hits;
        submitTimer.stop();
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::getNumHits:evaluate\t"
            << submitTimer.CPUTime() << "\t" << submitTimer.realTime()
            << "\t" << mpi_rank << "\t" << mpi_size;
    }
#ifdef FQ_PERFORMANCE_TEST
    ibis::horometer waitTimer;
    waitTimer.start();
    MPI_Barrier(mpi_comm);
    waitTimer.stop();
#endif
    ibis::horometer mpiTimer;
    mpiTimer.start();
    MPI_Allreduce(&numHits, &totalHits, 1, MPI_LONG, MPI_SUM, mpi_comm);
    mpiTimer.stop();
    LOGGER(FastQuery::reportTiming())
#ifdef FQ_PERFORMANCE_TEST
        << "\nStatistic\tQueryProcessor::getNumHits:wait\t"
        << waitTimer.CPUTime() << "\t" << waitTimer.realTime()
        << "\t" << mpi_rank
#endif
        << "\nStatistic\tQueryProcessor::getNumHits:MPI_Allreduce()\t"
        << mpiTimer.CPUTime() << "\t" << mpiTimer.realTime() << "\t"
        << mpi_rank;

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::getNumHits:results\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << totalHits;
#endif
    LOGGER(ibis::gVerbose > 2)
        << "QueryProcessor::getNumHits("
        << query << ", " << varPathStr.c_str() << "):"
        << " total hit count: " << totalHits;
    return totalHits;
} // QueryProcessor::getNumHits

uint64_t QueryProcessor::executeQuery
(const std::string &query, std::vector<uint64_t>& coords,
 const std::string &varPathStr, const FQ::SelectForm selectForm,
 uint64_t mpi_dim, uint64_t mpi_len,
 const bool bcast)
{
    if (! isValid("QueryProcessor::executeQuery")) return 0;

    std::string queryFB;
    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    if (! metadataMgr->parseQuery(query, varPathStr, queryFB, varInfoList,
                                  varSpaceList, mpi_dim) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::executeQuery("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }
#ifdef FQ_NOMPI
    FQ_Part part(varInfoList, varSpaceList, *dataFile, *indexFile);
    const char* tok = part.createQuery(queryFB);
    int64_t totalHits = part.submitQuery(tok);
    if (selectForm==FQ::POINTS_SELECTION) {
        if (totalHits > 0) part.getHitsAsPoints(tok, coords);
    } else if (selectForm==FQ::BOXES_SELECTION) {
        if (totalHits > 0) part.getHitsAsBoxes(tok, coords);
    }
#else
    int64_t totalHits = 0;
    int64_t numHits = 0;
    std::vector<uint64_t> splitCoords;
    const std::vector<uint64_t>& dims = varSpaceList[0].getCounts();
    uint64_t totalCount = varSpaceList[0].getSize();
    uint64_t mpi_max_iter = varSpaceList[0].getMpiMaxIter();
    uint64_t len = dims.size();
    if (selectForm==FQ::BOXES_SELECTION) len *= 2;

    if (mpi_len == 0)
	mpi_len = FQ_DEFAULT_MPI_LEN;
    if (mpi_len * mpi_size * mpi_max_iter < totalCount)
	mpi_max_iter = std::ceil((double)totalCount / (mpi_len*mpi_size));
    for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
        uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
        std::vector<VarSpace> varSplitSpaceList;
        for(unsigned int i=0; i<varSpaceList.size(); i++) {
            std::string mpi_text = varSpaceList[i].getText();
            std::vector<uint64_t> splitOffsets = varSpaceList[i].getOffsets();
            std::vector<uint64_t> splitCounts = varSpaceList[i].getCounts();
            std::vector<uint64_t> strides = varSpaceList[i].getStrides();

            splitOffsets[mpi_dim] += mpi_idx*mpi_len*strides[mpi_dim];
            if (mpi_idx*mpi_len >= totalCount) {
                // out of the range
                splitCounts[mpi_dim] = 0;
            } else if ((mpi_idx+1)*mpi_len>totalCount &&
                       (totalCount%mpi_len)!= 0) {
                // the last array may not have the full length
                splitCounts[mpi_dim] = totalCount%mpi_len;
            } else {
                splitCounts[mpi_dim] = mpi_len;
            }
            VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
            splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                     mpi_iter, mpi_max_iter, mpi_comm);
            varSplitSpaceList.push_back(splitVarSpace);
        }
        FQ_Part part(varInfoList, varSplitSpaceList, *dataFile, *indexFile);
        const char* tok = part.createQuery(queryFB);
        int64_t hits = part.submitQuery(tok);
        if (hits < 0) hits = 0;
        std::vector<uint64_t> tmpCoords;
        tmpCoords.reserve(hits*len);
        if (selectForm==FQ::POINTS_SELECTION) {
            if (hits > 0) part.getHitsAsPoints(tok, tmpCoords);
            for (unsigned int i=0; i<hits; i++) {
                tmpCoords[i*len+mpi_dim] += mpi_len*mpi_idx;
            }
        } else if (selectForm==FQ::BOXES_SELECTION) {
            if (hits > 0) part.getHitsAsBoxes(tok, tmpCoords);
#ifdef DEBUG
            for (int i=0; i<tmpCoords.size(); i++)
                std::cout << "coord:" << tmpCoords[i] << std::endl;
#endif
            for (unsigned int i=0; i<tmpCoords.size()/dims.size()/2; i++) {
                tmpCoords[i*len+mpi_dim] += mpi_len*mpi_idx;
                tmpCoords[i*len+dims.size()+mpi_dim] += mpi_len*mpi_idx;
            }
#ifdef DEBUG
            for (int i=0; i<tmpCoords.size(); i++)
                std::cout << "coord:" << tmpCoords[i] << std::endl;
#endif
        }
        splitCoords.insert(splitCoords.end(), tmpCoords.begin(),
                           tmpCoords.end());
        numHits += hits;
        LOGGER(ibis::gVerbose > 2)
            << "QueryProcessor::executeQuery("
            << query.c_str() << ", " << varPathStr.c_str() << "):"
            << " hit count[" << mpi_iter << ", " << mpi_idx <<"]: "
            << hits << "/" << numHits;
    }
    MPI_Allreduce(&numHits, &totalHits, 1, MPI_LONG, MPI_SUM, mpi_comm);

    int numCoords[mpi_size];
    int coordSize = splitCoords.size();
    MPI_Allgather(&coordSize, 1, MPI_INT, numCoords, 1, MPI_INT, mpi_comm);
    int displs[mpi_size];
    displs[0] = 0;
    for (int i=1; i<mpi_size; i++) {
        displs[i] = displs[i-1] + numCoords[i-1];
    }
    int totCoords = displs[mpi_size-1] + numCoords[mpi_size-1];
    if (coords.capacity() < totCoords) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::executeQuery("
            << query << ", " << varPathStr.c_str() << "):"
            << " buffer for retrieving coordinates is too small "
            << coords.capacity() << " < " << totCoords;
        return 0;
    }
    if (bcast) {
	coords.resize(totCoords);
	MPI_Allgatherv(&splitCoords[0], splitCoords.size(), MPI_UNSIGNED_LONG_LONG,
                   &coords[0], numCoords, displs, MPI_UNSIGNED_LONG_LONG, mpi_comm);
    } else {
        if (mpi_rank==0) coords.resize(totCoords);
        MPI_Gatherv(&splitCoords[0], splitCoords.size(), MPI_UNSIGNED_LONG_LONG,
                   &coords[0], numCoords, displs, MPI_UNSIGNED_LONG_LONG, 0, mpi_comm);
    }
#endif
    LOGGER(ibis::gVerbose > 2)
        << "QueryProcessor::executeQuery("
        << query << ", " << varPathStr.c_str() << "):"
        << " total hit count: " << totalHits;
    return totalHits;
} // QueryProcessor::executeQuery

uint64_t QueryProcessor::executeEqualitySelectionQuery
(const std::string &varNameStr,
 const std::vector<double>& identifiers, std::vector<uint64_t>& coords,
 const std::string &varPathStr, const FQ::SelectForm selectForm,
 uint64_t mpi_dim, uint64_t mpi_len, const bool bcast)
{
    if (! isValid("QueryProcessor::executeEqualitySelectionQuery")) return 0;

    std::vector<VarInfo> varInfoList;
    std::vector<VarSpace> varSpaceList;
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList, varSpaceList,varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::executeEqualitySelectionQuery("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
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
            << "Warning -- QueryProcessor::executeEqualitySelectionQuery("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList[idx].getPath() << "\"";
    }

    std::vector<VarInfo> selectInfoList;
    std::vector<VarSpace> selectSpaceList;
    selectInfoList.push_back(varInfoList[idx]);
    selectSpaceList.push_back(varSpaceList[idx]);
#ifdef FQ_NOMPI
    FQ_Part part(selectInfoList, selectSpaceList, *dataFile, *indexFile);
    // the name of first variable is always named to Var0
    const char* tok = part.createEqualitySelectionQuery("Var0", identifiers);
    int64_t totalHits = part.submitQuery(tok);
    if (selectForm==FQ::POINTS_SELECTION) {
        if (totalHits > 0) part.getHitsAsPoints(tok, coords);
    } else if (selectForm==FQ::BOXES_SELECTION) {
        if (totalHits > 0) part.getHitsAsBoxes(tok, coords);
    }
#else
    int64_t totalHits = 0;
    int64_t numHits = 0;
    std::vector<uint64_t> splitCoords;
    std::vector<uint64_t> offsets = selectSpaceList[0].getOffsets();
    std::vector<uint64_t> counts = selectSpaceList[0].getCounts();
    std::vector<uint64_t> strides = selectSpaceList[0].getStrides();
    unsigned int totalCount = selectSpaceList[0].getSize();
    uint64_t mpi_max_iter = selectSpaceList[0].getMpiMaxIter();
    uint64_t len = offsets.size();
    if (selectForm==FQ::BOXES_SELECTION) len *= 2;
    std::string mpi_text = selectSpaceList[0].getText();

    if (mpi_len == 0)
	mpi_len = FQ_DEFAULT_MPI_LEN;
    if (mpi_len * mpi_size * mpi_max_iter < totalCount)
	mpi_max_iter = std::ceil((double)totalCount / (mpi_len*mpi_size));
    for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
        uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
        std::vector<VarSpace> varSplitSpaceList;
        std::vector<uint64_t> splitOffsets = offsets;
        std::vector<uint64_t> splitCounts = counts;

        splitOffsets[mpi_dim] += mpi_idx*mpi_len*strides[mpi_dim];
        if (mpi_idx*mpi_len >= totalCount) {
            // out of the range
            splitCounts[mpi_dim] = 0;
        } else if ((mpi_idx+1)*mpi_len>totalCount && (totalCount%mpi_len)!= 0) {
            // the last array may not have the full length
            splitCounts[mpi_dim] = totalCount%mpi_len;
        } else {
            splitCounts[mpi_dim] = mpi_len;
        }
        VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
        splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                 mpi_iter, mpi_max_iter, mpi_comm);
        varSplitSpaceList.push_back(splitVarSpace);
        FQ_Part part(selectInfoList, varSplitSpaceList, *dataFile, *indexFile);
        const char* tok =
            part.createEqualitySelectionQuery("Var0", identifiers);
        int64_t hits = part.submitQuery(tok);
        if (hits < 0) hits = 0;
        std::vector<uint64_t> tmpCoords;
        tmpCoords.reserve(hits*len);
        if (selectForm==FQ::POINTS_SELECTION) {
            if (hits > 0) part.getHitsAsPoints(tok, tmpCoords);
        } else if (selectForm==FQ::BOXES_SELECTION) {
            if (hits > 0) part.getHitsAsBoxes(tok, tmpCoords);
        }
        for (unsigned int i=0; i<tmpCoords.size(); i++) {
            tmpCoords[i*len+mpi_dim] += mpi_len*mpi_idx;
        }
        splitCoords.insert(splitCoords.end(), tmpCoords.begin(),
                           tmpCoords.end());
        numHits += hits;
        LOGGER(ibis::gVerbose > 2)
            << "QueryProcessor::executeEqualitySelectionQuery("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " hit count[" << mpi_iter << ", " << mpi_idx <<"]: " << hits;
    }
    MPI_Allreduce(&numHits, &totalHits, 1, MPI_LONG, MPI_SUM, mpi_comm);

    int numCoords[mpi_size];
    int coordSize = splitCoords.size();
    MPI_Allgather(&coordSize, 1, MPI_INT, numCoords, 1, MPI_INT, mpi_comm);
    int displs[mpi_size];
    displs[0] = 0;
    for (int i=1; i<mpi_size; i++) {
        displs[i] = displs[i-1] + numCoords[i-1];
    }
    int totCoords = displs[mpi_size-1] + numCoords[mpi_size-1];
    if (coords.capacity() < totCoords) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::executeEqualitySelectionQuery("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " buffer for retrieving coordinates is too small "
            << coords.capacity() << " < " << totCoords;
        return 0;
    }
    if (bcast) {
        coords.resize(totCoords);
        MPI_Allgatherv(&splitCoords[0], splitCoords.size(), MPI_UNSIGNED_LONG_LONG,
                   &coords[0], numCoords, displs, MPI_UNSIGNED_LONG_LONG, mpi_comm);
   } else {
        if (mpi_rank==0) coords.resize(totCoords);
        MPI_Gatherv(&splitCoords[0], splitCoords.size(), MPI_UNSIGNED_LONG_LONG,
                   &coords[0], numCoords, displs, MPI_UNSIGNED_LONG_LONG, 0, mpi_comm);
   }
#endif
    LOGGER(ibis::gVerbose > 2)
        << "QueryProcessor::executeEqualitySelectionQuery("
        << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
        << " total hit count: " << totalHits;
    return totalHits;
} // QueryProcessor::executeEqualitySelectionQuery

//Added by Suren: Histogram functions

/*!
  \brief Get 1D histogram
  Gets the number of records falling in the regular bins defined by the
  <tt>begin:end:stride</tt> triplet.  The triplet defines
  <tt> 1 + floor((end-begin)/stride) </tt> bins:
  @code
  [begin, begin+stride)
  [begin+stride, begin+stride*2)
  ...
  [begin+stride*floor((end-begin)/stride), end].
  @endcode
  Note that the bins all have closed ends on the left, and open ends on
  the right, except the last bin where both ends are closed.

  \param varPathStr IN: A variable name string in the file structure

  \param query IN: SQL-like query conditions

  \param begin IN: begining of the bins

  \param end IN: end of the bins

  \param stride IN: stride

  \param counts OUT: output counts
*/
uint64_t QueryProcessor::get1DHistogram
(const char *query,
 const std::string& varNameStr,
 const std::string& varPathStr,
 double begin, double end, double stride,
 std::vector<uint32_t> &counts,
 uint64_t mpi_dim, uint64_t mpi_len)
{
    if (! isValid("QueryProcessor::get1DHistogram")) {
        return 0;
    }

    ibis::horometer timer;
    timer.start();
    //Converted variable name to fit FastBit names
    std::string varFBName = "Var";
    std::vector<VarInfo> varInfoList_var;
    std::vector<VarSpace> varSpaceList_var;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery(query, varPathStr, queryFB,
                                  varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList_var, varSpaceList_var, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var[i].getPath().size()) {
                idx = i;
                len = varInfoList_var[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var[idx]);
        varSpaceList_query.push_back(varSpaceList_var[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str[25];
    sprintf(str, "%u", varIdx);
    varFBName += str;
    unsigned int dims1 = static_cast<uint32_t>(1+floor((end-begin)/stride));
    counts.assign(dims1,0);

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile,
                 *indexFile);

    part.get1DHist_counts (queryFB, varFBName, begin, end, stride, counts);
    LOGGER(ibis::gVerbose > 0)
        << "Testing -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << "test the serial version";

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get1DHistoram:timer\t"
        << timer.CPUTime() << "\t" << timer.realTime()
        << "\t";
#else

    std::vector<uint64_t> dims = varSpaceList_query[0].getCounts();
    uint64_t totalCount = varSpaceList_query[0].getSize();
    uint64_t mpi_max_iter = varSpaceList_query[0].getMpiMaxIter();
    std::vector<uint64_t> tempcounts; // for collecting histogram
    tempcounts.assign(dims1,0);
    if (mpi_len == 0)
	mpi_len = FQ_DEFAULT_MPI_LEN;
    if (mpi_len * mpi_size * mpi_max_iter < totalCount)
	mpi_max_iter = std::ceil((double)totalCount / (mpi_len*mpi_size));
    for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
        ibis::horometer splitTimer;
        splitTimer.start();
        uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
        std::vector<VarSpace> varSplitSpaceList;

        for(unsigned int i=0; i<varSpaceList_query.size(); i++) {
            std::string mpi_text = varSpaceList_query[i].getText();
            std::vector<uint64_t> splitOffsets =
                varSpaceList_query[i].getOffsets();
            std::vector<uint64_t> splitCounts =
                varSpaceList_query[i].getCounts();
            std::vector<uint64_t> strides =
                varSpaceList_query[i].getStrides();

            splitOffsets[mpi_dim] += mpi_idx*mpi_len*strides[mpi_dim];
            if (mpi_idx*mpi_len >= totalCount) {
                // out of the range
                splitCounts[mpi_dim] = 0;
            } else if ((mpi_idx+1)*mpi_len>totalCount &&
                       (totalCount%mpi_len)!= 0) {
                // the last array may not have the full length
                splitCounts[mpi_dim] = totalCount%mpi_len;
            } else {
                splitCounts[mpi_dim] = mpi_len;
            }
            VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
            splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                     mpi_iter, mpi_max_iter, mpi_comm);
            varSplitSpaceList.push_back(splitVarSpace);
        }
        splitTimer.stop();
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::get1DHistoram:init\t"
            << splitTimer.CPUTime() << "\t" << splitTimer.realTime()
            << "\t" << mpi_rank << "\t" << mpi_size;
        ibis::horometer HistTimer;
        std::vector<uint32_t> itercounts; // for collecting iteration
        HistTimer.start();
        FQ_Part part(varInfoList_query, varSplitSpaceList, *dataFile,
                     *indexFile, mpi_rank);
        part.get1DHist_counts (queryFB, varFBName, begin, end, stride,
                               itercounts);
        HistTimer.stop();
        for (int i=0;i<itercounts.size();i++)
            tempcounts[i] += itercounts[i];
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::get1DHistoram:Historam\t"
            << HistTimer.CPUTime() << "\t" << HistTimer.realTime()
            << "\t" << mpi_rank;
    }
    ibis::horometer reduceTimer;
    reduceTimer.start();
    void* pt = static_cast<void*>(&counts[0]);
    MPI_Allreduce(&tempcounts[0], pt, (int)counts.size(), MPI_UNSIGNED,
                  MPI_SUM, mpi_comm);
    reduceTimer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get1DHistoram:Allreduce\t"
        << reduceTimer.CPUTime() << "\t" << reduceTimer.realTime()
        << "\t" << mpi_rank;
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get1DHistoram:timer\t"
        << timer.CPUTime() << "\t" << timer.realTime()
        << "\t" << mpi_rank;
#endif

    return 1;
}       // End QueryProcessor::get1DHistogram counts

// Compute weighted conditional 1D histogram with regularly spaced bins.
uint64_t QueryProcessor::get1DHistogram
(const char *query,
 const std::string& varNameStr,
 const std::string& varPathStr,
 double begin, double end, double stride,
 const std::string& wtNameStr,
 std::vector<double> &weights) {
    if (! isValid("QueryProcessor::get1DHistogram")) {
        return 0;
    }
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;
    //Converted variable name to fit FastBit names
    std::string varFBName = "Var";
    std::vector<VarInfo> varInfoList_var;
    std::vector<VarSpace> varSpaceList_var;
    //Converted weight name to fit FastBit names
    std::string weightFBName = "Var";
    std::vector<VarInfo> varInfoList_weight;
    std::vector<VarSpace> varSpaceList_weight;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList_var, varSpaceList_var, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var[i].getPath().size()) {
                idx = i;
                len = varInfoList_var[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var[idx]);
        varSpaceList_query.push_back(varSpaceList_var[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str[25];
    sprintf(str, "%u", varIdx);
    varFBName += str;

    // verify if weight name exists and
    // pick the best fit weight name if there are multiple
    unsigned int numWts = metadataMgr->getAllVariables
        (varInfoList_weight, varSpaceList_weight, varPathStr, wtNameStr);
    if (numWts == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " no valid weight is found";
        return 0;
    }

    idx = 0;
    if (numWts > 1) {
        // pick the best fit weight with the shortest length
        std::string wtName = varInfoList_weight[0].getPath();
        int wlen = wtName.size();
        for (unsigned int i = 1; i < numWts; i++) {
            if (wlen > varInfoList_weight[i].getPath().size()) {
                idx = i;
                wlen = varInfoList_weight[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " multiple weight variables are found. Use variable \""
            << varInfoList_weight[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int wtIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_weight[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_weight[idx])) {
                wtIdx = (int) i;
            }
        }
    }

    // if weight name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (wtIdx == -1) {
        wtIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_weight[idx]);
        varSpaceList_query.push_back(varSpaceList_weight[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", wtIdx);
    weightFBName += str2;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile,
                 *indexFile);

    part.get1DHist_weights (queryFB, varFBName, begin, end, stride,
                            weightFBName, weights);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 0;
}

// Compute 1D histogram with adaptive bins.
uint64_t QueryProcessor::get1DHistogram
(const std::string& varNameStr,
 uint32_t nbin,
 const std::string& varPathStr,
 std::vector<double> &bounds,
 std::vector<uint32_t> &counts) {
    return 1;
}


// Compute conditional 1D histogram with adaptive bins.
uint64_t QueryProcessor::get1DHistogram
(const char *query,
 const std::string& varNameStr,
 uint32_t nbin,
 const std::string& varPathStr,
 std::vector<double> &bounds,
 std::vector<uint32_t> &counts) {
    if (! isValid("QueryProcessor::get1DHistogram")) {
        return 0;
    }
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;
    //Converted variable name to fit FastBit names
    std::string varFBName = "Var";
    std::vector<VarInfo> varInfoList_var;
    std::vector<VarSpace> varSpaceList_var;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query))  {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList_var, varSpaceList_var, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var[i].getPath().size()) {
                idx = i;
                len = varInfoList_var[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var[idx]);
        varSpaceList_query.push_back(varSpaceList_var[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str[25];
    sprintf(str, "%u", varIdx);
    varFBName += str;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile,
                 *indexFile);

    part.get1DAdaptiveHist_cond_counts(queryFB, varFBName, nbin,
                                       bounds, counts);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif

    return 1;
}

// Partition values of the named variable into regularly spaced bins.
uint64_t QueryProcessor::get1DBins
(const char *query,
 const std::string& varNameStr,
 const std::string& varPathStr,
 double begin, double end, double stride,
 std::vector<ibis::bitvector> &bins) {
    if (! isValid("QueryProcessor::get1DHistogram")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName = "Var";
    std::vector<VarInfo> varInfoList_var;
    std::vector<VarSpace> varSpaceList_var;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query,
         varSpaceList_query)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList_var, varSpaceList_var, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var[i].getPath().size()) {
                idx = i;
                len = varInfoList_var[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var[idx]);
        varSpaceList_query.push_back(varSpaceList_var[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str[25];
    sprintf(str, "%u", varIdx);
    varFBName += str;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile,
                 *indexFile);

    part.get1DBins (queryFB, varFBName, begin, end, stride, bins);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
}       // End QueryProcessor::get1DBins

// Partition values of the named variable into regularly spaced bins.
// Passing <ibis::bitvector *> &bins
uint64_t QueryProcessor::get1DBins
(const char *query,
 const std::string& varNameStr,
 const std::string& varPathStr,
 double begin, double end, double stride,
 std::vector<ibis::bitvector *> &bins) {
    if (! isValid("QueryProcessor::get1DHistogram")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName = "Var";
    std::vector<VarInfo> varInfoList_var;
    std::vector<VarSpace> varSpaceList_var;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList_var, varSpaceList_var, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var[i].getPath().size()) {
                idx = i;
                len = varInfoList_var[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var[idx]);
        varSpaceList_query.push_back(varSpaceList_var[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str[25];
    sprintf(str, "%u", varIdx);
    varFBName += str;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile,
                 *indexFile);

    part.get1DBins (queryFB, varFBName, begin, end, stride, bins);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
}       // End QueryProcessor::get1DBins

//Partition values of the named variable into regularly spaced bins.
//This version returns a vector of pointers to bit vectors.
// Weights
uint64_t QueryProcessor::get1DBins
(const char *query,
 const std::string& varNameStr,
 const std::string& varPathStr,
 double begin, double end, double stride,
 const std::string& wtNameStr,
 std::vector<double> &weights,
 std::vector<ibis::bitvector *> &bins) {
    if (! isValid("QueryProcessor::get1DBins")) {
        return 0;
    }

    //Converted variable name to fit FastBit names
    std::string varFBName = "Var";
    std::vector<VarInfo> varInfoList_var;
    std::vector<VarSpace> varSpaceList_var;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;
    //Converted weight name to fit FastBit names
    std::string weightFBName = "Var";
    std::vector<VarInfo> varInfoList_weight;
    std::vector<VarSpace> varSpaceList_weight;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = metadataMgr->getAllVariables
        (varInfoList_var, varSpaceList_var, varPathStr, varNameStr);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    int idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var[i].getPath().size()) {
                idx = i;
                len = varInfoList_var[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var[idx]);
        varSpaceList_query.push_back(varSpaceList_var[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str[25];
    sprintf(str, "%u", varIdx);
    varFBName += str;

    // verify if weight name exists and
    // pick the best fit weight name if there are multiple
    unsigned int numWts = metadataMgr->getAllVariables
        (varInfoList_weight, varSpaceList_weight, varPathStr, wtNameStr);
    if (numWts == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " no valid weight is found";
        return 0;
    }

    idx = 0;
    if (numWts > 1) {
        // pick the best fit weight with the shortest length
        std::string wtName = varInfoList_weight[0].getPath();
        int wlen = wtName.size();
        for (unsigned int i = 1; i < numWts; i++)
            {
                if (wlen > varInfoList_weight[i].getPath().size())
                    {
                        idx = i;
                        wlen = varInfoList_weight[i].getPath().size();
                    }
            }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " multiple weight variables are found. use the best fit variable \""
            << varInfoList_weight[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int wtIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_weight[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_weight[idx])) {
                wtIdx = (int) i;
            }
        }
    }

    // if weight name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (wtIdx == -1) {
        wtIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_weight[idx]);
        varSpaceList_query.push_back(varSpaceList_weight[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", wtIdx);
    weightFBName += str2;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get1DBins_weights (queryFB, varFBName, begin, end, stride,
                            weightFBName, weights, bins);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
}       // End QueryProcessor::get1DBins

uint64_t QueryProcessor::get1DBins
(const char *query,
 const std::string& varNameStr,
 const std::string& varPathStr,
 uint32_t nbin,
 double begin, double end, double stride,
 std::vector<double> &bounds,
 std::vector<ibis::bitvector> &bins) {
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " Function is not implemented";
    return 1;
}       // End QueryProcessor::get1DBins () : Adaptive bins
// End of 1D Histogram functions


// 2D Histogram functions
// Compute conditional 2D Histogram
uint64_t QueryProcessor::get2DHistogram
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 std::vector<uint32_t> &counts,
 uint64_t mpi_dim, uint64_t mpi_len) {
    if (! isValid("QueryProcessor::get2DHistogram")) {
        return 0;
    }

    ibis::horometer timer;
    timer.start();
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;
    unsigned int dims1 = static_cast<uint32_t>(1+floor((end1-begin1)/stride1));
    unsigned int dims2 = static_cast<uint32_t>(1+floor((end2-begin2)/stride2));
    unsigned int countsSize = dims1*dims2;
    counts.assign(countsSize,0);

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);
    part.get2DHist_counts(queryFB, varFBName1, begin1, end1, stride1,
                          varFBName2, begin2, end2, stride2, counts);
#else
    std::vector<uint64_t> dims = varSpaceList_query[0].getCounts();
    uint64_t totalCount = varSpaceList_query[0].getSize();
    uint64_t mpi_max_iter = varSpaceList_query[0].getMpiMaxIter();
    std::vector<uint32_t> tempcounts; // for collecting histogram
    tempcounts.assign(countsSize, 0);
    if (mpi_len == 0)
	mpi_len = FQ_DEFAULT_MPI_LEN;
    if (mpi_len * mpi_size * mpi_max_iter < totalCount)
	mpi_max_iter = std::ceil((double)totalCount / (mpi_len*mpi_size));
    for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
        ibis::horometer splitTimer;
        splitTimer.start();
        uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
        std::vector<VarSpace> varSplitSpaceList;

        for(unsigned int i=0; i<varSpaceList_query.size(); i++) {
            std::string mpi_text = varSpaceList_query[i].getText();
            std::vector<uint64_t> splitOffsets =
                varSpaceList_query[i].getOffsets();
            std::vector<uint64_t> splitCounts =
                varSpaceList_query[i].getCounts();
            std::vector<uint64_t> strides =
                varSpaceList_query[i].getStrides();

            splitOffsets[mpi_dim] += mpi_idx*mpi_len*strides[mpi_dim];
            if (mpi_idx*mpi_len >= totalCount) {
                // out of the range
                splitCounts[mpi_dim] = 0;
            } else if ((mpi_idx+1)*mpi_len>totalCount &&
                       (totalCount%mpi_len)!= 0) {
                // the last array may not have the full length
                splitCounts[mpi_dim] = totalCount%mpi_len;
            } else {
                splitCounts[mpi_dim] = mpi_len;
            }
            VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
            splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                     mpi_iter, mpi_max_iter, mpi_comm);
            varSplitSpaceList.push_back(splitVarSpace);
        }
        splitTimer.stop();
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::get2DHistoram:init\t"
            << splitTimer.CPUTime() << "\t" << splitTimer.realTime()
            << "\t" << mpi_rank << "\t" << mpi_size;
        ibis::horometer HistTimer;
        // collecting histogram each iteration, add them into tempcounts
        std::vector<uint32_t> itercounts;
        HistTimer.start();
        FQ_Part part(varInfoList_query, varSplitSpaceList, *dataFile,
                     *indexFile, mpi_rank);
        part.get2DHist_counts (queryFB, varFBName1, begin1, end1, stride1,
                               varFBName2, begin2, end2, stride2, itercounts);
        HistTimer.stop();
        for (int i=0;i<itercounts.size();i++) {
            tempcounts[i] += itercounts[i];
        }
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::get2DHistoram:Historam\t"
            << HistTimer.CPUTime() << "\t" << HistTimer.realTime()
            << "\t" << mpi_rank;
    }
    MPI_Barrier(mpi_comm);
    ibis::horometer reduceTimer;
    reduceTimer.start();
    void* pt = static_cast<void*>(&counts[0]);
    MPI_Allreduce(&tempcounts[0], pt, (int)counts.size(), MPI_UNSIGNED,
                  MPI_SUM, mpi_comm);
    reduceTimer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get2DHistoram:Allreduce\t"
        << reduceTimer.CPUTime() << "\t" << reduceTimer.realTime()
        << "\t" << mpi_rank;
    timer.stop();

    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get2DHistoram:timer\t"
        << timer.CPUTime() << "\t" << timer.realTime()
        << "\t" << mpi_rank;
#endif
    return 1;
}

// Compute weights for conditional 2D histogram
uint64_t QueryProcessor::get2DHistogram
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 const std::string& wtNameStr,
 std::vector<double> &weights) {
    if (! isValid("QueryProcessor::get2DHistogram")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted weight name to fit FastBit names
    std::string weightFBName = "Var";
    std::vector<VarInfo> varInfoList_weight;
    std::vector<VarSpace> varSpaceList_weight;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

    // verify if weight name exists and
    // pick the best fit weight name if there are multiple
    unsigned int numWts = metadataMgr->getAllVariables
        (varInfoList_weight, varSpaceList_weight, varPathStr, wtNameStr);
    if (numWts == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " no valid weight is found";
        return 0;
    }

    idx = 0;
    if (numWts > 1) {
        // pick the best fit weight with the shortest length
        std::string wtName = varInfoList_weight[0].getPath();
        int wlen = wtName.size();
        for (unsigned int i = 1; i < numWts; i++) {
            if (wlen > varInfoList_weight[i].getPath().size()) {
                idx = i;
                wlen = varInfoList_weight[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " multiple weight variables are found.  Use the variable \""
            << varInfoList_weight[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int wtIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_weight[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_weight[idx])) {
                wtIdx = (int) i;
            }
        }
    }

    // if weight name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (wtIdx == -1) {
        wtIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_weight[idx]);
        varSpaceList_query.push_back(varSpaceList_weight[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str3[25];
    sprintf(str3, "%u", wtIdx);
    weightFBName += str3;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get2DHist_weights (queryFB, varFBName1, begin1, end1, stride1,
                            varFBName2, begin2, end2, stride2,
                            weightFBName, weights);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
}

// Compute 2D histogram with adaptive bins.
uint64_t QueryProcessor::get2DHistogram
(const std::string& varNameStr1,
 const std::string& varNameStr2,
 const std::string& varPathStr,
 uint32_t num_bins1,
 uint32_t num_bins2,
 std::vector<double> &bounds1,
 std::vector<double> &bounds2,
 std::vector<uint32_t> &counts,
 const char* const option) {
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get2DHistogram ("
        << varPathStr.c_str() << "):"
        << " Function is not implemented";
    return 1;
}  // QueryProcessor::get2DHistogram: adaptive bins

// Compute conditional 2D histogram with adaptive bins.
uint64_t QueryProcessor::get2DHistogram
(const char *query,
 const std::string& varNameStr1,
 const std::string& varNameStr2,
 const std::string& varPathStr,
 uint32_t num_bins1,
 uint32_t num_bins2,
 std::vector<double> &bounds1,
 std::vector<double> &bounds2,
 std::vector<uint32_t> &counts) {
    if (! isValid("QueryProcessor::get2DHistogram")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get2DAdaptiveHist_cond_counts
        (queryFB, varFBName1, varFBName2, num_bins1, num_bins2,
         bounds1, bounds2, counts);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
}  // QueryProcessor::get2DHistogram: conditional adaptive bins

// 2D Bins
uint64_t QueryProcessor::get2DBins
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 std::vector<ibis::bitvector> &bins) {
    if (! isValid("QueryProcessor::get2DBins")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get2DBins (queryFB, varFBName1, begin1, end1, stride1,
                    varFBName2, begin2, end2, stride2, bins);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
} // QueryProcessor::get2DBins

// 2D Bins --> <ibis::bitvector *> &bins
uint64_t QueryProcessor::get2DBins
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 std::vector<ibis::bitvector *> &bins) {
    if (! isValid("QueryProcessor::get2DBins")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get2DBins(queryFB, varFBName1, begin1, end1, stride1,
                   varFBName2, begin2, end2, stride2, bins);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
} // QueryProcessor::get2DBins

// get2DBins_weights
uint64_t QueryProcessor::get2DBins
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 const std::string& wtNameStr,
 std::vector<double> &weights,
 std::vector<ibis::bitvector *> &bins) {
    if (! isValid("QueryProcessor::get2DBins")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;

    //Converted weight name to fit FastBit names
    std::string weightFBName = "Var";
    std::vector<VarInfo> varInfoList_weight;
    std::vector<VarSpace> varSpaceList_weight;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DBins ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

    // verify if weight name exists and
    // pick the best fit weight name if there are multiple
    unsigned int numWts = metadataMgr->getAllVariables
        (varInfoList_weight, varSpaceList_weight, varPathStr, wtNameStr);
    if (numWts == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " no valid weight is found";
        return 0;
    }

    idx = 0;
    if (numWts > 1) {
        // pick the best fit weight with the shortest length
        std::string wtName = varInfoList_weight[0].getPath();
        int wlen = wtName.size();
        for (unsigned int i = 1; i < numWts; i++) {
            if (wlen > varInfoList_weight[i].getPath().size()) {
                idx = i;
                wlen = varInfoList_weight[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " multiple weight variables are found.  Use the variable \""
            << varInfoList_weight[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int wtIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_weight[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_weight[idx])) {
                wtIdx = (int) i;
            }
        }
    }

    // if weight name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (wtIdx == -1) {
        wtIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_weight[idx]);
        varSpaceList_query.push_back(varSpaceList_weight[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str3[25];
    sprintf(str3, "%u", wtIdx);
    weightFBName += str3;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get2DBins_weights (queryFB, varFBName1, begin1, end1, stride1,
                            varFBName2, begin2, end2, stride2,
                            weightFBName, weights, bins);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get1DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
} // QueryProcessor::get2DBins ()



// 3D Histogram functions
uint64_t QueryProcessor::get3DHistogram
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 const std::string& varNameStr3,
 double begin3, double end3, double stride3,
 std::vector<uint32_t> &counts,
 uint64_t mpi_dim, uint64_t mpi_len) {

    ibis::horometer timer;
    timer.start();
    //LOGGER(FastQuery::reportTiming()) << "In 3DHistoram..." << std::endl;
    if (! isValid("QueryProcessor::get3DHistogram")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted variable name to fit FastBit names
    std::string varFBName3 = "Var";
    std::vector<VarInfo> varInfoList_var3;
    std::vector<VarSpace> varSpaceList_var3;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    // Process the second variable
    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

    // Process the third variable
    numVar = metadataMgr->getAllVariables
        (varInfoList_var3, varSpaceList_var3, varPathStr, varNameStr3);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr3.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var3[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var3[i].getPath().size()) {
                idx = i;
                len = varInfoList_var3[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr3.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var3[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var3[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var3[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var3[idx]);
        varSpaceList_query.push_back(varSpaceList_var3[idx]);

    }
    // Convert the name of variable to fit FastBit convention
    char str3[25];
    sprintf(str3, "%u", varIdx);
    varFBName3 += str3;
    unsigned int dims1 = static_cast<uint32_t>(1+floor((end1-begin1)/stride1));
    unsigned int dims2 = static_cast<uint32_t>(1+floor((end2-begin2)/stride2));
    unsigned int dims3 = static_cast<uint32_t>(1+floor((end3-begin3)/stride3));
    unsigned int countsSize = dims1*dims2*dims3;
    counts.assign(countsSize, 0);
#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile, *indexFile);

    part.get3DHist_counts  (queryFB,
                            varFBName1, begin1, end1, stride1,
                            varFBName2, begin2, end2, stride2,
                            varFBName3, begin3, end3, stride3,
                            counts);
    timer.stop();
#else
    std::vector<uint64_t> dims = varSpaceList_query[0].getCounts();
    uint64_t totalCount = varSpaceList_query[0].getSize();
    uint64_t mpi_max_iter = varSpaceList_query[0].getMpiMaxIter();
    std::vector<uint32_t> tempcounts; // for collecting histogram
    tempcounts.assign(countsSize,0);
    if (mpi_len == 0)
	mpi_len = FQ_DEFAULT_MPI_LEN;
    if (mpi_len * mpi_size * mpi_max_iter < totalCount)
	mpi_max_iter = std::ceil((double)totalCount / (mpi_len*mpi_size));
    for (uint64_t mpi_iter=0; mpi_iter<mpi_max_iter; mpi_iter++) {
        ibis::horometer splitTimer;
        splitTimer.start();
        uint64_t mpi_idx = mpi_rank + mpi_iter*mpi_size;
        std::vector<VarSpace> varSplitSpaceList;
        for(unsigned int i=0; i<varSpaceList_query.size(); i++) {
            std::string mpi_text = varSpaceList_query[i].getText();
            std::vector<uint64_t> splitOffsets =
                varSpaceList_query[i].getOffsets();
            std::vector<uint64_t> splitCounts =
		varSpaceList_query[i].getCounts();
            std::vector<uint64_t> strides = varSpaceList_query[i].getStrides();
            splitOffsets[mpi_dim] += mpi_idx*mpi_len*strides[mpi_dim];
            if (mpi_idx*mpi_len >= totalCount) {// out of the range
                splitCounts[mpi_dim] = 0;
            } else if ((mpi_idx+1)*mpi_len>totalCount &&
                       (totalCount%mpi_len)!= 0) {
                // the last array may not have the full length
                splitCounts[mpi_dim] = totalCount%mpi_len;
            } else {
                splitCounts[mpi_dim] = mpi_len;
            }
            VarSpace splitVarSpace(splitOffsets, splitCounts, strides);
            splitVarSpace.setMpiInfo(mpi_text, mpi_dim, mpi_len, mpi_idx,
                                     mpi_iter, mpi_max_iter, mpi_comm);
            varSplitSpaceList.push_back(splitVarSpace);
        }
        splitTimer.stop();
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::get3DHistoram:init\t"
            << splitTimer.CPUTime() << "\t" << splitTimer.realTime()
            << "\t" << mpi_rank << "\t" << mpi_size;
        ibis::horometer HistTimer;
        std::vector<uint32_t> itercounts;
        HistTimer.start();
        FQ_Part part(varInfoList_query, varSplitSpaceList, *dataFile,
                     *indexFile, mpi_rank);
        part.get3DHist_counts (queryFB, varFBName1, begin1, end1, stride1,
                               varFBName2, begin2, end2, stride2,
                               varFBName3, begin3, end3, stride3, itercounts);
        HistTimer.stop();
        for (int i=0;i<itercounts.size();i++)
            tempcounts[i] += itercounts[i];
        LOGGER(FastQuery::reportTiming())
            << "Statistic\tQueryProcessor::get3DHistoram:Historam\t"
            << HistTimer.CPUTime() << "\t" << HistTimer.realTime()
            << "\t" << mpi_rank;
    }
    ibis::horometer reduceTimer;
    //counts.assign(tempcounts.size(), 0);
    void* pt = static_cast<void*>(&counts[0]);
    MPI_Allreduce(&tempcounts[0], pt, (int)counts.size(), MPI_UNSIGNED_LONG, MPI_SUM, mpi_comm);
    reduceTimer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get3DHistoram:counts.size()\t"
        << counts.size()
        << "\t" << mpi_rank;
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tQueryProcessor::get3DHistoram:Allreduce\t"
        << reduceTimer.CPUTime() << "\t" << reduceTimer.realTime()
        << "\t" << mpi_rank;
    timer.stop();

#endif
    return 1;
} //QueryProcessor::get3DHistogram: conditional histogram

uint64_t QueryProcessor::get3DHistogram
(const char *query,
 const std::string& varPathStr,
 const std::string& varNameStr1,
 double begin1, double end1, double stride1,
 const std::string& varNameStr2,
 double begin2, double end2, double stride2,
 const std::string& varNameStr3,
 double begin3, double end3, double stride3,
 const std::string& wtNameStr,
 std::vector<double> &weights) {
    if (! isValid("QueryProcessor::get3DHistogram")) {
        return 0;
    }
    //Converted variable name to fit FastBit names
    std::string varFBName1 = "Var";
    std::vector<VarInfo> varInfoList_var1;
    std::vector<VarSpace> varSpaceList_var1;
    //Converted variable name to fit FastBit names
    std::string varFBName2 = "Var";
    std::vector<VarInfo> varInfoList_var2;
    std::vector<VarSpace> varSpaceList_var2;
    //Converted variable name to fit FastBit names
    std::string varFBName3 = "Var";
    std::vector<VarInfo> varInfoList_var3;
    std::vector<VarSpace> varSpaceList_var3;
    //Converted weight name to fit FastBit names
    std::string weightFBName = "Var";
    std::vector<VarInfo> varInfoList_weight;
    std::vector<VarSpace> varSpaceList_weight;
    //Converted query to fit FastBit names
    std::string queryFB;
    std::vector<VarInfo> varInfoList_query;
    std::vector<VarSpace> varSpaceList_query;

    // Parse query and retrieve queryFB
    if (! metadataMgr->parseQuery
        (query, varPathStr, queryFB, varInfoList_query, varSpaceList_query)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << query << ", " << varPathStr.c_str() << "):"
            << " not a valid query";
        return 0;
    }

    // verify if variable name exists and
    // pick the best fit variable if there are multiple
    unsigned int numVar = 0;
    int idx = 0;
    int varIdx = -1;
    numVar = metadataMgr->getAllVariables
        (varInfoList_var1, varSpaceList_var1, varPathStr, varNameStr1);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var1[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var1[i].getPath().size()) {
                idx = i;
                len = varInfoList_var1[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr1.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var1[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var1[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var1[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var1[idx]);
        varSpaceList_query.push_back(varSpaceList_var1[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str1[25];
    sprintf(str1, "%u", varIdx);
    varFBName1 += str1;

    // Process the second variable
    numVar = metadataMgr->getAllVariables
        (varInfoList_var2, varSpaceList_var2, varPathStr, varNameStr2);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var2[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var2[i].getPath().size()) {
                idx = i;
                len = varInfoList_var2[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get2DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr2.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var2[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var2[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var2[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var2[idx]);
        varSpaceList_query.push_back(varSpaceList_var2[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str2[25];
    sprintf(str2, "%u", varIdx);
    varFBName2 += str2;

    // Process the third variable
    numVar = metadataMgr->getAllVariables
        (varInfoList_var3, varSpaceList_var3, varPathStr, varNameStr3);
    if (numVar == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr3.c_str() << "):"
            << " no valid variable is found";
        return 0;
    }

    idx = 0;
    if (numVar > 1) {
        // pick the best fit variable with the shortest length
        std::string varName = varInfoList_var3[0].getPath();
        int len = varName.size();
        for (unsigned int i=1; i<numVar; i++) {
            if (len > varInfoList_var3[i].getPath().size()) {
                idx = i;
                len = varInfoList_var3[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get3DHistogram ("
            << varPathStr.c_str() << ", " << varNameStr3.c_str() << "):"
            << " multiple variables are found. use the best fit variable \""
            << varInfoList_var3[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    varIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_var3[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_var3[idx])) {
                varIdx = (int) i;
            }
        }
    }

    // if variable name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (varIdx == -1) {
        varIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_var3[idx]);
        varSpaceList_query.push_back(varSpaceList_var3[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str3[25];
    sprintf(str3, "%u", varIdx);
    varFBName3 += str3;

    // verify if weight name exists and
    // pick the best fit weight name if there are multiple
    unsigned int numWts = metadataMgr->getAllVariables
        (varInfoList_weight, varSpaceList_weight, varPathStr, wtNameStr);
    if (numWts == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " no valid weight is found";
        return 0;
    }

    idx = 0;
    if (numWts > 1) {
        // pick the best fit weight with the shortest length
        std::string wtName = varInfoList_weight[0].getPath();
        int wlen = wtName.size();
        for (unsigned int i = 1; i < numWts; i++) {
            if (wlen > varInfoList_weight[i].getPath().size()) {
                idx = i;
                wlen = varInfoList_weight[i].getPath().size();
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- QueryProcessor::get1DHistogram ("
            << varPathStr.c_str() << ", " << wtNameStr.c_str() << "):"
            << " multiple weight variables are found.  Use the variable \""
            << varInfoList_weight[idx].getPath() << "\"";
    }

    // Check if variable name appeared in queryFB
    int wtIdx = -1;
    for (unsigned int i = 0; i < varInfoList_query.size(); i++) {
        if (varInfoList_query[i].getPath().compare
            (varInfoList_weight[idx].getPath()) == 0) {
            if (varSpaceList_query[i].isEqual(varSpaceList_weight[idx])) {
                wtIdx = (int) i;
            }
        }
    }

    // if weight name did not appear,
    // add the variable to varInfoList_query and to varSpaceList_query
    if (wtIdx == -1) {
        wtIdx = varInfoList_query.size();
        varInfoList_query.push_back(varInfoList_weight[idx]);
        varSpaceList_query.push_back(varSpaceList_weight[idx]);
    }

    // Convert the name of variable to fit FastBit convention
    char str4[25];
    sprintf(str4, "%u", wtIdx);
    weightFBName += str4;

#ifdef FQ_NOMPI
    // Initialize FQ_Part
    FQ_Part part(varInfoList_query, varSpaceList_query, *dataFile,
                 *indexFile);

    part.get3DHist_weights(queryFB, varFBName1, begin1, end1, stride1,
                           varFBName2, begin2, end2, stride2,
                           varFBName3, begin3, end3, stride3,
                           weightFBName, weights);
#else
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get3DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " MPI version not implemented";
#endif
    return 1;
}  //QueryProcessor::get3DHistogram: weights

uint64_t QueryProcessor::get3DHistogram
(const std::string& varNameStr1,
 const std::string& varNameStr2,
 const std::string& varNameStr3,
 const std::string& varPathStr,
 uint32_t nbin1,
 uint32_t nbin2,
 uint32_t nbin3,
 std::vector<double> &bounds1,
 std::vector<double> &bounds2,
 std::vector<double> &bounds3,
 std::vector<uint32_t> &counts) {
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get3DHistogram ("
        << varPathStr.c_str() << "):"
        << " Function is not implemented";
    return 1;
} //QueryProcessor::get3DHistogram: adaptive bins

uint64_t QueryProcessor::get3DHistogram
(const char *query,
 const std::string& varNameStr1,
 const std::string& varNameStr2,
 const std::string& varNameStr3,
 const std::string& varPathStr,
 uint32_t nbin1,
 uint32_t nbin2,
 uint32_t nbin3,
 std::vector<double> &bounds1,
 std::vector<double> &bounds2,
 std::vector<double> &bounds3,
 std::vector<uint32_t> &counts) {
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- QueryProcessor::get3DHistogram ("
        << query << ", " << varPathStr.c_str() << "):"
        << " Function is not implemented";
    return 1;
} //QueryProcessor::get3DHistogram: conditional adaptive bins

int QueryProcessor::recordRegions(const std::string &outputfile,
                                  const std::string &path,
                                  const std::string &cond,
                                  const std::string &col,
                                  const std::vector<double> &thr) const {
    if (col.empty() || thr.empty()) return 0; // nothing to do

    std::string fbcond;
    std::vector<VarInfo> infolist;
    std::vector<VarSpace> spacelist;
    int ierr;
    int cnt = 0;
    std::ofstream ofs(outputfile.c_str());
    std::ostream &outstr(ofs ? ofs : std::cout);

    LOGGER(ibis::gVerbose > 2)
        << "QueryProcessor::recordRegions attempting to identify connected "
        "regions with " << thr.size()
        << (thr.size()>1?" different thresholds":" threshold")
        << " on variable " << col;
    do { // while (moresteps)
        for (unsigned j = 0; j < thr.size(); ++ j) {
            // compose the query string
            std::string whereclause;
            {
                std::ostringstream oss;
                if (! cond.empty()) {
                    oss << cond << " && ";
                }
                oss << col << " >= " << thr[j];
                whereclause = oss.str();
            }
            ierr = metadataMgr->parseQuery(whereclause, path, fbcond,
                                           infolist, spacelist);
            if (ierr == 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- metadataMgr->parseQuery(" << cond
                    << ") failed while processing \"" << whereclause << '"';
                return -1;
            }
            LOGGER(ibis::gVerbose > 3)
                << "QueryProcessor::recordRegions converted \"" << whereclause
                << "\" into \"" << fbcond << '"';

            FQ_Part mypart(infolist, spacelist, *dataFile, *indexFile);
            const char *tok = mypart.createQuery(fbcond);
            ierr = mypart.submitQuery(tok);
            LOGGER(ibis::gVerbose > 0 && ierr < 0)
                << "Warning -- QueryProcessor::recordRegions failed to "
                "evaluate \"" << whereclause << '"';
            if (ierr <= 0) continue;
            LOGGER(ibis::gVerbose > 3)
                << "QueryProcessor::recordRegions -- \"" << whereclause
                << "\" selected " << ierr << " cell" << (ierr > 1 ? "s" : "");

            const ibis::meshQuery *mq = mypart.getQuery(tok);
            if (mq == 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- QueryProcessor::recordRegions failed to "
                    "retireve the query object";
                continue;
            }

            const unsigned nhits = mq->getNumHits();
            std::vector<uint32_t> lines;
            std::vector<uint32_t> labels;
            const std::vector<uint32_t> &dims = mypart.getMeshShape();
            const unsigned int nd = dims.size();
            if (nd == 0) return 0;

            ierr = mq->getHitsAsLines(lines, dims);
            if (ierr <= 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- QueryProcessor::recordRegions failed to "
                    "produce lines for \"" << whereclause << '"';
                continue;
            }
            if (ibis::gVerbose > 3) {
                ibis::util::logger lg;
                lg() << "QueryProcessor::recordRegions generated " << ierr
                     << " line" << (ierr>1?"s":"") << nhits << " hit"
                     << (nhits>1?"s":"") << " for " << whereclause
                     << " on a table with " << mypart.nRows()
                     << " rows defined on a mash of " << dims[0];
                for (unsigned ii = 1; ii < nd; ++ ii)
                    lg() << " x " << dims[ii];
            }
            ierr = mq->labelLines(nd, lines, labels);
            if (ierr <= 0) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- QueryProcessor::recordRegions failed to "
                    "produce connected components for \"" << whereclause << '"';
                continue;
            }
            if (labels.size()*(1+nd) != lines.size()) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- QueryProcessor::recordRegions produced "
                    "invalid connected components for \"" << whereclause << '"';
                continue;
            }

            cnt += ierr;
            outstr << "\n# QueryProcessor identified " << ierr
                   << " connected region" << (ierr>1?"s":"") << " containing "
                   << labels.size() << " line" << (labels.size()>1?"s":"")
                   << " and " << nhits << " hit" << (nhits>1?"s":"")
                   << " satisfying " << whereclause << "\n"
                   << ierr << ", " << dims[0];
            for (unsigned ii = 1; ii < nd; ++ ii)
                outstr << ", " << dims[ii];
            outstr << ", " << labels.size();
            for (unsigned i0 = 0; i0 < labels.size(); ++ i0) {
                outstr << "\n" << labels[i0];
                for (unsigned i1 = i0*(nd+1); i1 <= i0*(nd+1)+nd; ++ i1) {
                    outstr << ", " << lines[i1];
                }
            }

            LOGGER(ibis::gVerbose > 2)
                << "QueryProcessor::recordRegions processed regions \""
                << whereclause << "\" and wrote " << ierr << " region"
                << (ierr>1?"s":"") << " to "
                << (ofs ? outputfile.c_str() : "stdout");
        }
    } while (dataFile->nextStep());

    return cnt;
} // QueryProcessor::recordRegions

