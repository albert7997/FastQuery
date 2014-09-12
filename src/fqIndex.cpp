// $Id$
#include "fqIndex.h"
#include "fq.h"         // FastQuery::reportTiming()
#include <cmath>        // std::log

/// Constructor.
FQ_IndexUnbinned::FQ_IndexUnbinned(const FQ_Variable* c, bool readOnly)
    : ibis::relic(static_cast<ibis::column*>(0)), status(-1),
      varInfo(c->getVarInfo()), varSpace(c->getVarSpace()),
      mpi_idx(0), mpi_iter(0), mpi_max_iter(0), mpi_size(1), mpi_rank(0)
{
    if (c == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned: incomplete initialization.  "
            "The constructor needs a valid FQ_Variable pointer.";
        return;
    }
    col = c;
    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    nm = varInfo.getPath();
    if (! c->isAll()) {
        nm += varSpace.getText();
    }
    m_type = c->type();

#ifndef FQ_NOMPI
    mpi_idx = varSpace.getMpiIdx();
    mpi_iter = varSpace.getMpiIter();
    mpi_max_iter = varSpace.getMpiMaxIter();
    mpi_comm = varSpace.getMpiComm();
    MPI_Group mpi_group;
    MPI_Comm_group(mpi_comm, &mpi_group);
    MPI_Group_size(mpi_group, &mpi_size);
    MPI_Group_rank(mpi_group, &mpi_rank);
#endif // FQ_NOMPI

    std::string fbkey = "FastQuery.";
    fbkey += "forceIndexRebuild";

    bool err = ibis::gParameters().isTrue(fbkey.c_str());
    LOGGER(ibis::gVerbose > 2 && err)
        << "FQ_IndexUnbinned[" << nm.c_str() << "]: "
        << "forceIndexRebuild is set.";

    bool built = false;
    if (err == false) {
        built = readOld(indexFile);
    }
    if (built == false && readOnly == false) {
        buildNew();
    }
    if (ibis::gVerbose > 5) { // output a summary of the index
        ibis::util::logger lg;
        print(lg());
    }
} // FQ_IndexUnbinned::FQ_IndexUnbinned

/// Build a new index in-memory.
bool FQ_IndexUnbinned::buildNew() {
    ibis::horometer timer;
    ibis::horometer readTimer;
    ibis::horometer buildTimer;
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    ibis::horometer initIdleTimer;
    ibis::horometer readIdleTimer;
    ibis::horometer buildIdleTimer;

    initIdleTimer.start();
    MPI_Barrier(mpi_comm);
    initIdleTimer.stop();
#endif

    timer.start();
    readTimer.start();
    int ierr;
    switch (m_type) {
    case ibis::INT: {
        ibis::array_t<int32_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::DOUBLE: {
        ibis::array_t<double> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::BYTE: {
        ibis::array_t<signed char> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::buildNew: does not yet support data type "
            << ibis::TYPESTRING[(int)m_type];
        ierr = -3;
        break;}
    }
    if (ierr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::buildNew: failed to get values array with error code "
            << ierr;
        return false;
    }
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    buildIdleTimer.start();
    MPI_Barrier(mpi_comm);
    buildIdleTimer.stop();
#endif
    buildTimer.stop();

    // number of words in serialized bitmaps
    const size_t nobs = bits.size();
    if (nobs == 0) {
        status = 1;
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexUnbinned[" << nm.c_str() << "]::buildNew: "
            << "successfully build new unbinned indexes but size is 0";
    } else {
        offset64.resize(nobs+1);
        offset64[0] = 0;
        for (unsigned i = 0; i < nobs; ++ i) {
            offset64[i+1] = offset64[i];
            const ibis::bitvector* tmp = bits[i];
            if (tmp != 0) {
                const ibis::bitvector::word_t w = tmp->getSerialSize();
                if (w > 1) {
                    offset64[i+1] += (w >> 2);
                }
            }
        }
    }
    timer.stop();

    LOGGER(FastQuery::reportTiming())
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        << "\nStatistic\tFQ_IndexUnbinned::buildNew:init_idle\t"
        << initIdleTimer.CPUTime() << "\t" << initIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::buildNew:read_idle\t"
        << readIdleTimer.CPUTime() << "\t" << readIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::buildNew:build_idle\t"
        << buildIdleTimer.CPUTime() << "\t" << buildIdleTimer.realTime()
#endif
        << "\nStatistic\tFQ_IndexUnbinned::buildNew:getValuesArray()\t"
        << readTimer.CPUTime() << "\t" << readTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::buildNew:construct()\t"
        << buildTimer.CPUTime() << "\t" << buildTimer.realTime()
        << "\t" << nobs
        << "\nStatistic\tFQ_IndexUnbinned::buildNew\t"
        << timer.CPUTime() << "\t" << timer.realTime() ;
    status = 1;
    return true;
} // FQ_IndexUnbinned::buildNew

bool FQ_IndexUnbinned::readOld(ArrayIODriver& indexFile) {
    ibis::horometer timer;
    timer.start();

    uint64_t nobs = 0;
    bool berr = indexFile.getBitmapKeyLength(nm, &nobs, mpi_idx);
    if (berr == false || nobs == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "FQ_IndexUnbinned[" << nm.c_str() << "]::readOld: "
            << "no existing indexes to read in file \""
            << indexFile.getFileName() << "\"";
        return false;
    }
    uint64_t noffsets = 0;
    berr = indexFile.getBitmapOffsetLength(nm, &noffsets, mpi_idx);
    if (berr == false || noffsets == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexUnbinned[" << nm.c_str() << "]::readOld: "
            << "no existing bitmap offsets in file \""
            << indexFile.getFileName() << "\"";
        return false;
    }
    if ( (nobs+1) != noffsets ) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::readOld: the number of keys (" << nobs << ")"
            << " does not match to the number of offsets ("
            << noffsets << ")\n\t #offset should be equal to (#key+1)";
        nobs = 0;
        return false;
    }

    clear();
    // read in the keys (store the resulting values in vals)
    vals.resize(nobs);

    switch (m_type) {
    case ibis::INT: {
        ibis::array_t<int32_t> buf(nobs);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            vals.resize(nobs);
            for (unsigned i = 0; i < nobs; ++ i)
                vals[i] = buf[i];
        }
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> buf(nobs);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            vals.resize(nobs);
            for (unsigned i = 0; i < nobs; ++ i)
                vals[i] = buf[i];
        }
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> buf(nobs);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            vals.resize(nobs);
            for (unsigned i = 0; i < nobs; ++ i)
                vals[i] = buf[i];
        }
        break;}
    case ibis::DOUBLE: {
        vals.resize(nobs);
        berr = indexFile.getBitmapKeys(nm, vals.begin(), mpi_idx);
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> buf(nobs);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            vals.resize(nobs);
            for (unsigned i = 0; i < nobs; ++ i)
                vals[i] = buf[i];
        }
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> buf(nobs);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            vals.resize(nobs);
            for (unsigned i = 0; i < nobs; ++ i)
                vals[i] = buf[i];
        }
        break;}
    case ibis::BYTE: {
        ibis::array_t<char> buf(nobs);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            vals.resize(nobs);
            for (unsigned i = 0; i < nobs; ++ i)
                vals[i] = buf[i];
        }
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::readOld: does not yet support data type "
            << ibis::TYPESTRING[(int)m_type];
        berr = false;
    }
    }
    if (berr == false) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::readOld: failed to read BitmapKeys from file \""
            << indexFile.getFileName() << "\"";
        clear();
        return false;
    }

    // read in the offsets
    FQ::DataType offtype = indexFile.getBitmapOffsetType(nm);
    switch(offtype) {
    case FQ::FQT_INT: {
        offset64.clear();
        offset32.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset32.begin(), mpi_idx);
        break;}
    case FQ::FQT_LONG: {
        offset32.clear();
        offset64.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset64.begin(), mpi_idx);
        break;}
    case FQ::FQT_UINT: {
        offset64.clear();
        offset32.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset32.begin(), mpi_idx);
        break;}
    case FQ::FQT_ULONG: {
        offset32.clear();
        offset64.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset64.begin(), mpi_idx);
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]::readOld: "
            << "unsupported data type for BitmapOffsets";
        berr = false;
        break;}
    }
    if (berr == false) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]::readOld: "
            << "failed to read BitmapOffsets from file \""
            << indexFile.getFileName() << "\"";
        clear();
        return false;
    }

    // activate the first bitmap so that the size information is available
    bits.resize(nobs);
    for (unsigned i = 0; i < nobs; ++ i)
        bits[i] = 0;
    //activate(0);

    std::vector<uint64_t> dims = varSpace.getCounts();
    if (! dims.empty()) {
        nrows = 1;
        for (unsigned i = 0; i < dims.size(); ++ i)
            nrows *= dims[i];
    }
    else {
        nrows = 0;
    }
    fname = ibis::util::strnewdup(nm.c_str());
    status = 0;

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexUnbinned::readOld\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;
    return true;
} // FQ_IndexUnbinned::readOld

/// Write the index to an active HDF5 file.
bool FQ_IndexUnbinned::write(ArrayIODriver& outIndexFile) const
{
    ibis::horometer timer;
    ibis::horometer setupTimer;
    ibis::horometer metadataTimer;
    ibis::horometer bitmapTimer;
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    ibis::horometer initIdleTimer;
    ibis::horometer bitmapIdleTimer;
    ibis::horometer metadataIdleTimer;
    ibis::horometer setupIdleTimer;

    initIdleTimer.start();
    MPI_Barrier(mpi_comm);
    initIdleTimer.stop();
#endif

    timer.start();
    const size_t nobs = bits.size();
#ifdef FQ_NOMPI
    if (nobs == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexUnbinned[" << nm.c_str()
            << "]::write has no bitmaps to write";
        return true;
    }
#endif

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    bool berr = true;

    if (status <= 0) {
        if (&outIndexFile == &indexFile) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexUnbinned[" << nm.c_str()
                << "]::write skips writing to the same file \""
                << indexFile.getFileName() << "\"";
            return true;
        }
        activate();
    }

    unsigned int bitmapLen = 0;
    if (nobs > 0) {
        if (offset64.size() <= nobs) {
            offset64.resize(nobs+1);
            offset64[0] = 0;
            for (unsigned i = 0; i < nobs; ++ i) {
                offset64[i+1] = offset64[i];
                if (bits[i] != 0) {
                    const ibis::bitvector::word_t w = bits[i]->getSerialSize();
                    if (w > 1) offset64[i+1] += (w >> 2);
                }
            }
        }
        bitmapLen = offset64[nobs];
    }

    setupTimer.resume();
#ifndef FQ_NOMPI
    if (mpi_idx < mpi_size) {
        // in the first iteration, create offset dataset collectively
        berr = outIndexFile.createOffsetTable(nm, mpi_max_iter);
        if (! berr) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
                << "]::write failed to create offset dataset";
            return false;
        }
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexUnbinned[" << nm.c_str()
            << "]::write created offset dataset";
    }

    berr = outIndexFile.setBitmapLength(nm, bitmapLen, mpi_iter);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to set bitmapLen";
        return false;
    }
#endif

    berr = outIndexFile.createBitmap(nm, bitmapLen, mpi_iter, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to create bitmap dataset with bitmap length "
            << bitmapLen;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexUnbinned[" << nm.c_str()
        << "]::write created bitmap dataset with length " << bitmapLen;

    unsigned int bitmapOffsetLen = 0;
    if (nobs > 0) bitmapOffsetLen = nobs+1;

#ifndef FQ_NOMPI
    metadataTimer.resume();
    berr = outIndexFile.setBitmapOffsetLength(nm, bitmapOffsetLen, mpi_iter);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to set bitmapOffsetLen";
        return false;
    }
#endif

    berr = outIndexFile.createBitmapOffsets(nm, bitmapOffsetLen, FQ::FQT_LONG,
                                            mpi_iter, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to create bitmapOffset dataset with length "
            << bitmapOffsetLen;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexUnbinned[" << nm.c_str()
        << "]::write created bitmapOffset dataset with length "
        << bitmapOffsetLen;
    setupTimer.stop();

    unsigned int nkeys = 0;
    if (nobs > 0) nkeys = nobs;
#ifndef FQ_NOMPI
    metadataTimer.resume();
    berr = outIndexFile.setBitmapKeyLength(nm, nkeys, mpi_iter);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to set bitmapKeyLength";
        return false;
    }
#endif

    berr = outIndexFile.createBitmapKeys
        (nm, nkeys, varInfo.getType(), mpi_iter, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to create bitmapKey dataset with length "
            << nkeys;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexUnbinned[" << nm.c_str()
        << "]::write created bitmapKey dataset with length " << nkeys;
    setupTimer.stop();

#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    setupIdleTimer.start();
    MPI_Barrier(mpi_comm);
    setupIdleTimer.stop();
#endif
    setupTimer.stop();

    metadataTimer.start();
    // collectively/independently write bitmapOffsets
    if (nobs > 0) {
        berr = outIndexFile.setBitmapOffsets
            (nm, const_cast<int64_t*>(offset64.begin()),
             nobs+1, mpi_idx);
    } else {
        // to make sure it can be handled collectively with other MPI tasks
        berr = outIndexFile.setBitmapOffsets
            (nm, NULL, 0, mpi_idx);
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write failed to record the bitmap offsets for file \""
            << outIndexFile.getFileName() << "\"";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexUnbinned[" << nm.c_str()
        << "]::write recorded bitmapOffsets as destined file \""
        << outIndexFile.getFileName() << "\"";
    metadataTimer.stop();

    // collectively/independently write the keys, need return the min/max
    // values to their original type
    metadataTimer.resume();
    bool berr1 = true;
    bool berr2 = true;
    switch(m_type) {
    case ibis::DOUBLE: {
#ifdef FQ_NOMPI
        double range[2];
        range[0] = vals[0];
        range[1] = vals[nobs-1];
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_DOUBLE);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_DOUBLE);
#endif
        berr = outIndexFile.setBitmapKeys(nm, vals.begin(), nobs, mpi_idx);
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> tmp(vals.size());
        for (unsigned i = 0; i < vals.size(); ++ i)
            tmp[i] = static_cast<float>(vals[i]);
#ifdef FQ_NOMPI
        float range[2];
        range[0] = static_cast<float>(vals[0]);
        range[1] = static_cast<float>(vals[nobs-1]);
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_FLOAT);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_FLOAT);
#endif
        berr = outIndexFile.setBitmapKeys(nm, tmp.begin(), nobs, mpi_idx);
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> tmp(vals.size());
        for (unsigned i = 0; i < vals.size(); ++ i)
            tmp[i] = static_cast<int64_t>(vals[i]);
#ifdef FQ_NOMPI
        int64_t range[2];
        range[0] = static_cast<int64_t>(vals[0]);
        range[1] = static_cast<int64_t>(vals[nobs-1]);
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_LONG);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_LONG);
#endif
        berr = outIndexFile.setBitmapKeys(nm, tmp.begin(), nobs, mpi_idx);
        break;}
    case ibis::INT: {
        ibis::array_t<int32_t> tmp(vals.size());
        for (unsigned i = 0; i < vals.size(); ++ i) {
            tmp[i] = static_cast<int32_t>(vals[i]);
        }
#ifdef FQ_NOMPI
        int32_t range[2];
        range[0] = static_cast<int32_t>(vals[0]);
        range[1] = static_cast<int32_t>(vals[nobs-1]);
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_INT);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_INT);
#endif
        berr = outIndexFile.setBitmapKeys(nm, tmp.begin(), nobs, mpi_idx);
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> tmp(vals.size());
        for (unsigned i = 0; i < vals.size(); ++ i)
            tmp[i] = static_cast<uint64_t>(vals[i]);
#ifdef FQ_NOMPI
        uint64_t range[2];
        range[0] = static_cast<uint64_t>(vals[0]);
        range[1] = static_cast<uint64_t>(vals[nobs-1]);
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_ULONG);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_ULONG);
#endif
        berr = outIndexFile.setBitmapKeys(nm, tmp.begin(), nobs, mpi_idx);
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> tmp(vals.size());
        for (unsigned i = 0; i < vals.size(); ++ i)
            tmp[i] = static_cast<uint32_t>(vals[i]);
#ifdef FQ_NOMPI
        uint32_t range[2];
        range[0] = static_cast<uint32_t>(vals[0]);
        range[1] = static_cast<uint32_t>(vals[nobs-1]);
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_UINT);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_UINT);
#endif
        berr = outIndexFile.setBitmapKeys(nm, tmp.begin(), nobs, mpi_idx);
        break;}
    case ibis::BYTE: {
        ibis::array_t<char> tmp(vals.size());
        for (unsigned i = 0; i < vals.size(); ++ i)
            tmp[i] = static_cast<char>(vals[i]);
#ifdef FQ_NOMPI
        char range[2];
        range[0] = static_cast<char>(vals[0]);
        range[1] = static_cast<char>(vals[nobs-1]);
        berr1 = outIndexFile.setActualRange(nm, range, FQ::FQT_UINT);
        berr2 = outIndexFile.setExpectedRange(nm, range, FQ::FQT_UINT);
#endif
        berr = outIndexFile.setBitmapKeys(nm, tmp.begin(), nobs, mpi_idx);
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
            << "]::write does not yet support data type "
            << ibis::TYPESTRING[(int)m_type];
    }
        return false;
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]::write "
            << "failed to write the bitmap keys array[" << nobs <<"]";
        return false;
    }
    if (! berr1) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]::write "
            << "failed to write the actual key range";
        return false;
    }
    if (! berr2) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]::write "
            << "failed to write the expected key range";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexUnbinned[" << nm.c_str()
        << "]::write recorded bitmapKeys as destined for file \""
        << outIndexFile.getFileName() << "\"";
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    metadataIdleTimer.resume();
    MPI_Barrier(mpi_comm);
    metadataIdleTimer.stop();
#endif
    metadataTimer.stop();
    // go through the bitmaps and write them by chunk
    bitmapTimer.start();
#ifdef FQ_PACK_BITMAP
    ibis::array_t<ibis::bitvector::word_t> chunk;
    uint64_t head = 0;
    for (unsigned i = 0; i < nobs; ++ i) {
        ibis::array_t<ibis::bitvector::word_t> buf;
        bits[i]->write(buf);
        uint64_t pos = offset64[i] - head;
        if (offset64[i+1]-head > FQ_CHUNK_SIZE) {
            // break the bitmap into two chunks
            uint64_t len = FQ_CHUNK_SIZE-(offset64[i]-head);
            chunk.insert(chunk.begin()+pos, buf.begin(), buf.begin()+len);
            // write the current chunk
            if (! outIndexFile.writeBitmap
                (nm, head, head+chunk.size(),
                 reinterpret_cast<uint32_t*>(chunk.begin()), mpi_idx)) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
                    << "]::write failed to write bitvector chunk";
                return false;
            }
            // restart with a new chunk
            chunk.clear();
            chunk.insert(chunk.begin(), buf.begin()+len, buf.end());
            head += FQ_CHUNK_SIZE;
        } else {
            uint64_t len = offset64[i+1] - offset64[i];
            chunk.insert(chunk.begin()+pos, buf.begin(), buf.begin()+len);
        }
    }
    if (nobs > 0) {
        if (! outIndexFile.writeBitmap
            (nm, head, head+chunk.size(),
             reinterpret_cast<uint32_t*>(chunk.begin()), mpi_idx)) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
                << "]::write failed to write bitvector chunk";
            return false;
        }
    }
#else
    // go through the bitmaps to write each one separately
    for (unsigned i = 0; i < nobs; ++ i) {
        if (offset64[i] < offset64[i+1]) {
            ibis::array_t<ibis::bitvector::word_t> buf;
            bits[i]->write(buf);
            berr = outIndexFile.writeBitmap
                (nm, offset64[i], offset64[i+1],
                 reinterpret_cast<uint32_t*>(buf.begin()), mpi_idx);
            if (! berr) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
                    << "]::write failed to write bitvector "
                    << i << " (" << offset64[i+1]-offset64[i] << " words)";
                return false;
                // is there a way to rollback the write operations?
            }
        }
    }
#endif

#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    bitmapIdleTimer.start();
    MPI_Barrier(mpi_comm);
    bitmapIdleTimer.stop();
#endif
    bitmapTimer.stop();
    timer.stop();

    LOGGER(FastQuery::reportTiming())
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        << "\nStatistic\tFQ_IndexUnbinned::write:init_idle\t"
        << initIdleTimer.CPUTime() << "\t" << initIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::write:setup_idle\t"
        << setupIdleTimer.CPUTime() << "\t" << setupIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::write:metadata_idle\t"
        << metadataIdleTimer.CPUTime() << "\t" << metadataIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::write:bitmap_idle\t"
        << bitmapIdleTimer.CPUTime() << "\t" << bitmapIdleTimer.realTime()
#endif
        << "\nStatistic\tFQ_IndexUnbinned::write:setup\t"
        << setupTimer.CPUTime() << "\t" << setupTimer.realTime()
        << "\nStatistic\tFQ_IndexUnbinned::write:metadata\t"
        << metadataTimer.CPUTime() << "\t" << metadataTimer.realTime()
        << "\t" << nobs
        << "\nStatistic\tFQ_IndexUnbinned::write:bitmap\t"
        << bitmapTimer.CPUTime() << "\t" << bitmapTimer.realTime()
        << "\t" << bitmapLen
        << "Statistic\tFQ_IndexUnbinned::write\t"
        << timer.CPUTime() << "\t" << timer.realTime();
    return true;
} // FQ_IndexUnbinned::write

/// Activate all bitmaps at once by reading all of them into memory.
void FQ_IndexUnbinned::activate() const {
    ibis::horometer timer;
    timer.start();

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    // assumes that vals, bits and offsets have been initialized properly.
    const unsigned nobs = vals.size();
    if (nobs == 0 || nobs > bits.size() ||
        (offset64.size() <= nobs && offset32.size() <= nobs)) return;

    bool berr;
    if (bits[0] == 0 && str == 0) {
        // attempt to read all bitmaps into memory and store the bytes in str
        const size_t nwords = (offset64.size()>nobs ? offset64[nobs] :
                               (int64_t)offset32[nobs]);
        str = new ibis::fileManager::storage
            (nwords*sizeof(ibis::bitvector::word_t));
        berr = indexFile.readBitmap(nm, 0, nwords, (uint32_t*)(str->begin()),
                                    mpi_idx);
        if (berr == false) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexUnbinned[" << nm.c_str()
                << "]::activate: unable to read all bitmaps at once";
            delete str;
            str = 0;
        }
    }

    if (str) {
        if (offset64.size() > nobs) {
            for (unsigned i = 0; i < nobs; ++ i) {
                if (bits[i] == 0 && offset64[i+1] > offset64[i]) {
                    ibis::array_t<ibis::bitvector::word_t>
                        buf(str, offset64[i]*sizeof(ibis::bitvector::word_t),
                            offset64[i+1]*sizeof(ibis::bitvector::word_t));
                    bits[i] = new ibis::bitvector(buf);
                }
            }
        }
        else {
            for (unsigned i = 0; i < nobs; ++ i) {
                if (bits[i] == 0 && offset32[i+1] > offset32[i]) {
                    ibis::array_t<ibis::bitvector::word_t>
                        buf(str, offset32[i]*sizeof(ibis::bitvector::word_t),
                            offset32[i+1]*sizeof(ibis::bitvector::word_t));
                    bits[i] = new ibis::bitvector(buf);
                }
            }
        }
    }
    else { // read one bitmap at a time
        for (unsigned i = 0; i < nobs; ++ i) {
            uint64_t start, end;
            if (offset64.size() > nobs) {
                start = offset64[i];
                end = offset64[i+1];
            }
            else {
                start = offset32[i];
                end = offset32[i+1];
            }
            if (bits[i] == 0 && end > start) {
                ibis::array_t<ibis::bitvector::word_t> buf(end-start);
                berr = indexFile.readBitmap
                    (nm, start, end, (uint32_t*)(buf.begin()), mpi_idx);
                if (berr) {
                    bits[i] = new ibis::bitvector(buf);
                }
                else {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
                        << "]::activate: failed to read bitmap " << i
                        << " (offset " << start << ", size "
                        << end-start << ")";
                }
            }
        }
    }
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexUnbinned::active()\t"
        << timer.CPUTime() << "\t" << timer.realTime() ;
} // FQ_IndexUnbinned::activate

/// Activate the ith bitmap.
void FQ_IndexUnbinned::activate(uint32_t i) const {
    ibis::horometer timer;
    timer.start();

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    const unsigned nobs = vals.size();

    if (i >= nobs) return;
    if (bits[i] != 0) return;
    if (bits.size() != nobs) return;
    if (offset64.size() <= nobs && offset32.size() <= nobs) return;
    if (offset64.size() > nobs ? offset64[i] >= offset64[i+1] :
        offset32[i] >= offset32[i+1]) {
        return;
    }

    bool berr;
    const int64_t sz1 = (offset64.size()>nobs ? (offset64[1]-offset64[0]) :
                         (int64_t)(offset32[1]-offset32[0]));
    const int64_t sza = (offset64.size()>nobs ? (offset64[nobs]-offset64[0]) :
                         (int64_t)(offset32[nobs]-offset32[0]));
    if (i == 0 && (nobs == 1U || sz1*5/4 >= sza
                   || (int64_t)(sz1*std::log((double)nobs)) >= sza)) {
        // if the first bitmap takes up a majority of the total bytes, read
        // them all
        const size_t nwords = (offset64.size() > nobs ? offset64[nobs] :
                               (int64_t)offset32[nobs]);
        str = new ibis::fileManager::storage
            (nwords*sizeof(ibis::bitvector::word_t));
        berr = indexFile.readBitmap
            (nm, 0, nwords, (uint32_t*)(str->begin()), mpi_idx);
        if (berr == false) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexUnbinned[" << nm.c_str() << "]"
                << "::activate" << "(" << i << "): "
                << "unable to read all bitmaps at once";
            delete str;
            str = 0;
        }
    }

    if (offset64.size() > nobs) {
        if (str) { // already in memory
            ibis::array_t<ibis::bitvector::word_t> buf
                (str, offset64[i]*sizeof(ibis::bitvector::word_t),
                 offset64[i+1]*sizeof(ibis::bitvector::word_t));
            bits[i] = new ibis::bitvector(buf);
        }
        else {
            //ibis::array_t<ibis::bitvector::word_t>
            ibis::array_t<ibis::bitvector::word_t>
                buf(offset64[i+1]-offset64[i]);
            berr = indexFile.readBitmap(nm, offset64[i], offset64[i+1],
                                        (uint32_t*)(buf.begin()), mpi_idx);
            if (berr) {
                bits[i] = new ibis::bitvector(buf);
            }
            else {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]"
                    << "::activate" << "(" << i << "): "
                    << " failed to read bitmap " << i
                    << " (offset " << offset64[i] << ", size "
                    << offset64[i+1]-offset64[i] << ")";
            }
        }
    }
    else if (str) { // already in memory
        ibis::array_t<ibis::bitvector::word_t> buf
            //(str, offset32[i]*sizeof(ibis::bitvector::word_t),
            // offset32[i+1]*sizeof(ibis::bitvector::word_t));
            (str, offset32[i]*sizeof(uint32_t),
             offset32[i+1]*sizeof(uint32_t));
        bits[i] = new ibis::bitvector(buf);
    }
    else {
        ibis::array_t<ibis::bitvector::word_t> buf(offset32[i+1]-offset32[i]);
        berr = indexFile.readBitmap(nm, offset32[i], offset32[i+1],
                                    (uint32_t*)(buf.begin()), mpi_idx);
        if (berr) {
            bits[i] = new ibis::bitvector(buf);
        }
        else {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]"
                << "::activate" << "(" << i << "):"
                << " failed to read bitmap " << i
                << " (offset " << offset32[i] << ", size "
                << offset32[i+1]-offset32[i] << ")";
        }
    }
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexUnbinned::active(" << i << ")\t"
        << timer.CPUTime() << "\t" << timer.realTime() ;
} // FQ_IndexUnbinned::activate

/// Activate the ith through jth bitmap.
void FQ_IndexUnbinned::activate(uint32_t i, uint32_t j) const {
    ibis::horometer timer;
    timer.start();

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();

    const unsigned nobs = vals.size();
    if (i >= nobs || i >= j) return; // empty range
    if (bits.size() != nobs) return;
    if (offset64.size() <= nobs && offset32.size() <= nobs) return;

    bool berr;
    const int64_t sz1 = (offset64.size()>nobs ? (offset64[1]-offset64[0]) :
                         (int64_t)(offset32[1]-offset32[0]));
    const int64_t sza = (offset64.size()>nobs ? (offset64[nobs]-offset64[0]) :
                         (int64_t)(offset32[nobs]-offset32[0]));
    if (i == 0 && (nobs == 1U || sz1*5/4 >= sza
                   || (int64_t)(sz1*std::log((double)nobs)) >= sza)) {
        // if the first bitmap takes up a majority of the total bytes, read
        // them all
        const size_t nwords = (offset64.size() > nobs ? offset64[nobs] :
                               (int64_t)offset32[nobs]);
        str = new ibis::fileManager::storage
            (nwords*sizeof(ibis::bitvector::word_t));
        berr = indexFile.readBitmap
            (nm, 0, nwords, (uint32_t*)(str->begin()), mpi_idx);
        if (berr == false) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexUnbinned[" << nm.c_str() << "]"
                << "::activate" << "(" << i << "," << j << "):"
                << " unable to read all bitmaps at once";
            delete str;
            str = 0;
        }
    }

    if (offset64.size() > nobs) {
        if (str) { // already in memory
            while (i < j) {
                if (bits[i] == 0 && offset64[i+1] > offset64[i]) {
                    ibis::array_t<ibis::bitvector::word_t> buf
                        (str, offset64[i]*sizeof(ibis::bitvector::word_t),
                         offset64[i+1]*sizeof(ibis::bitvector::word_t));
                    bits[i] = new ibis::bitvector(buf);
                }
                ++ i;
            }
        }
        else {
            while (i < j) {
                // skip to next empty bit vector
                while (i < j && bits[i] != 0)
                    ++ i;
                // the last bitvector to activate. can not be larger
                // than j
                unsigned aj = (i<j ? i + 1 : j);
                while (aj < j && bits[aj] == 0)
                    ++ aj;
                if (offset64[aj] > offset64[i]) {
                    // read bitmaps into memory in one shot
                    const uint64_t start = offset64[i];
                    ibis::array_t<ibis::bitvector::word_t>
                        buf(offset64[aj]-start);
                    berr = indexFile.readBitmap
                        (nm, offset64[i], offset64[aj],
                         (uint32_t*)(buf.begin()), mpi_idx);
                    if (berr) {
                        while (i < aj) {
                            if (bits[i] == 0 && offset64[i+1] > offset64[i]) {
                                ibis::array_t<ibis::bitvector::word_t>
                                    tmp(buf, offset64[i]-start,
                                        offset64[i+1]-start);
                                bits[i] = new ibis::bitvector(tmp);
                            }
                            ++ i;
                        }
                    }
                    else {
                        LOGGER(ibis::gVerbose > 0)
                            << "Warning -- FQ_IndexUnbinned[" << nm.c_str()
                            << "]::activate" << "(" << i << "," << j << "):"
                            << " failed to read bitmaps("
                            << offset64[i] << ", " << offset64[aj] << ")";
                    }
                    i = aj;
                } // if (offset64[aj] > offset64[i])
            } // while (i < j)
        }
    }
    else if (str) { // already in memory
        while (i < j) {
            if (bits[i] == 0 && offset32[i+1] > offset32[i]) {
                ibis::array_t<ibis::bitvector::word_t> buf
                    (str, offset32[i]*sizeof(ibis::bitvector::word_t),
                     offset32[i+1]*sizeof(ibis::bitvector::word_t));
                bits[i] = new ibis::bitvector(buf);
            }
            ++ i;
        }
    }
    else {
        while (i < j) {
            // skip to next empty bit vector
            while (i < j && bits[i] != 0)
                ++ i;
            // the last bitvector to activate. can not be larger
            // than j
            unsigned aj = (i<j ? i + 1 : j);
            while (aj < j && bits[aj] == 0)
                ++ aj;
            if (offset32[aj] > offset32[i]) {
                // read bitmaps into memory in one shot
                const uint32_t start = offset32[i];
                ibis::array_t<ibis::bitvector::word_t> buf(offset32[aj]-start);
                berr = indexFile.readBitmap(nm, offset32[i], offset32[aj],
                                            (uint32_t*)(buf.begin()), mpi_idx);
                if (berr) {
                    while (i < aj) {
                        if (bits[i] == 0 && offset32[i+1] > offset32[i]) {
                            ibis::array_t<ibis::bitvector::word_t>
                                tmp(buf, offset32[i]-start,
                                    offset32[i+1]-start);
                            bits[i] = new ibis::bitvector(tmp);
                        }
                        ++ i;
                    }
                }
                else {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- FQ_IndexUnbinned[" << nm.c_str() << "]"
                        << "::activate" << "(" << i << "," << j << "):"
                        << " failed to read bitmaps("
                        << offset32[i] << ", " << offset32[aj] << ")";
                }
                i = aj;
            } // if (offset32[aj] > offset32[i])
        } // while (i < j)
    }
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexUnbinned::active(" << i << ", " << j << ")\t"
        << timer.CPUTime() << "\t" << timer.realTime();
} // FQ_IndexUnbinned::activate

FQ_IndexBinned::FQ_IndexBinned(const FQ_Variable *c, const char *binning,
                               bool readOnly)
    : ibis::bin(static_cast<ibis::column*>(0), static_cast<const char*>(0)),
      status(-1), varInfo(c->getVarInfo()), varSpace(c->getVarSpace()),
      mpi_idx(0), mpi_iter(0), mpi_max_iter(0), mpi_size(1), mpi_rank(0)
{
    if (c == 0) { // nothing can be done if c == 0
        LOGGER(ibis::gVerbose > 0)
            << "Warning FQ_IndexBinned::ctor incomplete initialization.  "
            "The constructor needs a valid FQ_Variable pointer";
        return;
    }
    col = c;
    ArrayIODriver& indexFile = c->getIndexFile();
    nm = varInfo.getPath();
    if (! c->isAll()) {
        nm += varSpace.getText();
    }
    m_type = c->type();

#ifndef FQ_NOMPI
    mpi_idx = varSpace.getMpiIdx();
    mpi_iter = varSpace.getMpiIter();
    mpi_max_iter = varSpace.getMpiMaxIter();
    mpi_comm = varSpace.getMpiComm();
    MPI_Group mpi_group;
    MPI_Comm_group(mpi_comm, &mpi_group);
    MPI_Group_size(mpi_group, &mpi_size);
    MPI_Group_rank(mpi_group, &mpi_rank);
#endif

    std::string fbkey = "FastQuery.";
    fbkey += "forceIndexRebuild";
    bool err = ibis::gParameters().isTrue(fbkey.c_str());
    LOGGER(ibis::gVerbose > 2 && err)
        << "FQ_IndexBinned[" << nm.c_str() << "]: "
        << "forceIndexRebuild is set.";

    bool built = false;
    if (err == false) {
        // reading an existing index, no use for the binning option
        built = readOld(indexFile);
    }
    if (built == false && readOnly == false) {
        // building a new index, need to pass the binning option to the
        // actual worker function through calling ibis::column::indexSpec
        if (binning != 0 && *binning != 0) {
            // skip leading space
            while (isspace(*binning) && *binning != 0)
                ++ binning;
            if (*binning != 0) // replace the index specificiation
                const_cast<ibis::column*>(col)->indexSpec(binning);
        }
        buildNew();
    }
    if (ibis::gVerbose > 5) { // output a summary of the index
        ibis::util::logger lg;
        print(lg());
    }
} // FQ_IndexBinned::FQ_IndexBinned

bool FQ_IndexBinned::buildNew() {
    ibis::horometer timer;
    ibis::horometer readTimer;
    ibis::horometer buildTimer;
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    ibis::horometer initIdleTimer;
    ibis::horometer readIdleTimer;
    ibis::horometer buildIdleTimer;

    initIdleTimer.start();
    MPI_Barrier(mpi_comm);
    initIdleTimer.stop();
#endif

    timer.start();
    readTimer.start();
    int ierr;
    uint64_t nElements=0;
    switch (m_type) {
    case ibis::INT: {
        ibis::array_t<int32_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0) {
            construct(arr);
            nElements = arr.size();
        }
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::DOUBLE: {
        ibis::array_t<double> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    case ibis::BYTE: {
        ibis::array_t<signed char> arr;
        ierr = reinterpret_cast<const FQ_Variable*>(col)->getValuesArray(&arr);
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        readIdleTimer.start();
        MPI_Barrier(mpi_comm);
        readIdleTimer.stop();
#endif
        readTimer.stop();
        buildTimer.start();
        if (ierr >= 0)
            construct(arr);
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::buildNew: "
            << "does not yet support data type "
            << ibis::TYPESTRING[(int)m_type];
        ierr = -3;
        break;}
    }
    if (ierr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::buildNew: "
            << "failed to get values array with error code " << ierr;
        return false;
    }
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    buildIdleTimer.start();
    MPI_Barrier(mpi_comm);
    buildIdleTimer.stop();
#endif
    buildTimer.stop();

    // number of words in serialized bitmaps
    const size_t nobs = bits.size();
    if (nobs == 0) {
        status = 1;
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexBinned[" << nm.c_str() << "]::buildNew: "
            << "successfully build new binned indexes but size is 0";
    } else {
        offset64.resize(nobs+1);
        offset64[0] = 0;
        for (unsigned i = 0; i < nobs; ++ i) {
            offset64[i+1] = offset64[i];
            const ibis::bitvector* tmp = bits[i];
            if (tmp != 0) {
                const ibis::bitvector::word_t w = tmp->getSerialSize();
                if (w > 1) {
                    offset64[i+1] += (w >> 2);
                }
            }
        }
    }
    timer.stop();

    LOGGER(FastQuery::reportTiming())
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        << "\nStatistic\tFQ_IndexBinned::buildNew:init_idle\t"
        << initIdleTimer.CPUTime() << "\t" << initIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexBinned::buildNew:read_idle\t"
        << readIdleTimer.CPUTime() << "\t" << readIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexBinned::buildNew:build_idle\t"
        << buildIdleTimer.CPUTime() << "\t" << buildIdleTimer.realTime()
#endif
        << "\nStatistic\tFQ_IndexBinned::buildNew:getValuesArray()\t"
        << readTimer.CPUTime() << "\t" << readTimer.realTime() << "\t"
        << nElements
        << "\nStatistic\tFQ_IndexBinned::buildNew:construct()\t"
        << buildTimer.CPUTime() << "\t" << buildTimer.realTime() << "\t"
        << nobs
        << "Statistic\tFQ_IndexBinned::buildNew\t"
        << timer.CPUTime() << "\t" << timer.realTime() ;

    status = 1;
    return true;
} // FQ_IndexBinned::buildNew

bool FQ_IndexBinned::readOld(ArrayIODriver &indexFile) {
    ibis::horometer timer;
    timer.start();

    uint64_t nkeys = 0;
    bool berr = indexFile.getBitmapKeyLength(nm, &nkeys, mpi_idx);
    if (berr == false || nkeys == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "FQ_IndexBinned[" << nm.c_str() << "]::readOld: "
            << "no existing bitmap keys in file \""
            << indexFile.getFileName() << "\"";
        return false;
    }
    uint64_t noffsets = 0;
    berr = indexFile.getBitmapOffsetLength(nm, &noffsets, mpi_idx);
    if (berr == false || noffsets == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "FQ_IndexBinned[" << nm.c_str() << "]::readOld: "
            << "no existing bitmap offsets in file \""
            << indexFile.getFileName() << "\"";
        return false;
    }
    if ( nkeys != (noffsets-1)*2 ) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- FQ_IndexBinned[" << nm.c_str()
            << "]::readOld: the number of keys (" << nkeys << ")"
            << " does not match to the number of offsets ("
            << noffsets << ")"
            << "\n\t #key should be equal to (#offset-1)*2";
        nobs = 0;
        return false;
    }

    clear();
    nobs = nkeys / 2;
    bounds.resize(nobs);
    maxval.resize(nobs);
    minval.resize(nobs);

    // read in the keys (store the reading result in a temporary buf)
    switch (m_type) {
    case ibis::INT: {
        ibis::array_t<int32_t> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    case ibis::DOUBLE: {
        ibis::array_t<double> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    case ibis::BYTE: {
        ibis::array_t<char> buf(nkeys);
        berr = indexFile.getBitmapKeys(nm, buf.begin(), mpi_idx);
        if (berr) {
            for (unsigned i = 0; i < nobs; ++ i) {
                maxval[i] = buf[i+nobs];
                minval[i] = buf[i];
                if (i > 0)
                    bounds[i-1] =
                        ibis::util::compactValue(maxval[i-1], buf[i]);
            }
            bounds[nobs-1] = DBL_MAX;
        }
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str()
            << "]::readOld: does not yet support data type "
            << ibis::TYPESTRING[(int)m_type];
        berr = false;
        break;}
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str()
            << "]::readOld: failed to read BitmapKeys from file \""
            << indexFile.getFileName() << "\"";
        clear();
        return false;
    }

    // read in the offsets
    FQ::DataType type = indexFile.getBitmapOffsetType(nm);
    switch(type) {
    case FQ::FQT_INT: {
        offset64.clear();
        offset32.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset32.begin(), mpi_idx);
        break;}
    case FQ::FQT_LONG: {
        offset32.clear();
        offset64.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset64.begin(), mpi_idx);
        break;}
    case FQ::FQT_UINT: {
        // there is a potential overflow problem because unsigned int is
        // used in previous fastbit version
        offset64.clear();
        offset32.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset32.begin(), mpi_idx);
        break;}
    case FQ::FQT_ULONG: {
        // there is a potential overflow problem because unsigned int is
        // used in previous fastbit version
        offset32.clear();
        offset64.resize(nobs+1);
        berr = indexFile.getBitmapOffsets(nm, offset64.begin(), mpi_idx);
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::readOld: "
            << "unsupported data type for BitmapOffsets";
        berr = false;
        break;}
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str()
            << "]::readOld: failed to read BitmapOffsets from file \""
            << indexFile.getFileName() << "\"";
        clear();
        return false;
    }

    // activate the first bitmap so that the size information is available
    bits.resize(nobs);
    for (unsigned i = 0; i < nobs; ++ i)
        bits[i] = 0;

    //activate(0);
    std::vector<uint64_t> dims = varSpace.getCounts();
    if (! dims.empty()) {
        nrows = 1;
        for (unsigned i = 0; i < dims.size(); ++ i)
            nrows *= dims[i];
    }
    else {
        nrows = 0;
    }
    fname = ibis::util::strnewdup(nm.c_str());
    status = 0;

    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexBinned::readOld\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;
    return true;
} // FQ_IndexBinned::readOld

/// Write the index to an active HDF5 file.
bool FQ_IndexBinned::write(ArrayIODriver& outIndexFile) const
{
    ibis::horometer timer;
    ibis::horometer setupTimer;
    ibis::horometer metadataTimer;
    ibis::horometer bitmapTimer;
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    ibis::horometer initIdleTimer;
    ibis::horometer bitmapIdleTimer;
    ibis::horometer metadataIdleTimer;
    ibis::horometer setupIdleTimer;

    initIdleTimer.start();
    MPI_Barrier(mpi_comm);
    initIdleTimer.stop();
#endif

    timer.start();
    const size_t nobs = bits.size();
#ifdef FQ_NOMPI
    if (nobs == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "no bitmaps need to be written";
        return true;
    }
#endif

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    bool berr = true;

    if (status <= 0) {
        if (&outIndexFile == &indexFile) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
                << "no need to write to the same file \""
                << indexFile.getFileName() << "\"";
            return true;
        }
        activate();
    }

    unsigned int bitmapLen = 0;
    if (nobs > 0) {
        if (offset64.size() <= nobs) {
            offset64.resize(nobs+1);
            offset64[0] = 0;
            for (unsigned i = 0; i < nobs; ++ i) {
                offset64[i+1] = offset64[i];
                if (bits[i] != 0) {
                    const ibis::bitvector::word_t w = bits[i]->getSerialSize();
                    if (w > 1) offset64[i+1] += (w >> 2);
                }
            }
        }
        bitmapLen = offset64[nobs];
    }

    setupTimer.start();
#ifndef FQ_NOMPI
    if (mpi_idx < mpi_size) {
        // in the first iteration, create offset dataset collectively
        berr = outIndexFile.createOffsetTable(nm, mpi_max_iter);
        if (! berr) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
                << "failed to create offset dataset";
            return false;
        }
        LOGGER(ibis::gVerbose > 2)
            << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "created offset dataset";
    }

    berr = outIndexFile.setBitmapLength(nm, bitmapLen, mpi_iter);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to set bitmapLen";
        return false;
    }
#endif

    berr = outIndexFile.createBitmap(nm, bitmapLen, mpi_iter, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to create bitmap dataset with bitmap length "
            << bitmapLen;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
        << "bitmap dataset is created with bitmap length " << bitmapLen;

    unsigned int bitmapOffsetLen = 0;
    if (nobs > 0) bitmapOffsetLen = nobs+1;

#ifndef FQ_NOMPI
    berr = outIndexFile.setBitmapOffsetLength(nm, bitmapOffsetLen, mpi_iter);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to set bitmapOffsetLen";
        return false;
    }
#endif

    berr = outIndexFile.createBitmapOffsets(nm, bitmapOffsetLen,
                                            FQ::FQT_LONG, mpi_iter, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str()
            << "]::write: failed to create bitmapOffset dataset with length "
            << bitmapOffsetLen;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
        << "bitmapOffset dataset is created with bitmap offset length "
        << bitmapOffsetLen;

    unsigned int nkeys = 2*nobs;
#ifndef FQ_NOMPI
    berr = outIndexFile.setBitmapKeyLength(nm, nkeys, mpi_iter);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to set bitmapKeyLength";
        return false;
    }
#endif

    berr = outIndexFile.createBitmapKeys(nm, nkeys, varInfo.getType(),
                                         mpi_iter, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to create bitmapKey dataset with bitmap key length "
            << nkeys;
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
        << "bitmapKey dataset is created with bitmap key length " << nkeys;

#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    setupIdleTimer.start();
    MPI_Barrier(mpi_comm);
    setupIdleTimer.stop();
#endif
    setupTimer.stop();

    metadataTimer.start();
    // write the bitmap offsets
    if (nobs > 0) {
        berr = outIndexFile.setBitmapOffsets
            (nm, const_cast<int64_t*>(offset64.begin()),
             nobs+1, mpi_idx);
    } else {
        // to make sure it can be handled collectively with other MPI tasks
        berr = outIndexFile.setBitmapOffsets
            (nm, NULL, 0, mpi_idx);
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to write the bitmap offset to file \""
            << outIndexFile.getFileName() << "\"";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexBinned[" << nm.c_str() << "]::write: "
        << "recorded bitmapOffsets as destined for file \""
        << outIndexFile.getFileName() << "\"";

    // collectively/independently write the keys, need return the min/max
    // values to their original type
    bool berr1 = true;
    bool berr2 = true;
    switch(m_type) {
    case ibis::DOUBLE: {
        ibis::array_t<double> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = minval[i];
            buf[i+nobs] = maxval[i];
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<double> range(2);
        range[0] = getMin();
        range[1] = getMax();
        berr1 = outIndexFile.setActualRange
            (nm, range.begin(), FQ::FQT_DOUBLE);
        berr2 = outIndexFile.setExpectedRange
            (nm, range.begin(), FQ::FQT_DOUBLE);
#endif
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = static_cast<float>(minval[i]);
            buf[i+nobs] = static_cast<float>(maxval[i]);
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<float> range(2);
        range[0] = static_cast<float>(getMin());
        range[1] = static_cast<float>(getMax());
        berr1 = outIndexFile.setActualRange(nm, range.begin(), FQ::FQT_FLOAT);
        berr2 = outIndexFile.setExpectedRange(nm, range.begin(), FQ::FQT_FLOAT);
#endif
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = static_cast<int64_t>(minval[i]);
            buf[i+nobs] = static_cast<int64_t>(maxval[i]);
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<int64_t> range(2);
        range[0] = static_cast<int64_t>(getMin());
        range[1] = static_cast<int64_t>(getMax());
        berr1 = outIndexFile.setActualRange(nm, range.begin(), FQ::FQT_LONG);
        berr2 = outIndexFile.setExpectedRange(nm, range.begin(), FQ::FQT_LONG);
#endif
        break;}
    case ibis::INT: {
        ibis::array_t<int32_t> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = static_cast<int32_t>(minval[i]);
            buf[i+nobs] = static_cast<int32_t>(maxval[i]);
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<int32_t> range(2);
        range[0] = static_cast<int32_t>(getMin());
        range[1] = static_cast<int32_t>(getMax());
        berr1 = outIndexFile.setActualRange(nm, range.begin(), FQ::FQT_INT);
        berr2 = outIndexFile.setExpectedRange(nm, range.begin(), FQ::FQT_INT);
#endif
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = static_cast<uint64_t>(minval[i]);
            buf[i+nobs] = static_cast<uint64_t>(maxval[i]);
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<uint64_t> range(2);
        range[0] = static_cast<uint64_t>(getMin());
        range[1] = static_cast<uint64_t>(getMax());
        berr1 = outIndexFile.setActualRange(nm, range.begin(), FQ::FQT_ULONG);
        berr2 = outIndexFile.setExpectedRange(nm, range.begin(), FQ::FQT_ULONG);
#endif
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = static_cast<uint32_t>(minval[i]);
            buf[i+nobs] = static_cast<uint32_t>(maxval[i]);
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<uint32_t> range(2);
        range[0] = static_cast<uint32_t>(getMin());
        range[1] = static_cast<uint32_t>(getMax());
        berr1 = outIndexFile.setActualRange(nm, range.begin(), FQ::FQT_UINT);
        berr2 = outIndexFile.setExpectedRange(nm, range.begin(), FQ::FQT_UINT);
#endif
        break;}
    case ibis::BYTE: {
        ibis::array_t<char> buf(nkeys);
        for (unsigned i = 0; i < nobs; ++ i) {
            buf[i] = static_cast<char>(minval[i]);
            buf[i+nobs] = static_cast<char>(maxval[i]);
        }
        berr = outIndexFile.setBitmapKeys(nm, buf.begin(), nkeys, mpi_idx);
#ifdef FQ_NOMPI
        ibis::array_t<char> range(2);
        range[0] = static_cast<char>(getMin());
        range[1] = static_cast<char>(getMax());
        berr1 = outIndexFile.setActualRange(nm, range.begin(), FQ::FQT_UINT);
        berr2 = outIndexFile.setExpectedRange(nm, range.begin(), FQ::FQT_UINT);
#endif
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str()
            << "]::write: does not yet support data type "
            << ibis::TYPESTRING[(int)m_type];
        return false;}
    }
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to write the bitmap keys array[" << nobs <<"]";
        return false;
    }
    if (! berr1) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to write the actual key range";
        return false;
    }
    if (! berr2) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]::write: "
            << "failed to write the expected key range";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "FQ_IndexBinned[" << nm.c_str() << "]::write "
        << "recorded bitmapKeys as destined for file \""
        << outIndexFile.getFileName() << "\"";

#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    metadataIdleTimer.resume();
    MPI_Barrier(mpi_comm);
    metadataIdleTimer.stop();
#endif
    metadataTimer.stop();
    // go through the bitmaps and write them by chunk
    bitmapTimer.start();
#ifdef FQ_PACK_BITMAP
    ibis::array_t<ibis::bitvector::word_t> chunk;
    uint64_t head = 0;
    for (unsigned i = 0; i < nobs; ++ i) {
        ibis::array_t<ibis::bitvector::word_t> buf;
        bits[i]->write(buf);
        uint64_t pos = offset64[i] - head;
        if (offset64[i+1]-head > FQ_CHUNK_SIZE) {
            // break the bitmap into two chunks
            uint64_t len = FQ_CHUNK_SIZE-(offset64[i]-head);
            chunk.insert(chunk.begin()+pos, buf.begin(), buf.begin()+len);
            // write the current chunk
            if (! outIndexFile.writeBitmap
                (nm, head, head+chunk.size(),
                 reinterpret_cast<uint32_t*>(chunk.begin()), mpi_idx)) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_IndexBinned[" << nm.c_str()
                    << "]::write: failed to write bitvector chunk";
                return false;
            }
            // restart with a new chunk
            chunk.clear();
            chunk.insert(chunk.begin(), buf.begin()+len, buf.end());
            head += FQ_CHUNK_SIZE;
        } else {
            uint64_t len = offset64[i+1] - offset64[i];
            chunk.insert(chunk.begin()+pos, buf.begin(), buf.begin()+len);
        }
    }
    if (nobs > 0) {
        if (! outIndexFile.writeBitmap
            (nm, head, head+chunk.size(),
             reinterpret_cast<uint32_t*>(chunk.begin()), mpi_idx)) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- FQ_IndexBinned[" << nm.c_str()
                << "]::write: failed to write bitvector chunk";
            return false;
        }
    }
#else
    // go through the bitmaps to write each one separately
    for (unsigned i = 0; i < nobs; ++ i) {
        if (offset64[i] < offset64[i+1]) {
            ibis::array_t<ibis::bitvector::word_t> buf;
            bits[i]->write(buf);
            berr = outIndexFile.writeBitmap
                (nm, offset64[i], offset64[i+1],
                 reinterpret_cast<uint32_t*>(buf.begin()), mpi_idx);
            if (! berr) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_IndexBinned[" << nm.c_str()
                    << "]::write: failed to write bitvector "
                    << i << " (" << offset64[i+1]-offset64[i] << " words)";
                return false;
                // is there a way to rollback the write operations?
            }
        }
    }
#endif

#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
    bitmapIdleTimer.start();
    MPI_Barrier(mpi_comm);
    bitmapIdleTimer.stop();
#endif
    bitmapTimer.stop();
    timer.stop();

    LOGGER(FastQuery::reportTiming())
#if (! defined(FQ_NOMPI)) && defined(FQ_PERFORMANCE_TEST)
        << "\nStatistic\tFQ_IndexBinned::write:init_idle\t"
        << initIdleTimer.CPUTime() << "\t" << initIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexBinned::write:setup_idle\t"
        << setupIdleTimer.CPUTime() << "\t" << setupIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexBinned::write:metadata_idle\t"
        << metadataIdleTimer.CPUTime() << "\t"
        << metadataIdleTimer.realTime()
        << "\nStatistic\tFQ_IndexBinned::write:bitmap_idle\t"
        << bitmapIdleTimer.CPUTime() << "\t" << bitmapIdleTimer.realTime()
#endif
        << "\nStatistic\tFQ_IndexBinned::write:setup\t"
        << setupTimer.CPUTime() << "\t" << setupTimer.realTime()
        << "\nStatistic\tFQ_IndexBinned::write:metadata\t"
        << metadataTimer.CPUTime() << "\t" << metadataTimer.realTime()
        << "\t" << nobs
        << "\nStatistic\tFQ_IndexBinned::write:bitmap\t"
        << bitmapTimer.CPUTime() << "\t" << bitmapTimer.realTime() << "\t"
        << bitmapLen
        << "\nStatistic\tFQ_IndexBinned::write\t"
        << timer.CPUTime() << "\t" << timer.realTime() ;
    return true;
} // FQ_IndexBinned::write

/// Activate all bitmaps at once by reading all of them into memory.
void FQ_IndexBinned::activate() const {
    ibis::horometer timer;
    timer.start();

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    // assumes that bounds, minval, maxval, bits and offsets have been
    // initialized properly.
    if (nobs == 0 || nobs != bits.size() ||
        (offset32.size() <= nobs && offset64.size() <= nobs)) return;

    bool berr;
    if ((bits.empty() || bits[0] == 0) && str == 0) {
        // read all bitmaps into memory and store the bytes in str
        const size_t nwords = (offset64.size() > nobs ? offset64[nobs] :
                               (int64_t) offset32[nobs]);
        str = new
            ibis::fileManager::storage(nwords*sizeof(ibis::bitvector::word_t));
        berr = indexFile.readBitmap(nm, 0, nwords,
                                    (uint32_t*)(str->begin()), mpi_idx);
        if (berr == false) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexBinned[" << nm.c_str() << "]::activate: "
                << "unable to read all bitmaps at once";
            delete str;
            str = 0;
        }
    }

    if (str) {
        if (offset64.size() > nobs) {
            for (unsigned i = 0; i < nobs; ++ i) {
                if (bits[i] == 0 && offset64[i+1] > offset64[i]) {
                    ibis::array_t<ibis::bitvector::word_t>
                        buf(str, offset64[i]*sizeof(ibis::bitvector::word_t),
                            offset64[i+1]*sizeof(ibis::bitvector::word_t));
                    bits[i] = new ibis::bitvector(buf);
                }
            }
        }
        else {
            for (unsigned i = 0; i < nobs; ++ i) {
                if (bits[i] == 0 && offset32[i+1] > offset32[i]) {
                    ibis::array_t<ibis::bitvector::word_t>
                        buf(str, offset32[i]*sizeof(ibis::bitvector::word_t),
                            offset32[i+1]*sizeof(ibis::bitvector::word_t));
                    bits[i] = new ibis::bitvector(buf);
                }
            }
        }
    }
    else { // need to read one bitmap at a time
        for (unsigned i = 0; i < nobs; ++ i) {
            uint64_t start, end;
            if (offset64.size() > nobs) {
                start = offset64[i];
                end = offset64[i+1];
            }
            else {
                start = offset32[i];
                end = offset32[i+1];
            }
            if (bits[i] == 0 && end > start) {
                ibis::array_t<ibis::bitvector::word_t> buf(end-start);
                berr = indexFile.readBitmap(nm, start, end,
                                            (uint32_t*)(buf.begin()), mpi_idx);
                if (berr) {
                    bits[i] = new ibis::bitvector(buf);
                }
                else {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- FQ_IndexBinned[" << nm.c_str()
                        << "]::activate: failed to read bitmap " << i
                        << " (offset " << start << ", size " << end-start
                        << ")";
                }
            }
        }
    }
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexBinned::active()\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;
} // FQ_IndexBinned::activate

/// Activate the ith bitmap
void FQ_IndexBinned::activate(uint32_t i) const {
    ibis::horometer timer;
    timer.start();

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();

    if (i >= nobs) return;
    if (bits[i] != 0) return;
    if (bits.size() != nobs) return;
    if (offset64.size() <= nobs && offset32.size() <= nobs) return;
    std::string evt = "FQ_IndexBinned::activate";
    if (ibis::gVerbose > 2) {
        std::ostringstream oss;
        oss << "(" << i << ")";
        evt += oss.str();
    }
    if (offset64.size() > nobs ? offset64[i] >= offset64[i+1] :
        offset32[i] >= offset32[i+1]) {
        return;
    }

    bool berr;
    const int64_t sz1 = (offset64.size()>nobs ? (offset64[1]-offset64[0]) :
                         (int64_t)(offset32[1]-offset32[0]));
    const int64_t sza = (offset64.size()>nobs ? (offset64[nobs]-offset64[0]) :
                         (int64_t)(offset32[nobs]-offset32[0]));
    if (i == 0 && (nobs == 1U || sz1*5/4 >= sza
                   || (int64_t)(sz1*std::log((double)nobs)) >= sza)) {
        // if the first bitmap takes up a majority of the total bytes, read
        // them all
        const size_t nwords = (offset64.size() > nobs ? offset64[nobs] :
                               (int64_t)offset32[nobs]);
        str = new ibis::fileManager::storage
            (nwords*sizeof(ibis::bitvector::word_t));
        berr = indexFile.readBitmap(nm, 0, nwords,
                                    (uint32_t*)(str->begin()), mpi_idx);
        if (berr == false) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexBinned[" << nm.c_str() << "]"
                << "::activate" << "(" << i << "): "
                << "unable to read all bitmaps at once";
            delete str;
            str = 0;
        }
    }

    if (offset64.size() > nobs) {
        if (str) { // already in memory
            ibis::array_t<ibis::bitvector::word_t> buf
                (str, offset64[i]*sizeof(ibis::bitvector::word_t),
                 offset64[i+1]*sizeof(ibis::bitvector::word_t));
            bits[i] = new ibis::bitvector(buf);
        }
        else { // read the specified bitmap based on the offsets array
            ibis::array_t<ibis::bitvector::word_t>
                buf(offset64[i+1]-offset64[i]);
            berr = indexFile.readBitmap(nm, offset64[i], offset64[i+1],
                                        (uint32_t*)(buf.begin()), mpi_idx);
            if (berr) {
                bits[i] = new ibis::bitvector(buf);
            }
            else {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]"
                    << "::activate" << "(" << i << "):"
                    << " failed to read bitmap " << i
                    << " (offset " << offset64[i] << ", size "
                    << offset64[i+1]-offset64[i] << ")";
            }
        }
    }
    else if (str) { // already in memory
        ibis::array_t<ibis::bitvector::word_t> buf
            (str, offset32[i]*sizeof(ibis::bitvector::word_t),
             offset32[i+1]*sizeof(ibis::bitvector::word_t));
        bits[i] = new ibis::bitvector(buf);
    }
    else { // read the specified bitmap based on the offsets array
        ibis::array_t<ibis::bitvector::word_t> buf(offset32[i+1]-offset32[i]);
        berr = indexFile.readBitmap(nm, offset32[i], offset32[i+1],
                                    (uint32_t*)(buf.begin()), mpi_idx);
        if (berr) {
            bits[i] = new ibis::bitvector(buf);
        }
        else {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]"
                << "::activate" << "(" << i << "):"
                << " failed to read bitmap " << i
                << " (offset " << offset32[i] << ", size "
                << offset32[i+1]-offset32[i] << ")";
        }
    }
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexBinned::active(" << i << ")\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;
} // FQ_IndexBinned::activate

/// Activate the ith through jth bitmap
void FQ_IndexBinned::activate(uint32_t i, uint32_t j) const {
    ibis::horometer timer;
    timer.start();

    ArrayIODriver& indexFile =
        reinterpret_cast<const FQ_Variable*>(col)->getIndexFile();
    if (i >= nobs || i >= j) return; // empty range
    if (bits.size() != nobs) return;
    if (offset64.size() <= nobs && offset32.size() <= nobs) return;
    std::string evt = "FQ_IndexBinned::activate";
    if (ibis::gVerbose > 2) {
        std::ostringstream oss;
        oss << "(" << i << ", " << j << ")";
        evt += oss.str();
    }
    ibis::util::timer mytimer(evt.c_str(), 3);

    bool berr;
    const int64_t sz1 = (offset64.size()>nobs ? (offset64[1]-offset64[0]) :
                         (int64_t)(offset32[1]-offset32[0]));
    const int64_t sza = (offset64.size()>nobs ? (offset64[nobs]-offset64[0]) :
                         (int64_t)(offset32[nobs]-offset32[0]));
    if (i == 0 && (nobs == 1U || sz1*5/4 >= sza
                   || (int64_t)(sz1*std::log((double)nobs)) >= sza)) {
        // if the first bitmap takes up a majority of the total bytes, read
        // them all
        const size_t nwords = (offset64.size() > nobs ? offset64[nobs] :
                               (int64_t)offset32[nobs]);
        str = new ibis::fileManager::storage
            (nwords*sizeof(ibis::bitvector::word_t));
        berr = indexFile.readBitmap(nm, 0, nwords,
                                    (uint32_t*)(str->begin()), mpi_idx);
        if (berr == false) {
            LOGGER(ibis::gVerbose > 2)
                << "FQ_IndexBinned[" << nm.c_str() << "]"
                << "::activate" << "(" << i << "," << j << "):"
                << " unable to read all bitmaps at once";
            delete str;
            str = 0;
        }
    }

    if (offset64.size() > nobs) {
        if (str) { // already in memory
            while (i < j) {
                if (bits[i] == 0 && offset64[i+1] > offset64[i]) {
                    ibis::array_t<ibis::bitvector::word_t> buf
                        (str, offset64[i]*sizeof(ibis::bitvector::word_t),
                         offset64[i+1]*sizeof(ibis::bitvector::word_t));
                    bits[i] = new ibis::bitvector(buf);
                }
                ++ i;
            }
        }
        else {
            while (i < j) {
                // skip to next empty bit vector
                while (i < j && bits[i] != 0)
                    ++ i;
                // the last bitvector to activate. can not be larger
                // than j
                unsigned aj = (i<j ? i + 1 : j);
                while (aj < j && bits[aj] == 0)
                    ++ aj;
                if (offset64[aj] > offset64[i]) {
                    // read bitmaps into memory in one shot
                    const unsigned start = offset64[i];
                    ibis::array_t<ibis::bitvector::word_t>
                        buf(offset64[aj]-start);
                    berr = indexFile.readBitmap
                        (nm, offset64[i], offset64[aj],
                         (uint32_t*)(buf.begin()), mpi_idx);
                    if (berr) {
                        while (i < aj) {
                            if (bits[i] == 0 && offset64[i+1] > offset64[i]) {
                                ibis::array_t<ibis::bitvector::word_t>
                                    tmp(buf, offset64[i]-start,
                                        offset64[i+1]-start);
                                bits[i] = new ibis::bitvector(tmp);
                            }
                            ++ i;
                        }
                    }
                    else {
                        LOGGER(ibis::gVerbose > 0)
                            << "Warning -- FQ_IndexBinned[" << nm.c_str()
                            << "]::activate" << "(" << i << "," << j << "):"
                            << " failed to read bitmaps("
                            << offset64[i] << ", " << offset64[aj] << ")";
                    }
                    i = aj;
                } // if (offset64[aj] > offset64[i])
            } // while (i < j)
        }
    }
    else if (str) { // already in memory
        while (i < j) {
            if (bits[i] == 0 && offset32[i+1] > offset32[i]) {
                ibis::array_t<ibis::bitvector::word_t> buf
                    (str, offset32[i]*sizeof(ibis::bitvector::word_t),
                     offset32[i+1]*sizeof(ibis::bitvector::word_t));
                bits[i] = new ibis::bitvector(buf);
            }
            ++ i;
        }
    }
    else {
        while (i < j) {
            // skip to next empty bit vector
            while (i < j && bits[i] != 0)
                ++ i;
            // the last bitvector to activate. can not be larger
            // than j
            unsigned aj = (i<j ? i + 1 : j);
            while (aj < j && bits[aj] == 0)
                ++ aj;
            if (offset32[aj] > offset32[i]) {
                // read bitmaps into memory in one shot
                const unsigned start = offset32[i];
                ibis::array_t<ibis::bitvector::word_t> buf(offset32[aj]-start);
                berr = indexFile.readBitmap
                    (nm, offset32[i], offset32[aj],
                     (uint32_t*)(buf.begin()), mpi_idx);
                if (berr) {
                    while (i < aj) {
                        if (bits[i] == 0 && offset32[i+1] > offset32[i]) {
                            ibis::array_t<ibis::bitvector::word_t>
                                tmp(buf, offset32[i]-start,
                                    offset32[i+1]-start);
                            bits[i] = new ibis::bitvector(tmp);
                        }
                        ++ i;
                    }
                }
                else {
                    LOGGER(ibis::gVerbose > 0)
                        << "Warning -- FQ_IndexBinned[" << nm.c_str() << "]"
                        << "::activate" << "(" << i << "," << j << "):"
                        << " failed to read bitmaps("
                        << offset32[i] << ", " << offset32[aj] << ")";
                }
                i = aj;
            } // if (offset32[aj] > offset32[i])
        } // while (i < j)
    }
    timer.stop();
    LOGGER(FastQuery::reportTiming())
        << "Statistic\tFQ_IndexBinned::active("<< i << ", " << j << ")\t"
        << timer.CPUTime() << "\t" << timer.realTime() << "\t" << mpi_rank;
} // FQ_IndexBinned::activate

