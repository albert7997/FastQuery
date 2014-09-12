#ifndef _FASTQUERY_INDEX_H
#define _FASTQUERY_INDEX_H

#include "arrayIODriver.h"  // the underline file model operations
#include "fqColumn.h"

#include <index.h>      // FastBit abstract index class
#include <irelic.h>     // FastBit unbinned index classes
#include <ibin.h>       // FastBit binned index classes

#include <vector>
#include <string>

/**
   This file defines the data structure of the index classes specialized
   for the HDF5 files.  They are required for reading from and writing to
   HDF5 files.  Currently, only two specializations are implemented,
   FQ_IndexUnbinned and FQ_IndexBinned.  The first of which inherents
   from ibis::relic, the basic bitmap index that indexes every distinct
   value.  The second version implements the binned version of the bitmap
   index, which potentially requires less storage but may require access of
   the raw data to answer some queries.
*/
class FQ_IndexUnbinned : public ibis::relic {
public:
    // Build index in memory if readOnly is false.
    FQ_IndexUnbinned(const FQ_Variable* c, bool readOnly);
    ~FQ_IndexUnbinned() {};

    /// Returns negative number(-1) if index is not built
    /// Returns 0 if index is read from file
    /// Returns postive number(1) if index is built in memory
    int getStatus() const {return status;}

    using ibis::relic::write;
    bool write(ArrayIODriver& outIndexFile) const;

protected:
    virtual void activate() const;
    virtual void activate(uint32_t i) const;
    virtual void activate(uint32_t i, uint32_t j) const;

    virtual void clear() {ibis::relic::clear(); status=-1;}

private:
    int status; //-1: no index, 0: read old index, 1: built new index
    ibis::TYPE_T m_type;
    std::string nm;
    VarInfo &varInfo;
    VarSpace &varSpace;

    int mpi_idx;
    int mpi_iter;
    int mpi_max_iter;
    int mpi_size;
    int mpi_rank;
#ifndef FQ_NOMPI
    MPI_Comm mpi_comm;
#endif
    /// Load the content of the index from file
    bool readOld(ArrayIODriver& indexFile);

    /// Create new index
    bool buildNew();

    FQ_IndexUnbinned(); // no default constructor
    FQ_IndexUnbinned(const FQ_IndexUnbinned&); // no copy constructor
    const FQ_IndexUnbinned& operator=(const FQ_IndexUnbinned&);
}; // class FQ_IndexUnbinned


/**
   Bitmap indices with bins.
*/
class FQ_IndexBinned : public ibis::bin {
public:
    // Build index in memory if readOnly is false.
    FQ_IndexBinned(const FQ_Variable* c, const char *binning, bool readOnly);
    ~FQ_IndexBinned() {};

    /// Returns negative number(-1) if index is not built
    /// Returns 0 if index is read from file
    /// Returns postive number(1) if index is built in memory
    int getStatus() const {return status;}

    using ibis::bin::write;
    bool write(ArrayIODriver& outIndexFile) const;

protected:
    virtual void activate() const;
    virtual void activate(uint32_t i) const;
    virtual void activate(uint32_t i, uint32_t j) const;

    virtual void clear() {ibis::bin::clear(); status=-1;}

private:
    int status;
    ibis::TYPE_T m_type;
    std::string nm;
    VarInfo &varInfo;
    VarSpace &varSpace;

    int mpi_idx;
    int mpi_iter;
    int mpi_max_iter;
    int mpi_size;
    int mpi_rank;
#ifndef FQ_NOMPI
    MPI_Comm mpi_comm;
#endif
    
    /// Load the content of the index from file
    bool readOld(ArrayIODriver &indexFile);
    // vidcina added it 
    void* Construct(void* arg);
    /// Create new index
    bool buildNew();

    FQ_IndexBinned(); // no default constructor
    FQ_IndexBinned(const FQ_IndexBinned&); // no copy constructor
    const FQ_IndexBinned& operator=(const FQ_IndexBinned&);
}; // class FQ_IndexBinned
#endif
