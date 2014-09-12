#ifndef _FQ_COLUMN_H
#define _FQ_COLUMN_H

#include <column.h>     // ibis::column
#include <algorithm>
#include <sstream>      // std::ostringstream
#include <typeinfo>     // typeid
#include <vector>       // std::vector

#include "fqPart.h"     // ibis::column
#include "fqVar.h"
/**
   The class FQ_variable is a thin wrapper on ibis::column.  
   It maps each variable to a column of FastBit data partition.
*/
class FQ_Variable : public ibis::column {
public:
    FQ_Variable(const FQ_Part *tbl, const std::string &variableName,
		const std::string &variablePath,
                const VarInfo &varInfo_, const VarSpace &varSpace_,
                const ArrayIODriver &dataFile_,
		const ArrayIODriver &indexFile_);

    virtual ~FQ_Variable() {};

    bool isValid (const std::string &func) const;

    virtual int getValuesArray(void* arr) const;
    template <typename E>
    int getPointValues(ibis::array_t<E>& arr,
                       const std::vector<uint64_t>& coords) const;

    template <typename T> ibis::array_t<T>*
    selectData(const ibis::bitvector& mask, 
	       ibis::array_t<T>* array) const;
    virtual ibis::array_t<double>*
    selectDoubles(const ibis::bitvector& mask) const;
    virtual ibis::array_t<float>*
    selectFloats(const ibis::bitvector& mask) const;
    virtual ibis::array_t<int32_t>*
    selectInts(const ibis::bitvector& mask) const;
    virtual ibis::array_t<int64_t>*
    selectLongs(const ibis::bitvector& mask) const;
    virtual ibis::array_t<signed char>*
    selectBytes(const ibis::bitvector& mask) const;
    virtual ibis::array_t<uint32_t>*
    selectUInts(const ibis::bitvector& mask) const;
    virtual ibis::array_t<uint64_t>*
    selectULongs(const ibis::bitvector& mask) const;

    virtual void loadIndex(const char*, int) const throw ();
    virtual long indexSize() const;

#ifdef FQ_NOMPI
    // these functions are not implemented under MPI mode
    virtual double getActualMin() const;
    virtual double getActualMax() const;
#endif

    ArrayIODriver& getArrayIODriver() const {return this->dataFile;}
    ArrayIODriver& getIndexFile() const {return this->indexFile;}
    VarInfo& getVarInfo() const {return this->varInfo;}
    VarSpace& getVarSpace() const {return this->varSpace;}
    bool isAll() const {return useAll;}

protected:
    using ibis::column::searchSorted;
    /// Resolve a continuous range condition on a sorted column.
    virtual int searchSorted(const ibis::qContinuousRange&,
                             ibis::bitvector&) const;
    /// Resolve a discrete range condition on a sorted column.
    virtual int searchSorted(const ibis::qDiscreteRange&,
                             ibis::bitvector&) const;

private:
    // member variables
    bool useAll;
    ArrayIODriver &dataFile;
    ArrayIODriver &indexFile;
    VarInfo &varInfo;
    VarSpace &varSpace;
    bool valid;

    FQ_Variable(); // no default constructor

    /// Attempt to read a bitmap index from the HDF5.  If successful,
    /// return the index read, otherwise return 0 (null pointer).
    ibis::index* readIndex() const;
}; // class FQ_Variable
#endif
