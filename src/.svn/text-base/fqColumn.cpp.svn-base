#include "fqColumn.h"
#include "fqIndex.h"    
#include <memory>	// std::auto_ptr

/*******************************************
 * Constructor
 ********************************************/
FQ_Variable::FQ_Variable(const FQ_Part *tbl,
			 const std::string &variableName, 
			 const std::string &variablePath, 
			 const VarInfo &varInfo_,
			 const VarSpace &varSpace_,
			 const ArrayIODriver &dataFile_,
			 const ArrayIODriver &indexFile_)
    : ibis::column(tbl, ibis::OID, variableName.c_str(), variablePath.c_str()), 
      varInfo(const_cast<VarInfo&>(varInfo_)), 
      varSpace(const_cast<VarSpace&>(varSpace_)),  
      dataFile(const_cast<ArrayIODriver&>(dataFile_)), 
      indexFile(const_cast<ArrayIODriver&>(indexFile_))
{
    useAll = true;
    if (varInfo.getSize() != varSpace.getSize()) {
	useAll = false;
    }
    FQ::DataType type = varInfo.getType();
    switch (type) {
    case FQ::FQT_FLOAT:
        m_type = ibis::FLOAT; break;
    case FQ::FQT_DOUBLE:
        m_type = ibis::DOUBLE; break;
    case FQ::FQT_BYTE:
        m_type = ibis::BYTE; break;
    case FQ::FQT_UBYTE:
        m_type = ibis::UBYTE; break;
    case FQ::FQT_SHORT:
        m_type = ibis::SHORT; break;
    case FQ::FQT_USHORT:
        m_type = ibis::USHORT; break;
    case FQ::FQT_INT:
        m_type = ibis::INT; break;
    case FQ::FQT_UINT:
        m_type = ibis::UINT; break;
    case FQ::FQT_LONG:
        m_type = ibis::LONG; break;
    case FQ::FQT_ULONG:
        m_type = ibis::ULONG; break;
    default:
	LOGGER(ibis::gVerbose > 0)
	    << "Warnning -- FQ_Variable(" << varInfo.getPath()
	    << "): does not yet support data type " << type;
	valid = false;
    }

    // find out what field is sorted, and if this is the varaible, set
    // m_sorted field to true
    std::string sorted_var = dataFile.getSortedFieldName();
    isSorted(strcmp(sorted_var.c_str(), varInfo.getPath().c_str())==0);
    valid = true;
} // FQ_Variable::FQ_Variable

/*******************************************
 * Public Functions
 ********************************************/

bool FQ_Variable::isValid (const std::string &func) const {
    if (! valid) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << func.c_str()
            << ": partition was not initialized correctly";
        return false;
    } else {
        return true;
    }
} // FQ_Variable::isValid

#ifdef FQ_NOMPI
double FQ_Variable::getActualMin() const{
    if (! isValid("FQ_Variable::getActualMin")) return 0;

    std::string path = varInfo.getPath();
    if (! useAll) {
    	path += varSpace.getText();
    }
	
    switch (m_type) {
    case ibis::BYTE: {
        char range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    case ibis::DOUBLE: {
        double range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    case ibis::FLOAT: {
        float range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    case ibis::INT: {
        int32_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    case ibis::LONG: {
        int64_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    case ibis::UINT: {
        uint32_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    case ibis::ULONG: {
        uint64_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[0];
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
	    << "Warning -- FQ_Variable::"
            << "getActualMin[" << path.c_str() <<"]:"
	    << " does not yet support column type "
            << ibis::TYPESTRING[(int)m_type];
        break;}
    }
    LOGGER(ibis::gVerbose > 0)
	<< "Warning -- FQ_Variable::"
        << "getActualMin[" << path.c_str() <<"]:"
	<< " failed to read actualRange";
    return 0.0;
} // FQ_Variable::getActualMin

double FQ_Variable::getActualMax() const{
    if (! isValid("FQ_Variable::getActualMax")) return 0;

    std::string path = varInfo.getPath();
    if (! useAll) {
    	path += varSpace.getText();
    }

    switch (m_type) {
    case ibis::BYTE: {
        char range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    case ibis::DOUBLE: {
        double range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    case ibis::FLOAT: {
        float range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    case ibis::INT: {
        int32_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    case ibis::LONG: {
        int64_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    case ibis::UINT: {
        uint32_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    case ibis::ULONG: {
        uint64_t range[2];
        bool berr = indexFile.getActualRange(path, (void*)range);
	if (berr) return range[1];
        break;}
    default: {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- FQ_Variable"
            << "::getActualMax: does not yet support column type "
            << ibis::TYPESTRING[(int)m_type];
        break;}
    }
    LOGGER(ibis::gVerbose > 0)
	<< "Warning -- FQ_Variable"
	<< "::getActualMax: failed to read actualRange";
    return 0.0;
} // FQ_Variable::getActualMax
#endif

/// Read all values into an array_t object.  The argument arr must be
/// an array_t<Type>* with the correct Type.  Currently, the 
/// supported types mapped to C++ elementary types as follows:
///
/// - Float: array_t<float>*
/// - Double: array_t<double>*
/// - Int64: array_t<int64_t>*
/// - UInt64: array_t<uint64_t>*
/// - Int: array_t<int32_t>*
/// - UInt: array_t<uint32_t>*
/// - Short: array_t<int16_t>*
/// - UShort: array_t<uint16_t>*
/// - Byte: array_t<signed char>*
/// - UByte: array_t<unsigned char>*
///
/// It returns 0 to indicate success and a negative number to indicate
/// error.
///
/// @note The array is resized to have enough space to store the bytes.  No
/// type checking is possible, this function relies on the caller to pass
/// in an array of the correct type.
int FQ_Variable::getValuesArray(void* arr) const {
    if (! isValid("FQ_Variable::getValuesArray")) return -3;

    std::string evt = "FQ_Variable::getValuesArray";
    if (ibis::gVerbose > 1) {
        std::ostringstream oss;
        oss << '(' << (thePart->name() ? thePart->name() : "?") << '.'
            << name() << ", " << ibis::TYPESTRING[(int)m_type]
            << ", nrows=" << thePart->nRows() << ')';
        evt += oss.str();
    }
    uint64_t nElements = varSpace.getSize();
    
    LOGGER(ibis::gVerbose > 0)
	<< "FQ_Variable::getValuesArray the nElements size is " << nElements ;

    if (nElements == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- " << "FQ_Variable::getValuesArray:" 
	    << " no data needs to be read from file";
	return 0;
    }

    bool berr = false;
    switch (m_type) {
    case ibis::BYTE: {
        ibis::array_t<signed char> &vals =
	    *static_cast<ibis::array_t<signed char>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    case ibis::FLOAT: {
        ibis::array_t<float> &vals =
	    *static_cast<ibis::array_t<float>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    case ibis::DOUBLE: {
        ibis::array_t<double> &vals =
	    *static_cast<ibis::array_t<double>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    case ibis::INT: {
        ibis::array_t<int32_t> &vals =
	    *static_cast<ibis::array_t<int32_t>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> &vals =
	    *static_cast<ibis::array_t<uint32_t>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> &vals =
	    *static_cast<ibis::array_t<int64_t>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> &vals =
	    *static_cast<ibis::array_t<uint64_t>*>(arr);
        vals.resize(nElements);
	if (useAll) {
    	    berr = dataFile.getData(varInfo.getPath(), vals.begin());
	} else {
    	    berr = dataFile.getArrayData
		(varInfo.getPath(),
		 varSpace.getOffsets(), varSpace.getCounts(), 
		 varSpace.getStrides(), vals.begin());
	}
        break;}
    default: {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- FQ_Variable["
            << (thePart ? thePart->name() : "?") << '.'
            << m_name << "]::getValuesArray does not yet support column type "
            << ibis::TYPESTRING[(int)m_type];
        return -3;}
    } // switch (m_type)

    if (berr == false) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "FQ_Variable::getValuesArray:" 
	    << " fail to read data from file";
	return -2;
    } else {
	return 0;
    }
} // FQ_Variable::getValuesArray

template <typename E>
int FQ_Variable::getPointValues(ibis::array_t<E>& arr,
                                const std::vector<uint64_t>& coords) const {
    if (! isValid("FQ_Variable::getPointValues")) return 0;

    std::string evt = "FQ_Variable::getPointValues";
    if (ibis::gVerbose > 1) {
        std::ostringstream oss;
        oss << '(' << (thePart->name() ? thePart->name() : "?") << '.'
            << name() << ", " << ibis::TYPESTRING[(int)m_type]
            << ", coords[" << coords.size() << "])";
        evt += oss.str();
    }
    arr.resize(coords.size() / varInfo.getNDims());

    LOGGER(ibis::gVerbose > 0)
        << "FQ_Variable::getPointValues the arr size is " << arr.size() ;

    bool ret;
    if (useAll) {
    	ret = dataFile.getPointData(varInfo.getPath(), coords,
				    static_cast<void*>(arr.begin()));
    } else {
	// convert to absolute coords based on space info.
	std::vector<uint64_t> absCoords;
	absCoords.resize(coords.size());
	std::vector<uint64_t> offsets = varSpace.getOffsets();
	std::vector<uint64_t> counts = varSpace.getCounts();
	std::vector<uint64_t> strides = varSpace.getStrides();
	unsigned int idx = 0;
	for (unsigned int i=0; i<arr.size(); i++) {
	    for (unsigned int j=0; j<varInfo.getNDims(); j++) {
		absCoords[idx] = coords[idx]* strides[j] + offsets[j];
		idx++;
	    }
	}
    	ret = dataFile.getPointData(varInfo.getPath(), absCoords,
				    static_cast<void*>(arr.begin()));
    }
    if (ret)
        return arr.size();
    else
        return -1;
} // FQ_Variable::getPointValues

long FQ_Variable::indexSize() const{
    if (! isValid("FQ_Variable::indexSize")) return 0;

    std::string path = varInfo.getPath();
    if (! useAll) {
    	path += varSpace.getText();
    }
    uint64_t size = 0;
    int mpi_idx = -1;
#ifndef FQ_NOMPI
    mpi_idx = varSpace.getMpiIdx();
#endif
    bool berr = indexFile.getBitmapLength(path, &size, mpi_idx);

    if (ibis::gVerbose > 3) {
        if (!berr || size<=0)
            logWarning("indexSize", "failed to determine the bitmap length "
                       "for variable %s in file %s", name(),
                       indexFile.getFileName().c_str());

        else
            logMessage("indexSize", "found bitmap length "
                       "for variable %s in file %s to be %lu",
                       name(), indexFile.getFileName().c_str(), size);
    }
    return size;
} // FQ_Variable::indexSize

/// This function is required to have two arguments in order to be a proper
/// virtual function, however, it does not make use any of them.
/// Therefore, they are never named.
void FQ_Variable::loadIndex(const char*, int) const throw () {
    if (! isValid("FQ_Variable::loadIndex")) return;

    writeLock lock(this, "loadIndex");
    if (idx == 0 && thePart->nRows() > 0) {
        try { // if an index is not available, create one
            if (ibis::gVerbose > 7) {
                std::string fname = indexFile.getFileName();
                logMessage("loadIndex", "loading an index from %s",
                           fname.c_str());
            }
            if (idx == 0) {
	    	// try to read index first
	    	bool readOnly = true;
		std::auto_ptr<FQ_IndexBinned>
		    binnedIdx(new FQ_IndexBinned(this, "", readOnly));
		if (binnedIdx.get() != 0 && binnedIdx->getStatus() >= 0 ) {
		    idx = binnedIdx.release();
		}
		else {
		    std::auto_ptr<FQ_IndexUnbinned>
			unbinnedIdx(new FQ_IndexUnbinned(this, readOnly));
		    if (unbinnedIdx.get() != 0 &&
			unbinnedIdx->getStatus() >= 0) {
			idx = unbinnedIdx.release();
		    }
#ifdef FQ_ALWAYS_INDEX
		    else { // allow a new index to be created
			readOnly = false;
			idx = new FQ_IndexUnbinned(this, readOnly);
		    }
#endif
		}
            }
	    if (idx == 0) {
		logWarning("loadIndex", "failed to create an index object");
	    }
            else if (ibis::gVerbose > 8) {
                ibis::util::logger lg;
                idx->print(lg());
            }
        }
        catch (const char *s) {
            logWarning("loadIndex", "ibis::index::ceate(%s) throw "
                       "the following exception\n%s", name(), s);
            idx = 0;
        }
        catch (const std::exception& e) {
            logWarning("loadIndex", "ibis::index::create(%s) failed "
                       "to create a new index -- %s", name(), e.what());
            idx = 0;
        }
        catch (...) {
            logWarning("loadIndex", "ibis::index::create(%s) failed "
                       "to create a new index -- unknown error", name());
            idx = 0;
        }
    }
} // FQ_Variable::loadIndex

int FQ_Variable::searchSorted(const ibis::qContinuousRange& rng,
			      ibis::bitvector& hits) const {
    if (! isValid("FQ_Variable::searchSorted")) return 0;

    int ierr;
    LOGGER(ibis::gVerbose >= 5)
        << "... entering FQ_Variable::searchSorted to resolve " << rng;
    switch (m_type) {
    case ibis::BYTE: {
        ibis::array_t<signed char> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
	/*
	  case ibis::UBYTE: {
	  ibis::array_t<unsigned char> vals;
	  ierr = getValuesArray(&vals);
	  if (ierr >= 0) {
	  ierr = ibis::column::searchSortedICC(vals, rng, hits);
	  }
	  break;}
	  case ibis::SHORT: {
	  ibis::array_t<int16_t> vals;
	  ierr = getValuesArray(&vals);
	  if (ierr >= 0) {
	  ierr = ibis::column::searchSortedICC(vals, rng, hits);
	  }
	  break;}
	  case ibis::USHORT: {
	  ibis::array_t<uint16_t> vals;
	  ierr = getValuesArray(&vals);
	  if (ierr >= 0) {
	  ierr = ibis::column::searchSortedICC(vals, rng, hits);
	  }
	  break;}
	*/
    case ibis::FLOAT: {
        ibis::array_t<float> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
    case ibis::DOUBLE: {
        ibis::array_t<double> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
    case ibis::INT: {
        ibis::array_t<int32_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICC(vals, rng, hits);
        }
        break;}
    default: {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- FQ_Variable["
            << (thePart ? thePart->name() : "?") << '.'
            << m_name << "]::searchSorted(" << rng
            << ") does not yet support column type "
            << ibis::TYPESTRING[(int)m_type];
        ierr = -5;
        break;}
    } // switch (m_type)
    return (ierr < 0 ? ierr : 0);
} // FQ_Variable::searchSorted

int FQ_Variable::searchSorted(const ibis::qDiscreteRange& rng,
			      ibis::bitvector& hits) const {
    if (! isValid("FQ_Variable::searchSorted")) return 0;

    int ierr;
    LOGGER(ibis::gVerbose >= 5)
        << "... entering FQ_Variable::searchSorted to resolve " << rng;
    switch (m_type) {
    case ibis::BYTE: {
        ibis::array_t<signed char> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
	/*
	  case ibis::UBYTE: {
	  ibis::array_t<unsigned char> vals;
	  ierr = getValuesArray(&vals);
	  if (ierr >= 0) {
	  ierr = ibis::column::searchSortedICD(vals, rng, hits);
	  }
	  break;}
	  case ibis::SHORT: {
	  ibis::array_t<int16_t> vals;
	  ierr = getValuesArray(&vals);
	  if (ierr >= 0) {
	  ierr = ibis::column::searchSortedICD(vals, rng, hits);
	  }
	  break;}
	  case ibis::USHORT: {
	  ibis::array_t<uint16_t> vals;
	  ierr = getValuesArray(&vals);
	  if (ierr >= 0) {
	  ierr = ibis::column::searchSortedICD(vals, rng, hits);
	  }
	  break;}
	*/
    case ibis::FLOAT: {
        ibis::array_t<float> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
    case ibis::DOUBLE: {
        ibis::array_t<double> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
    case ibis::INT: {
        ibis::array_t<int32_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
    case ibis::UINT: {
        ibis::array_t<uint32_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
    case ibis::LONG: {
        ibis::array_t<int64_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
    case ibis::ULONG: {
        ibis::array_t<uint64_t> vals;
        ierr = getValuesArray(&vals);
        if (ierr >= 0) {
            ierr = ibis::column::searchSortedICD(vals, rng, hits);
        }
        break;}
    default: {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- FQ_Variable["
            << (thePart ? thePart->name() : "?") << '.'
            << m_name << "]::searchSorted(" << rng.colName() << " IN ...) "
            << "does not yet support column type "
            << ibis::TYPESTRING[(int)m_type];
        ierr = -5;
        break;}
    } // switch (m_type)
    return (ierr < 0 ? ierr : 0);
} // FQ_Variable::searchSorted

ibis::array_t<double>*
FQ_Variable::selectDoubles(const ibis::bitvector& mask) const {
    ibis::array_t<double>* array = new ibis::array_t<double>;
    return selectData(mask, array);
} // FQ_Variable::selectDoubles

ibis::array_t<float>*
FQ_Variable::selectFloats(const ibis::bitvector& mask) const {
    ibis::array_t<float>* array = new ibis::array_t<float>;
    return selectData(mask, array);
} // FQ_Variable::selectFloats

ibis::array_t<signed char>*
FQ_Variable::selectBytes(const ibis::bitvector& mask) const {
    ibis::array_t<signed char>* array = new ibis::array_t<signed char>;
    return selectData(mask, array);
} // FQ_Variable::selectBytes

ibis::array_t<int32_t>*
FQ_Variable::selectInts(const ibis::bitvector& mask) const {
    ibis::array_t<int32_t>* array = new ibis::array_t<int32_t>;
    return selectData(mask, array);
} // FQ_Variable::selectInts
ibis::array_t<int64_t>*
FQ_Variable::selectLongs(const ibis::bitvector& mask) const {
    ibis::array_t<int64_t>* array = new ibis::array_t<int64_t>;
    return selectData(mask, array);
} // FQ_Variable::selectLongs
ibis::array_t<uint32_t>*
FQ_Variable::selectUInts(const ibis::bitvector& mask) const {
    ibis::array_t<uint32_t>* array = new ibis::array_t<uint32_t>;
    return selectData(mask, array);
} // FQ_Variable::selectUInts

ibis::array_t<uint64_t>*
FQ_Variable::selectULongs(const ibis::bitvector& mask) const {
    ibis::array_t<uint64_t>* array = new ibis::array_t<uint64_t>;
    return selectData(mask, array);
} // FQ_Variable::selectULongs

template <typename T> ibis::array_t<T>*
FQ_Variable::selectData(const ibis::bitvector& mask,
			ibis::array_t<T>* array) const {
    if (! isValid("FQ_Variable::selectData")) return 0;

    ibis::array_t<T> prop;
    uint32_t i = 0;
    uint32_t tot = mask.cnt();
    ibis::horometer timer;
    if (ibis::gVerbose > 3) {
        LOGGER(ibis::gVerbose > 4)
            << "FQ_Variable[" << (thePart->name() ? thePart->name() : "?")
            << "." << name() << "]::selectData starting timer..";
        timer.start();
    }
#ifdef DEBUG
    LOGGER(1) << "Debug -- reading " << name() << ", mask.cnt() = " << tot
              << ", mask.size() = " << mask.size()
              << ", mask.bytes() = " << mask.bytes()
              << ", mask.size()*8/pagesize = "
              << mask.size()*8/ibis::fileManager::pageSize()
              << ", read all = "
              << (mask.bytes()/240 > mask.size()/ibis::fileManager::pageSize() ?
                  "yes" : "no");
#endif
    if (mask.size() == mask.cnt()) {
        getValuesArray(array);
        i = array->size();
        LOGGER(ibis::gVerbose > 1)
            << "FQ_Variable[" << (thePart->name() ? thePart->name() : "?")
            << "." << name() << "]::selectData using getValuesArray to retrieve "
            << i;
    }
    else if (mask.size() < 1048576 || tot+tot > mask.size() ||
             mask.bytes()/240 > mask.size()/ibis::fileManager::pageSize()) {
        // read all values than extract the ones marked with 1 in mask
        getValuesArray(&prop); // retrieving all values of this variable
        array->resize(tot);
        if (tot > prop.size()) tot = prop.size();
        const uint32_t nprop = prop.size();
        ibis::bitvector::indexSet index = mask.firstIndexSet();
        if (nprop >= mask.size()) {
            while (index.nIndices() > 0) {
                const ibis::bitvector::word_t *idx0 = index.indices();
                if (index.isRange()) {
                    for (uint32_t j = *idx0; j<idx0[1]; ++j, ++i) {
                        (*array)[i] = (prop[j]);
                    }
                }
                else {
                    for (uint32_t j = 0; j<index.nIndices(); ++j, ++i) {
                        (*array)[i] = (prop[idx0[j]]);
                    }
                }
                ++ index;
            }
        }
        else {
            while (index.nIndices() > 0) {
                const ibis::bitvector::word_t *idx0 = index.indices();
                if (*idx0 >= nprop) break;
                if (index.isRange()) {
                    for (uint32_t j = *idx0;
                         j<(idx0[1]<=nprop ? idx0[1] : nprop);
                         ++j, ++i) {
                        (*array)[i] = (prop[j]);
                    }
                }
                else {
                    for (uint32_t j = 0; j<index.nIndices(); ++j, ++i) {
                        if (idx0[j] < nprop)
                            (*array)[i] = (prop[idx0[j]]);
                        else
                            break;
                    }
                }
                ++ index;
            }
        }
        LOGGER(ibis::gVerbose > 1)
            << "FQ_Variable[" << (thePart->name() ? thePart->name() : "?")
            << "." << name()
            << "]::select using getValuesArray and extracted "
            << i;
    }
    else {
        // generate the coordinates and extract their values
        std::vector<uint64_t> coord;
        coord.reserve(tot);
        for (ibis::bitvector::indexSet ix = mask.firstIndexSet();
             ix.nIndices() > 0; ++ ix) {
            const ibis::bitvector::word_t *ind = ix.indices();
            if (ix.isRange()) {
                for (unsigned int j = ind[0]; j < ind[1]; ++ j)
                    coord.push_back(static_cast<int32_t>(j));
            }
            else {
                for (unsigned int j = 0; j < ix.nIndices(); ++j)
                    coord.push_back(static_cast<int32_t>(ind[j]));
            }
        }
	uint64_t nElements = coord.size()/varInfo.getNDims();
        array->resize(nElements);
        i = getPointValues(*array, coord);
        LOGGER(ibis::gVerbose > 1)
            << "FQ_Variable[" << (thePart->name() ? thePart->name() : "?")
            << "." << name() << "]::select using getPointValues. i = "
            << i;
    }
    if (i != tot) {
        array->resize(i);
        logWarning("select", "expects to retrieve %lu elements "
                   "but only got %lu", static_cast<long unsigned>(tot),
                   static_cast<long unsigned>(i));
    }
    else if (ibis::gVerbose > 3) {
        timer.stop();
        LOGGER(ibis::gVerbose >= 0)
            << "FQ_Variable[" << (thePart->name() ? thePart->name() : "?")
            << "." << name() << "]::select extracted " << tot << " value"
            << (tot > 1 ? "s" : "") << " out of " << mask.size() << " took "
            << timer.CPUTime() << " sec (CPU) and " << timer.realTime()
            << " sec (elapsed) time";
    }
    return array;
} // FQ_Variable::selectDoubles
