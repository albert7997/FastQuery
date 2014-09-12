#include "pnetCDFfile.h"
#include "fq.h"		// FastQuery::reportTiming()
#include <stack>

/*******************************************
 * Constructor & De-Constructor
 ********************************************/
PNETCDF::PNETCDF(const std::string fileName, const bool readOnly,
		 const std::string indexPath, const MPI_Comm comm) 
{
    _fileName = fileName;
    _indexPath = indexPath;

#ifndef FQ_NOMPI
    mpi_comm = comm;
    MPI_Group mpi_group;
    MPI_Comm_group(comm, &mpi_group);
    MPI_Group_size(mpi_group, &mpi_size);
    MPI_Group_rank(mpi_group, &mpi_rank);
#endif

    if (! __openFile(readOnly)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF:"
	    << " failed to open file \""
	    << _fileName.c_str() << "\"";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF:"
	    << " successfully opened file \""
	    << _fileName.c_str() << "\"";
    }
}

PNETCDF::~PNETCDF() 
{
    if (! __closeFile()) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::~PNETCDF:"
	    << " failed to close file \""
	    << _fileName.c_str() << "\"";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::~PNETCDF:"
	    << " successfully closed file \""
	    << _fileName.c_str() << "\"";
    }
}

/*******************************************
 * Public Functions
 ********************************************/
bool PNETCDF::getAttribute(const std::string &variable,
			   const std::string &attrName, void *values)
{
    return __getAttribute(variable, attrName, values);
} // PNETCDF::getAttribute

bool PNETCDF::setAttribute(const std::string &variable,
			   const std::string &attrName,
			   const void *values, const uint64_t len,
			   const FQ::DataType fqType)
{
    return __setAttribute(variable, attrName, values, len, fqType);
} // PNETCDF::setAttribute

bool PNETCDF::getAttributeInfo(const std::string &variable,
			       const std::string &attrName,
			       uint64_t *length, FQ::DataType *type)
{
    if (! __getAttributeInfo(variable, attrName, length, type)) {
	return false;
    }
    return true;
} // PNETCDF::getAttributeInfo

bool PNETCDF::createDataset(const std::string& variable,
			    const std::vector<uint64_t > dims,
			    const FQ::DataType fqType)
{
    if (dims.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::createDataset("
            << variable.c_str() << "):"
            << " the number of dataset dimensions is 0";
        return false;
    }
    if (! __createDataset(variable, dims, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::createDataset("
            << variable.c_str() << "):"
            << " failed to create the dataset";
        return false;
    }
    return true;
} // PNETCDF::createDataset

bool PNETCDF::setData(const std::string &variable, const void *data)
{
    if (! __writeData(variable, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::setData("
	    << variable.c_str() << "):"
	    << " failed to set data";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::setData("
	    << variable.c_str() << "):"
	    << " successfully set data";
	return true;
    }
} // PNETCDF::setData

bool PNETCDF::getData(const std::string &variable, void *data)
{
    if (! __readData(variable, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getData("
	    << variable.c_str() << "):"
	    << " failed to get data";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getData("
	    << variable.c_str() << "):"
	    << " successfully got data";
	return true;
    }
} // PNETCDF::getData

bool PNETCDF::setArrayData(const std::string &variable,
			   const std::vector<uint64_t> &offsets,
			   const std::vector<uint64_t> &counts,
			   const std::vector<uint64_t> &strides,
			   const void *data)
{
    if (! __writeArrayData(variable, offsets, counts, strides, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::setArrayData("
            << variable.c_str() << "):"
            << " failed to set data";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::setArrayData("
            << variable.c_str() << "):"
            << " successfully set data";
        return true;
    }
} // PNETCDF::setDataSection

bool PNETCDF::getArrayData(const std::string &variable, 
			   const std::vector<uint64_t> &offsets,
			   const std::vector<uint64_t> &counts, 
			   const std::vector<uint64_t> &strides, void *data)
{
    if (! __readArrayData(variable, offsets, counts, strides, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getArrayData("
	    << variable.c_str() << "):"
	    << " failed to get data";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getArrayData("
	    << variable.c_str() << "):"
	    << " successfully got data";
	return true;
    }
} // PNETCDF::getData

bool PNETCDF::getPointData(const std::string &variable,
			   const std::vector<uint64_t> &coords, void *data)
{
    std::vector <uint64_t> dims;
    if (! __getDatasetDimension(variable, dims)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getPointData("
	    << variable.c_str() << "):"
	    << " failed to get dataset dimension";
	return false;
    }
    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getPointData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }
    
    unsigned int nSelects = coords.size()/dims.size();
    unsigned int nElements = dims[0];
    for (unsigned int i=1; i<dims.size(); i++) nElements *= dims[i];

    if (nSelects > nElements/2) {
	// read all then return selected points
    	switch(fqType) {
	case FQ::FQT_FLOAT: {
	    float* ptr = (float*)data;
	    std::vector<float> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_DOUBLE: {
	    double* ptr = (double*)data;
	    std::vector<double> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_BYTE: {
	    signed char* ptr = (signed char*)data;
	    std::vector<signed char> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_UBYTE: {
	    unsigned char* ptr = (unsigned char*)data;
	    std::vector<unsigned char> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_SHORT: {
	    int16_t* ptr = (int16_t*)data;
	    std::vector<int16_t> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_USHORT: {
	    uint16_t* ptr = (uint16_t*)data;
	    std::vector<uint16_t> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_INT: {
	    int32_t* ptr = (int32_t*)data;
	    std::vector<int32_t> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_UINT: {
	    uint32_t* ptr = (uint32_t*)data;
	    std::vector<uint32_t> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_LONG: {
	    int64_t* ptr = (int64_t*)data;
	    std::vector<int64_t> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	case FQ::FQT_ULONG: {
	    uint64_t* ptr = (uint64_t*)data;
	    std::vector<uint64_t> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- PNETCDF::getPointData("
			    << variable.c_str() << "):"
			    << " invalid coords out of range";
			return false;
		    }
		    ptr[i] = values[pos];
		}
	    }
	    break;}
	default: {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- PNETCDF::getPointData("
		<< variable.c_str() << "):"
		<< " unknown FQ data type " << fqType;
	    return false;
	}
	}
    } else {
	// read point by point
        for (unsigned int i=0; i<nSelects; i++) {
	    std::vector<uint64_t> offsets;
	    offsets.clear();
	    std::vector<uint64_t> counts;
	    counts.clear();
	    std::vector<uint64_t> strides;
	    strides.clear();
	    for (unsigned int dim=0; dim<dims.size(); dim++) {
	    	offsets.push_back(coords[i*dims.size()+dim]);
		counts.push_back(1);
		strides.push_back(1);
	    }
	    int ierr;
	    switch(fqType){
	    case FQ::FQT_FLOAT: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((float*)data)[i]));
		break;}
	    case FQ::FQT_DOUBLE: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((double*)data)[i]));
		break;}
	    case FQ::FQT_BYTE: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((signed char*)data)[i]));
		break;}
	    case FQ::FQT_UBYTE: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((unsigned char*)data)[i]));
		break;}
	    case FQ::FQT_SHORT: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((int16_t*)data)[i]));
		break;}
	    case FQ::FQT_USHORT: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((uint16_t*)data)[i]));
		break;}
	    case FQ::FQT_INT: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((int32_t*)data)[i]));
		break;}
	    case FQ::FQT_UINT: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((uint32_t*)data)[i]));
		break;}
	    case FQ::FQT_LONG: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((int64_t*)data)[i]));
		break;}
	    case FQ::FQT_ULONG: {
	        ierr = __readArrayData(variable, offsets, counts, strides,
				       &(((uint64_t*)data)[i]));
		break;}
	    default:
		LOGGER(ibis::gVerbose > 0)
		    << "Warning -- PNETCDF::getPointData("
		    << variable.c_str() << "):"
		    << " unknown FQ data type " << fqType;
		return false;
	    }
	    if (! ierr) {
	        LOGGER(ibis::gVerbose > 0)
		    << "Warning -- PNETCDF::getPointData("
		    << variable.c_str() << "):"
		    << " failed to read array data";
	        return false;
	    }
	}
    }

    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::getPointData("
	<< variable.c_str() << "):"
	<< " successfully got data";
    return true;
} // PNETCDF::getPointData

bool PNETCDF::getBitmapKeys(const std::string &variable,
			    void *keys, const uint64_t mpi_idx)
{
    int i;
    bool berr = true;
    uint64_t start, end;
    int t1=0, t2;

    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::getBitmapKeys("
	<< variable.c_str() << "):"
	<< "_mpi_iter=" << _mpi_iter << ", mpi_idx=" << mpi_idx;

    if (! __getOffsets(variable, BitmapKeyColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getBitmapKeys("
	    << variable.c_str() << "):"
	    << " failed to get bitmap keys length";
	return false;
    }

    for (i=0; i<_mpi_iter-1; i++) {
	std::string datasetName = _indexPath;
	datasetName += variable;
	std::stringstream ss;
	ss << i;
	//datasetName += "."+ss.str();
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
	if (! berr) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- PNETCDF::getBitmapKeys("
		<< variable.c_str() << "):"
		<< " failed to get bitmap keys";
	    return false;
	}
    }
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::getBitmapKeys("
	<< variable.c_str() << "):"
	<< " successfully got bitmap keys";
    return true;
} // PNETCDF::getBitmapKeys

bool PNETCDF::setBitmapKeys(const std::string &variable,
			    const void *keys, 
			    const uint64_t nkeys,
			    const uint64_t mpi_idx)
{
    bool berr = true;
    uint64_t start, end;

    if (! __getOffsets(variable, BitmapKeyColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::setBitmapKeys("
	    << variable.c_str() << "):"
	    << " failed to get bitmap keys length";
	return false;
    }

    std::string datasetName = _indexPath;
    datasetName += variable;
    std::stringstream ss;
    ss << _mpi_iter;
    datasetName += "."+ss.str();
    datasetName += ".MPIbitmapKeys";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.resize(1);
    counts.resize(1);
    strides.resize(1);
    //offsets[0] = start;
    offsets[0] = 0; // NOTE: offset is always 0
    counts[0] = end-start;
    strides[0] = 1;
    berr = __writeArrayData(datasetName, offsets, counts, strides, keys);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::setBitmapKeys("
	    << variable.c_str() << "):"
	    << " failed to set bitmap keys";
	return false;
    } 
    LOGGER(ibis::gVerbose > 2)
        << "PNETCDF::setBitmapKeys("
	<< variable.c_str() << "):"
	<< " successfully set bitmap keys";
    return true;
} // PNETCDF::setBitmapKeys

bool PNETCDF::getBitmapKeyLength(const std::string &variable,
				 uint64_t *nkeys,
				 const uint64_t mpi_idx)
{
    bool berr = true;
    berr = __getOffsetLength(variable, BitmapKeyColIdx, nkeys, mpi_idx);
    if (! berr) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::getBitmapKeyLength("
	    << variable.c_str() << "):"
	    << " failed to get bitmap keys length";
	return false;
    } 
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::getBitmapKeyLength("
	<< variable.c_str() << "):"
	<< " successfully got bitmap keys length " << (uint64_t)(*nkeys);
    return true;
} // PNETCDF::getBitmapKeyLength

bool PNETCDF::getBitmapOffsets(const std::string &variable,
			       void *bitmapOffsets,
			       const uint64_t mpi_idx)
{
    bool berr = true;
    uint64_t start, end;

    if (! __getOffsets(variable, BitmapOffsetColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getBitmapOffsets("
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
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getBitmapOffsets("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offsets";
	return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "PNETCDF::getBitmapOffsets("
	<< variable.c_str() << "):"
	<< " successfully got bitmap offsets";
    return true;
} // PNETCDF::getBitmapOffsets

bool PNETCDF::setBitmapOffsets(const std::string &variable,
			       const void *bitmapOffsets, 
			       const uint64_t noffsets,
			       const uint64_t mpi_idx)
{
    bool berr = true;
    uint64_t start, end;

    if (! __getOffsets(variable, BitmapOffsetColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::setBitmapOffsets("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offsets length";
	return false;
    }

    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::setBitmapOffsets("
	<< variable.c_str() << "):"
	<< " start(" << start << "), end(" << end << ")";

    std::string datasetName = _indexPath;
    datasetName += variable;
    std::stringstream ss;
    ss << _mpi_iter;
    datasetName += "."+ss.str();
    datasetName += ".MPIbitmapOffsets";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.resize(1);
    counts.resize(1);
    strides.resize(1);
    //offsets[0] = start;
    offsets[0] = 0; // NOTE: it should be always 0
    counts[0] = end-start;
    strides[0] = 1;
    berr = __writeArrayData(datasetName, offsets, counts, strides,
			    bitmapOffsets);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::setBitmapOffsets("
	    << variable.c_str() << "):"
	    << " failed to set bitmap offsets";
	return false;
    } 
    LOGGER(ibis::gVerbose > 2)
        << "PNETCDF::setBitmapOffsets("
	<< variable.c_str() << "):"
	<< " successfully set bitmap offsets with size " << noffsets;
    return true;
} // PNETCDF::setBitmapOffsets

bool PNETCDF::getBitmapOffsetLength(const std::string &variable,
				    uint64_t *noffsets,
				    const uint64_t mpi_idx)
{
    bool berr = true;
    berr = __getOffsetLength(variable, BitmapOffsetColIdx, noffsets, mpi_idx);

    if (! berr) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::getBitmapOffsetsLength("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offsets length";
	return false;
    }
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::getBitmapOffsetsLength("
	<< variable.c_str() << "):"
	<< " successfully got bitmap offsets length "
	<< (uint64_t)(*noffsets);
    return true;
} // PNETCDF::getBitmapOffsetLength

FQ::DataType PNETCDF::getBitmapOffsetType(const std::string& variable)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmapOffsets";
  
    FQ::DataType fqType;   
    if (!  __getDatasetType(datasetName, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getBitmapOffsetsType("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offset type";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getBitmapOffsetsType("
	    << variable.c_str() << "):"
	    << " successfully got bitmap offset FQ type " << fqType;
    }
    return fqType;
} // PNETCDF::getBitmapOffsetsType

bool PNETCDF::getBitmapLength(const std::string &variable, uint64_t *nElements,
			      const uint64_t mpi_idx)
{
    bool berr = true;

    berr = __getOffsetLength(variable, BitmapColIdx, nElements, mpi_idx);
    if (! berr) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::getBitmapLength("
	    << variable.c_str() << "):"
	    << " failed to get bitmap length";
	return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "PNETCDF::getBitmapLength("
	<< variable.c_str() << "):"
	<< " successfully got bitmap length " << *nElements;
    return true;
} // PNETCDF::getBitmapLength


bool PNETCDF::readBitmap(const std::string& variable,
			 const uint64_t startoffset, 
			 const uint64_t endoffset,
			 uint32_t *data, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIbitmap";
    uint64_t count = endoffset - startoffset;

    if (count == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::readBitmap(" 
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

    uint64_t start, end;
    if (! __getOffsets(variable, BitmapColIdx, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::readBitmap("
	    << variable.c_str() << "):"
	    << " failed to get bitmap length";
	return false;
    }
    offsets[0] += start;

    if (! __readArrayData(datasetName, offsets, counts, strides, (void*)data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::readBitmap("
	    << variable.c_str() << "):"
	    << " failed to read bitmap with size " << count
	    << "[" << startoffset << " - " << endoffset << "]";
	return false;
    } 
    LOGGER(ibis::gVerbose > 2)
        << "PNETCDF::readBitmap("
	<< variable.c_str() << "):"
	<< " successfully read bitmap with size " << count
	<< "[" << startoffset << " - " << endoffset << "]";
    return true;
} // PNETCDF::readBitmap

bool PNETCDF::writeBitmap(const std::string& variable,
			  const uint64_t startoffset, 
			  const uint64_t endoffset,
			  const uint32_t *data,
			  const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    std::stringstream ss;
    ss << _mpi_iter;
    datasetName += "."+ss.str();
    datasetName += ".MPIbitmap";

    uint64_t count = endoffset - startoffset;
    if (count == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::writeBitmap "
	    << "(" << variable.c_str() << "):"
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

    if (! __writeArrayData(datasetName, offsets, counts, strides, (void*)data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::writeBitmap("
	    << variable.c_str() << "):"
	    << " failed to write bitmap with size " << (endoffset-startoffset)
	    << "[" << startoffset << " - " << endoffset << "]";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::writeBitmap("
	    << variable.c_str() << "):"
	    << " successfully wrote bitmap with size " << (endoffset-startoffset)
	    << "[" << startoffset << " - " << endoffset << "]";
    	return true;
    }
} // PNETCDF::writeBitmap

bool PNETCDF::createBitmapKeys(const std::string &variable,
			       const uint64_t nkeys, 
			       const FQ::DataType fqType,
			       const uint64_t mpi_iter,
			       const uint64_t mpi_idx)
{
    // get bitmapKey length
    uint64_t curLen = nkeys;
    _mpi_iter = mpi_iter;

    if (mpi_rank == 0 && mpi_iter > 0) {
	uint64_t start, end;
    	__getOffsets(variable, BitmapKeyColIdx, &start, &end,
		     mpi_iter*mpi_size-1);
	curLen = end;
    }
    MPI_Bcast(&curLen, 1, MPI_UNSIGNED_LONG, 0, mpi_comm);
    if (mpi_iter > 0 && curLen == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::createBitmapKeys("
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
    std::stringstream ss;
    ss << mpi_iter;
    datasetName += "."+ss.str();
    datasetName += ".MPIbitmapKeys";
    std::vector<uint64_t> dims;
    dims.push_back(curLen);
    uint64_t totLen;
    // create a dataset with the chunk size = nkeys from the first
    // processor
    totLen = extraLen;
    if (! __createDataset(datasetName, dims, fqType)) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::createBitmapKeys("
	    << variable.c_str() << "):"
	    << " failed to create the bitmap keys dataset";
	return false;
    }
    if (FastQuery::reportTiming()) {
	LOGGER(true) << "Statistic\tBitmapKeySize\t" << extraLen;
    } 

    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::createBitmapKeys("
    	<< variable.c_str() << "):"
    	<< " successfully create bitmap keys dataset";
    return true;
} // PNETCDF::createBitmapKeys

bool PNETCDF::createBitmapOffsets(const std::string &variable,
				  const uint64_t noffsets,
				  const FQ::DataType fqType,
				  const uint64_t mpi_iter,
				  const uint64_t mpi_idx)
{
    // get bitmapOffset length
    uint64_t curLen = noffsets;
    _mpi_iter = mpi_iter;
    if (mpi_rank == 0 && mpi_iter > 0) {
	uint64_t start, end;
    	__getOffsets(variable, BitmapOffsetColIdx, &start, &end,
		     mpi_iter*mpi_size-1);
	curLen = end;
    }
    MPI_Bcast(&curLen, 1, MPI_UNSIGNED_LONG, 0, mpi_comm);
    if (mpi_iter > 0 && curLen == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::createBitmapOffsets("
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
    std::stringstream ss;
    ss << mpi_iter;
    datasetName += "."+ss.str();
    datasetName += ".MPIbitmapOffsets";
    std::vector<uint64_t> dims;
    dims.push_back(curLen);
    uint64_t totLen;

    totLen = extraLen;
    if (! __createDataset(datasetName, dims, fqType)) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::createBitmapOffsets("
	    << variable.c_str() << "):"
	    << " failed to create the bitmap offsets dataset";
	return false;
    }
    if (FastQuery::reportTiming()) {
	LOGGER(true) << "Statistic\tBitmapOffsetSize\t" << extraLen;
    } 
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::createBitmapOffsets("
	<< variable.c_str() << "):"
	<< " successfully created the bitmap offsets dataset("
        << "mpi_iter=" << mpi_iter << ")";
    return true;
} // PNETCDF::createBitmapOffsets

bool PNETCDF::createBitmap(const std::string& variable, uint64_t nElements,
			   const uint64_t mpi_iter, uint64_t mpi_idx)
{
    // get bitmap length
    uint64_t curLen = nElements;
    _mpi_iter = mpi_iter;
    if (mpi_rank == 0 && mpi_iter > 0) {
	uint64_t start, end;
    	__getOffsets(variable, BitmapColIdx, &start, &end,
		     mpi_iter*mpi_size-1);
	curLen = end;
    }
    MPI_Bcast(&curLen, 1, MPI_UNSIGNED_LONG, 0, mpi_comm);
    if (mpi_iter > 0 && curLen == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::createBitmap("
	    << variable.c_str() << "):"
	    << " failed to get bitmap length";
	return false;
    }
    uint64_t extraLen = 0;
    uint64_t val = nElements;
    MPI_Allreduce(&val, &extraLen, 1, MPI_UNSIGNED_LONG, MPI_SUM, mpi_comm);

    // create the bitmap dataset collectively
    std::string datasetName = _indexPath;
    datasetName += variable;
    std::stringstream ss;
    ss << mpi_iter;
    datasetName += "."+ss.str();
    datasetName += ".MPIbitmap";
    std::vector<uint64_t> dims;
    dims.push_back(curLen);
    uint64_t totLen;

    totLen = extraLen;

    if (! __createDataset(datasetName, dims, FQ::FQT_UINT)) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::createBitmap("
	    << variable.c_str() << "):"
	    << " failed to create bitmap dataset";
	return false;
    }
    if (FastQuery::reportTiming()) {
	LOGGER(true) << "Statistic\tBitmapSize\t" << extraLen;
    } 

    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::createBitmap("
    	<< variable.c_str() << "):"
    	<< " successfully created bitmap dataset";
    return true;
} // PNETCDF::createBitmap

#ifndef FQ_NOMPI
bool PNETCDF::getActualRange(const std::string &variable, void *range)
{
    std::string attrName = "actualRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getAttribute(datasetName, attrName, range)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::getActualRange(" 
	    << variable.c_str() << "):"
 	    << " failed to get actual range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getActualRange(" 
	    << variable.c_str() << "):"
 	    << " successfully got actual range";
    	return true;
    }
} // PNETCDF::getActualRange

bool PNETCDF::setActualRange(const std::string &variable,
			     const void *range, const FQ::DataType fqType)
{
    std::string attrName = "actualRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __setAttribute(datasetName, attrName, range, 2, fqType)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::setActualRange(" 
	    << variable.c_str() << "):"
 	    << " failed to set actual range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::setActualRange(" 
	    << variable.c_str() << "):"
 	    << " successfully set actual range";
    	return true;
    }
} // PNETCDF::setActualRange

bool PNETCDF::getExpectedRange(const std::string &variable, void *range)
{
    std::string attrName = "expectedRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getAttribute(datasetName, attrName, range)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::getExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " failed to get expected range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " successfully got expected range";
    	return true;
    }
} // PNETCDF::getExpectedRange

bool PNETCDF::setExpectedRange(const std::string &variable,
			       const void *range, const FQ::DataType fqType)
{
    std::string attrName = "expectedRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __setAttribute(datasetName, attrName, range, 2, fqType)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::setExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " failed to set expected range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::setExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " successfully set expected range";
    	return true;
    }
} // PNETCDF::setExpectedRange
#endif

bool PNETCDF::getVariableInfo(const std::string &variable,
			      std::vector <uint64_t> &dims,
			      FQ::DataType *fqType)
{
    if (! __getDatasetType(variable, fqType)) {
    	LOGGER(ibis::gVerbose > 0)
    	    << "Warning -- " << "PNETCDF::getVariableInfo"
	    << "(" << variable.c_str() << "):"
    	    << " failed to get variable type";
	return false;
    } 
    if (! __getDatasetDimension(variable, dims)) {
    	LOGGER(ibis::gVerbose > 0)
    	    << "Warning -- " << "PNETCDF::getVariableInfo"
	    << "(" << variable.c_str() << "):"
    	    << " failed to get variable dimension";
	return false;
    } 
    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::getVariableInfo"
	<< "(" << variable.c_str() << "):"
    	<< " successfully get variable information"
	<< " dimension " << dims.size() << " FQ type " << *fqType;
    return true;
} // PNETCDF::getVariableInfo

bool PNETCDF::getAllVariables(const std::string &path,
			      std::vector<std::string> &variables)
{
    if (! __getAllVariables(path, variables)) {
    	LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getAllVariables(" 
	    << path.c_str() << "):"
	    << " failed to find variables"; 
	return false;
    } else {
    	LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::getAllVariables(" << path.c_str()
	    << "): successfully found number of variables: "
	    << variables.size(); 
	return true;
    }
}

std::string PNETCDF::getSortedFieldName()
{
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::getSortedFieldName"
    	<< " it is not implemented yet";
    return "";
} // PNETCDF::getSortedFieldName

#ifndef FQ_NOMPI
bool PNETCDF::setBitmapKeyLength(const std::string &variable,
				 const uint64_t nkeys,
				 const uint64_t mpi_iter)
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
            	    << "Warning -- PNETCDF::setBitmapKeyLength("
            	    << variable.c_str() << "):"
            	    << " failed to get bitmap keys length";
	    } else {
	    	curLen = end;
	    }
    	}
	if (berr) {
	    nkeysArray[0] += curLen;
	    for (unsigned int i=1; i<mpi_size; i++) {
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
	    << "Warning -- PNETCDF::setBitmapKeyLength("
	    << variable.c_str() << "):"
	    << " failed to set bitmap keys length";
    }
    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::setBitmapKeyLength("
    	<< variable.c_str() << "):"
    	<< " successfully set bitmap keys length";
    return true;
} // PNETCDF::setBitmapKeyLength

bool PNETCDF::setBitmapOffsetLength(const std::string &variable,
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
            	    << "Warning -- PNETCDF::setBitmapOffsetLength("
            	    << variable.c_str() << "):"
            	    << " failed to get bitmap offsets length";
	    } else {
	    	curLen = end;
	    }
    	}
	if (berr) {
	    noffsetsArray[0] += curLen;
	    for (unsigned int i=1; i<mpi_size; i++) {
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
	    << "Warning -- PNETCDF::setBitmapOffsetLength("
	    << variable.c_str() << "):"
	    << " failed to set bitmap offsets length";
	return false;
    }
    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::setBitmapOffsetLength("
    	<< variable.c_str() << "):"
    	<< " successfully set bitmap offsets length";
    return true;
} // PNETCDF::setBitmapOffsetLength

bool PNETCDF::setBitmapLength(const std::string &variable,
			      const uint64_t nElements,
			      const uint64_t mpi_iter)
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
            	    << "Warning -- PNETCDF::setBitmapLength("
            	    << variable.c_str() << "):"
            	    << " failed to set bitmap length";
	    } else {
	    	curLen = end;
	    }
    	}
	if (berr) {
	    nElementsArray[0] += curLen;
	    for (unsigned int i=1; i<mpi_size; i++) {
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
	    << "Warning -- PNETCDF::setBitmapLength("
	    << variable.c_str() << "):"
	    << " failed to set bitmap length";
    }
    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::setBitmapLength("
    	<< variable.c_str() << "):"
    	<< " successfully set bitmap length";
    return true;
}

bool PNETCDF::createOffsetTable(const std::string &variable,
				const uint64_t mpi_max_iter)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIoffsetTable";
    uint64_t nElements = mpi_max_iter * (uint64_t)mpi_size + 1;
    std::vector<uint64_t> dims;

    // the first element is always 0
    dims.push_back(nElements);
    // 3 columns to store the offset of bitmap, bitmapoffsets and bitmapkeys
    dims.push_back(3);
    
    // use FQT_INT instead of FQT_ULONG until pnetcdf supports it
    if (! __createDataset(datasetName, dims, FQ::FQT_ULONG)) {
        LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::createOffsetTable(" << variable.c_str()
	    << "): failed to create offset table dataset";
	return false;
    }
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::createOffsetTable(" << variable.c_str()
	<< "): successfully created offset table dataset with length "
	<< nElements;
    return true;
} // PNETCDF::createOffsetTable
#endif


/*******************************************
 * Private Functions
 ********************************************/
bool PNETCDF::__openFile(const bool readOnly) 
{
    //make sure the file is a valid PNETCDF file and that is exists...
    int ierr;

    if ( readOnly ) {
    	ierr = ncmpi_open(mpi_comm, _fileName.c_str(), NC_NOWRITE,
			  MPI_INFO_NULL, &_ncid);
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- " << "PNETCDF::openFile" 
		<< " cannot open non-exist file " << _fileName.c_str() 
		<< " with read permission"
		<< " (error code:" << ierr << "):";
	    _ncid = -1;	
	    return false;
	}
    } else {
    	ierr = ncmpi_open(mpi_comm, _fileName.c_str(), NC_WRITE,
			  MPI_INFO_NULL, &_ncid);
	/* if file doesn't exist, create one */
	if (ierr != NC_NOERR) {
	    ierr = ncmpi_create(mpi_comm, _fileName.c_str(),
				NC_CLOBBER|NC_64BIT_DATA, MPI_INFO_NULL,
				&_ncid);
	}
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
		<< "PNETCDF::openFile" 
		<< " cannot create file " << _fileName.c_str() 
		<< " with write permission"
		<< " (error code:" << ierr << "):";
	    _ncid = -1;
	    return false;
	}
    }
    ncmpi_enddef(_ncid); // leave define mode
    return true;
} // PNETCDF::__openFile

bool PNETCDF::__closeFile(){
    int ierr = ncmpi_close(_ncid);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::closeFile"
            << " fail to close file " << _fileName.c_str()
	    << " (error code:" << ierr << "):";
        return false;
    }
    _ncid = -1;
    return true;
} // PNETCDF::closeFile

nc_type PNETCDF::__getPNETCDFDataType(const FQ::DataType fqType)
{
    // no corresponding type for NC_CHAR and NC_SHORT? 
    switch (fqType) {
    case FQ::FQT_FLOAT:
	return NC_FLOAT;
    case FQ::FQT_DOUBLE:
	return NC_DOUBLE;
    case FQ::FQT_BYTE:
	return NC_BYTE;
    case FQ::FQT_UBYTE:
	return NC_UBYTE;
    case FQ::FQT_SHROT:
	return NC_SHORT;
    case FQ::FQT_USHORT: 
	return NC_USHORT;
    case FQ::FQT_INT:
	return NC_INT;
    case FQ::FQT_UINT: 
	return NC_UINT;
    case FQ::FQT_LONG: 
	return NC_INT64;
    case FQ::FQT_ULONG: 
	return NC_UINT64;
    default:
	return NC_NAT;
    }
} // PNETCDF::__getPNETCDFDataType

FQ::DataType PNETCDF::__getFQDataType(const nc_type type)
{
    switch(type) {
    case NC_FLOAT:
	return FQ::FQT_FLOAT;
    case NC_DOUBLE:
	return FQ::FQT_DOUBLE;
    case NC_CHAR:
    case NC_BYTE:
	return FQ::FQT_BYTE;
    case NC_UBYTE:
	return FQ::FQT_UBYTE;
    case NC_SHORT:
	return FQ::FQT_SHORT;
    case NC_USHORT:
	return FQ::FQT_USHORT;
    case NC_INT:
	return FQ::FQT_INT;
    case NC_UINT:
	return FQ::FQT_UINT;
    case NC_INT64:
	return FQ::FQT_LONG;
    case NC_UINT64:
	return FQ::FQT_ULONG;
    default:
	return FQ::FQT_UNKNOWN;
    }
} // PNETCDF::__getFQDataType

bool PNETCDF::__getDatasetType(const std::string &variable,
			       FQ::DataType *fqType)
{
    *fqType = FQ::FQT_UNKNOWN;

    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetType::ncmpi_inq_varid("
            << variable.c_str() << "):"
            << " failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    nc_type type;
    ierr = ncmpi_inq_vartype(_ncid, variableId, &type);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetType("
	    << variable.c_str() << "):"
	    << " failed to get the PNETCDF data type"
	    << " (error code:" << ierr << "):";
	return false;
    }
    
    switch (type) {
    case NC_CHAR:
    case NC_BYTE:
        *fqType = FQ::FQT_BYTE;
	break;
    case NC_UBYTE:
        *fqType = FQ::FQT_UBYTE;
	break;
    case NC_FLOAT:
        *fqType = FQ::FQT_FLOAT;
	break;
    case NC_DOUBLE:
        *fqType = FQ::FQT_DOUBLE;
	break;
    case NC_INT:
        *fqType = FQ::FQT_INT;
	break;
    case NC_UINT:
        *fqType = FQ::FQT_UINT;
	break;
    case NC_INT64:
        *fqType = FQ::FQT_LONG;
	break;
    case NC_UINT64:
        *fqType = FQ::FQT_ULONG;
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetType("
	    << variable.c_str() << "):"
            << " FQ does not yet support the PNETCDF data type";
	return false;
    }
    return true;
} // PNETCDF::__gtDatasetType

bool PNETCDF::__getDatasetLength(const std::string &variable, uint64_t *len)
{
    *len = 0;
    std::vector <uint64_t> dims;
    if (! __getDatasetDimension(variable, dims)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetLength"
	    << "(" << variable.c_str() << "):"
	    << " failed to get dataset dimension";
	return false;
    }
    *len = dims[0];
    for (unsigned int i=1; i<dims.size(); i++) {
	*len *= dims[i];
    }
    return true;
} // PNETCDF::__getDatasetSize

bool PNETCDF::__getDatasetDimension(const std::string &variable,
				    std::vector <uint64_t> &dims)
{
    dims.clear();
    int ndims, ierr, variableId;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetDimension::ncmpi_inq_varid("
            << variable.c_str() << "):"
            << " failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    ierr = ncmpi_inq_varndims(_ncid, variableId, &ndims);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetDimension("
            << variable.c_str() << "):"
            << " failed to get number of dimension"
	    << " (error code:" << ierr << "):";
        return false;
    }

    if (ndims == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetDimension("
            << variable.c_str() << "):"
            << " number of dimension is 0";
        return false;
    }

    int dimIds[ndims];
    ierr = ncmpi_inq_vardimid(_ncid, variableId, dimIds);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getDatasetDimension("
            << variable.c_str() << "):"
            << " failed to get dimension ids"
	    << " (error code:" << ierr << "):";
        return false;
    }

    for (int i=0; i<ndims; i++) {
        MPI_Offset len;
        ierr = ncmpi_inq_dimlen(_ncid, dimIds[i], &len);
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
	    	<< "Warning -- PNETCDF::__getDatasetDimension("
		<< variable.c_str() << "):"
		<< " failed to get dimension length"
	        << " (error code:" << ierr << "):";
	    return false;
	}
        dims.push_back(len);
    }
    return true;
} // PNETCDF::__getDatasetDimension

bool PNETCDF::__getAttribute(const std::string &variable,
			     const std::string &attrName, void *values)
{
    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getAttribute::ncmpi_inq_varid("
            << variable.c_str() << "): failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getAttribute("
	    << variable.c_str() << "): failed to get dataset type";
	return false;
    }

    switch(fqType) {
    case FQ::FQT_FLOAT:
	ierr = ncmpi_get_att_float
	    (_ncid, variableId, attrName.c_str(), (float *)values);
	break;
    case FQ::FQT_DOUBLE:
	ierr = ncmpi_get_att_double
	    (_ncid, variableId, attrName.c_str(), (double *)values);
	break;
    case FQ::FQT_BYTE:
	ierr = ncmpi_get_att_uchar
	    (_ncid, variableId, attrName.c_str(), (signed char *)values);
	break;
    case FQ::FQT_UBYTE:
	ierr = ncmpi_get_att_uchar
	    (_ncid, variableId, attrName.c_str(), (unsigned char *)values);
	break;
    case FQ::FQT_SHORT:
	ierr = ncmpi_get_att_short
	    (_ncid, variableId, attrName.c_str(), (int16_t *)values);
	break;
    case FQ::FQT_USHORT:
	ierr = ncmpi_get_att_ushort
	    (_ncid, variableId, attrName.c_str(), (uint16_t*)values);
	break;
    case FQ::FQT_INT:
	ierr = ncmpi_get_att_int
	    (_ncid, variableId, attrName.c_str(), (int32_t *)values);
	break;
    case FQ::FQT_UINT:
	ierr = ncmpi_get_att_uint
	    (_ncid, variableId, attrName.c_str(), (uint32_t *)values);
	break;
    case FQ::FQT_LONG:
	ierr = ncmpi_get_att_long
	    (_ncid, variableId, attrName.c_str(), (int64_t *)values);
	break;
    case FQ::FQT_ULONG:
	ierr = ncmpi_get_att_ulonglong
	    (_ncid, variableId, attrName.c_str(), (uint64_t*)values);
	break;
    default:
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__getAttribute"
	    << "Unsupported type";
	return false;
    }
    if (ierr != NC_NOERR ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getAttribute("
            << variable.c_str() << "," << attrName.c_str() << "):"
            << " fail to get attribute values"
	    << " (error code:" << ierr << "):";
        return false;
    }
    return true;

} // PNETCDF::__getAttribute

bool PNETCDF::__setAttribute(const std::string &variable,
			     const std::string &attrName, 
			     const void *values, const size_t len,
			     const FQ::DataType fqType)
{
    if (len == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__setAttribute"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " cannot create an attribute with size 0";
	return false;
    }

    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__setAttribute::ncmpi_inq_varid("
            << variable.c_str() << "):"
            << " failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    ncmpi_redef(_ncid); // enter define mode

    int attId;
    if (ncmpi_inq_attid(_ncid, variableId, attrName.c_str(), &attId)
	== NC_NOERR) {
        // delete the existing attribute
        int ierr = ncmpi_del_att(_ncid, variableId, attrName.c_str());
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- PNETCDF::__setAttribute("
		<< variable.c_str() << "," << attrName.c_str() << "):"
		<< " failed to remove existing attribute"
		<< " (error code:" << ierr << "):";
	    ncmpi_enddef(_ncid); // leave define mode
	    return false;
        }
    }

    // create attribute
    nc_type attrType;
    attrType = __getPNETCDFDataType(fqType);
    if (attrType <= 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__setAttribute"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " not yet supported FQ data type " << fqType;
	ncmpi_enddef(_ncid); // leave define mode
	return false;
    }

    switch(attrType) {
    case NC_FLOAT:
	ierr = ncmpi_put_att_float(_ncid, variableId, attrName.c_str(),
				   attrType, len, (float *)values);
	break;
    case NC_DOUBLE:
	ierr = ncmpi_put_att_double(_ncid, variableId, attrName.c_str(),
				    attrType, len, (double *)values);
	break;
    case NC_CHAR:
    case NC_BYTE:
	ierr = ncmpi_put_att_uchar(_ncid, variableId, attrName.c_str(),
				   attrType, len, (signed char *)values);
	break;
    case NC_UBYTE:
	ierr = ncmpi_put_att_uchar(_ncid, variableId, attrName.c_str(),
				   attrType, len, (unsigned char *)values);
	break;
    case NC_SHORT:
	ierr = ncmpi_put_att_short(_ncid, variableId, attrName.c_str(),
				   attrType, len, (int16_t*)values);
	break;
    case NC_USHORT:
	ierr = ncmpi_put_att_ushort(_ncid, variableId, attrName.c_str(),
				    attrType, len, (uint16_t*)values);
	break;
    case NC_INT:
	ierr = ncmpi_put_att_int(_ncid, variableId, attrName.c_str(),
				 attrType, len, (int32_t *)values);
	break;
    case NC_UINT:
	ierr = ncmpi_put_att_int(_ncid, variableId, attrName.c_str(),
				 attrType, len, (uint32_t *)values);
	break;
    case NC_LONG:
	ierr = ncmpi_put_att_long(_ncid, variableId, attrName.c_str(),
				  attrType, len, (int64_t*)values);
	break;
    case NC_ULONG:
	ierr = ncmpi_put_att_ulong(_ncid, variableId, attrName.c_str(),
				   attrType, len, (uint64_t*)values);
	break;
    default:
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__setAttribute"
	    << "Unsupported type";
	ncmpi_enddef(_ncid); // leave define mode
	return false;
    }
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__setAttribute"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " failed to set attribute"
	    << " (error code:" << ierr << "):";
	ncmpi_enddef(_ncid); // leave define mode
	return false;
    }
    ncmpi_enddef(_ncid); // leave define mode
    return true;
} // PNETCDF::__setAttribute

bool PNETCDF::__getAttributeInfo(const std::string &variable,
				 const std::string &attrName,
				 uint64_t *length, FQ::DataType *fqType)
{
    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getAttributeInfo::ncmpi_inq_varid("
            << variable.c_str() << "):"
            << " failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    nc_type type;
    ierr = ncmpi_inq_atttype(_ncid, variableId, attrName.c_str(), &type);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__getAttributeInfo"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " failed to get attribute type"
 	    << " (error code:" << ierr << "):";
	return false;
    }
    *fqType = __getFQDataType(type);
    if (*fqType == FQ::FQT_UNKNOWN) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getAttributeInfo("
            << "(" << variable.c_str() << "," << attrName.c_str() << "):"
            << " failed to get the attribute type";
        return false;
    }
    MPI_Offset ltmp;
    ierr = ncmpi_inq_attlen(_ncid, variableId, attrName.c_str(), &ltmp);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__getAttributeInfo"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " failed to get attribute length"
 	    << " (error code:" << ierr << "):";
	return false;
    }
    *length = ltmp;
    return true;
} // PNETCDF::__getAttributeInfo

bool PNETCDF::__readData(const std::string &variable, void *data)
{
    if (data == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readData("
	    << variable.c_str() << "):"
	    << " no data needs to be read";
	return false;
    }

    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readData::ncmpi_inq_varid("
            << variable.c_str() << "):"
            << " failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    switch (fqType) {
    case FQ::FQT_BYTE:
        ierr = ncmpi_get_var_schar_all(_ncid, variableId, (signed char*)data);
	break;
    case FQ::FQT_FLOAT:
        ierr = ncmpi_get_var_float_all(_ncid, variableId, (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = ncmpi_get_var_double_all(_ncid, variableId, (double*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = ncmpi_get_var_short_all(_ncid, variableId, (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
        ierr = ncmpi_get_var_ushort_all(_ncid, variableId, (int64_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = ncmpi_get_var_int_all(_ncid, variableId, (int32_t*)data);
	break;
    case FQ::FQT_UINT:
        ierr = ncmpi_get_var_uint_all(_ncid, variableId, (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
        ierr = ncmpi_get_var_long_all(_ncid, variableId, (int64_t*)data);
	break;
    case FQ::FQT_ULONG:
        ierr = ncmpi_get_var_ulonglong_all(_ncid, variableId, (uint64_t*)data);
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the PNETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readData("
            << variable.c_str() << "):"
            << " failed to read data (error code:" << ierr << "):";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::__readData("
            << variable.c_str() << "):"
            << " successfully read data";
        return true;
    }
} // PNETCDF::__readData

bool PNETCDF::__writeData(const std::string &variable, const void *data)
{
    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeData::ncmpi_inq_varid("
            << variable.c_str() << "):"
            << " failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    switch (fqType) {
    case FQ::FQT_FLOAT:
        ierr = ncmpi_put_var_float_all(_ncid, variableId, (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = ncmpi_put_var_double_all(_ncid, variableId, (double*)data);
	break;
    case FQ::FQT_BYTE:
        ierr = ncmpi_put_var_schar_all(_ncid, variableId, (signed char*)data);
	break;
    case FQ::FQT_UBYTE:
        ierr = ncmpi_put_var_uchar_all(_ncid, variableId, (unsigned char*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = ncmpi_put_var_short_all(_ncid, variableId, (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
        ierr = ncmpi_put_var_ushort_all(_ncid, variableId, (uint64_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = ncmpi_put_var_int_all(_ncid, variableId, (int32_t*)data);
	break;
    case FQ::FQT_UINT:
        ierr = ncmpi_put_var_uint_all(_ncid, variableId, (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
        ierr = ncmpi_put_var_long_all(_ncid, variableId, (int64_t*)data);
	break;
    case FQ::FQT_ULONG:
        ierr = ncmpi_put_var_ulong_all(_ncid, variableId, (uint64_t*)data);
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the PNETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeData("
            << variable.c_str() << "):"
            << " failed to write data (error code:" << ierr << "):";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::__writeData("
            << variable.c_str() << "):"
            << " successfully wrote data";
        return true;
    }
} // PNETCDF::__writeData

bool PNETCDF::__createDataset(const std::string& variable, 
			      const std::vector<uint64_t > dims,
			      const FQ::DataType fqType)
{
    if (_ncid < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__createDataset("
            << variable.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false; 
    }

    if (dims.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__createDataset("
	    << variable.c_str() << "):"
	    << " number of dimension is 0";
	return false;
    }

    if (variable.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::__createDataset("
            << variable.c_str() << "):"
            << " variable is an empty string";
        return false;
    }

    if (variable[0] != '/') {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::__createDataset("
            << variable.c_str() << "):"
            << " variable must be given as an absolute path starting with '/'";
        return false;
    }

    //check if dataset already exist
    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
#if 0
    // pnetcdf doesn't allow '[' or ']' or ',' as string, 
    // so convert them into '(', ')', '_', respectively
    pos = varName.find('[');
    while( pos != varName.npos ) {
	varName[pos] = '('; 	
	pos = varName.find('[');
    }
    pos = varName.find(']');
    while( pos != varName.npos ) {
	varName[pos] = ')'; 	
	pos = varName.find(']');
    }
    pos = varName.find(',');
    while( pos != varName.npos ) {
	varName[pos] = '_'; 	
	pos = varName.find(',');
    }
    // end conversion
#endif
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr == NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- ncmpi_inq_varid("
            << variable.c_str() << "): dataset already exists and "
	    "PNETCDF does not allow it to be re-created";
        return false;
    }

    // get PNETCDF data type
    nc_type type = __getPNETCDFDataType(fqType);
    if (type <= 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__createDataset("
	    << variable.c_str() << "):"
	    << " not yet supported FQ data type " << fqType;
	return false;
    }

    ncmpi_redef(_ncid); // enter define mode

    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::__createDataset("
	<< variable.c_str() << "):" << "dims.size=" << dims.size(); 

    // create dimension
    int dimIds[dims.size()];
    for (unsigned int i=0; i<dims.size(); i++) {
	// create dimension name by replacing '/' with _ and add dimension
	// rank at the end
	std::string dimName = varName;
	int pos = dimName.find('/');
	while( pos != dimName.npos ) {
	    dimName[pos] = '_'; 	
	    pos = dimName.find('/');
	}

	std::ostringstream oss;
	oss << dimName << "__" << i;
	dimName = oss.str();
	int ierr = ncmpi_def_dim(_ncid, dimName.c_str(), dims[i], &(dimIds[i]));
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- PNETCDF::__createDataset("
		<< variable.c_str() << "):"
		<< " failed to create dimension"
		<< "(error code:" << ierr << "):" << dimName.c_str();
	    ncmpi_enddef(_ncid); // leave define mode
	    return false;
	}
    }

    ierr = ncmpi_def_var(_ncid, varName.c_str(), type, dims.size(),
			 dimIds, &variableId);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__createDataset("
	    << variable.c_str() << "):"
	    << " failed to create the dataset"
	    << " (error code:" << ierr << "):";
        ncmpi_enddef(_ncid); // leave define mode
	return false;
    }
    ncmpi_enddef(_ncid); // leave define mode
    return true;
} // PNETCDF::__createDataset

bool PNETCDF::__readArrayData(const std::string& variable, 
			      const std::vector<uint64_t> &offsets, 
			      const std::vector<uint64_t> &counts, 
			      const std::vector<uint64_t> &strides, void *data)
{
    if (offsets.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data offset is empty";
	return false;
    }

    if (counts.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data count is empty";
	return false;
    }

    if (counts.size() != offsets.size()) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " size of data count and offset does not match";
	return false;
    }

    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readArrayData::ncmpi_inq_varid("
            << variable.c_str() << "): failed to get id"
	    << " (error code:" << ierr << "):";
        return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readArrayData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    switch (fqType) {
    case FQ::FQT_FLOAT:
        ierr = ncmpi_get_vars_float_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (float*)data);
	break;
    case FQ::FQT_DOUBLE:
	ierr = ncmpi_get_vars_double_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (double*)data);
	break;
    case FQ::FQT_BYTE:
	ierr = ncmpi_get_vars_schar_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset *)(&strides[0]),
	     (signed char*)data);
	break;
    case FQ::FQT_UBYTE:
	ierr = ncmpi_get_vars_uchar_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset *)(&strides[0]),
	     (signed char*)data);
	break;
    case FQ::FQT_SHORT:
	ierr = ncmpi_get_vars_long_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
	ierr = ncmpi_get_vars_ulonglong_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (uint16_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = ncmpi_get_vars_int_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (int32_t*)data);
	break;
    case FQ::FQT_UINT:
	ierr = ncmpi_get_vars_uint_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
	ierr = ncmpi_get_vars_long_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (int64_t*)data);
	break;
    case FQ::FQT_ULONG:
	ierr = ncmpi_get_vars_ulonglong_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (uint64_t*)data);
	break;
    default:
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__readArrayData("
	    << variable.c_str() << "):"
	    << " FQ does not yet support the PNETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__readArrayData("
            << variable.c_str() << "):"
            << " failed to read data (error code:" << ierr << "):";
        return false;
    } 
    return true;
} // PNETCDF::__readArrayData

bool PNETCDF::__writeArrayData(const std::string& variable, 
			       const std::vector<uint64_t> &offsets, 
			       const std::vector<uint64_t> &counts, 
			       const std::vector<uint64_t> &strides, 
			       const void *data)
{
    if (counts.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data count is empty";
	return true;
    }

    uint64_t totalElements = counts[0];
    for (unsigned int i=1; i<counts.size(); i++) {
	totalElements*=counts[i];
    }

    if (totalElements == 0) {
	LOGGER(ibis::gVerbose > 2)
	    << "Warning -- PNETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " size of total data count is 0";
	return true;
    }

    if (counts.size() != offsets.size()) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " size of data count and offset does not match";
	return false;
    }
    
    LOGGER(ibis::gVerbose > 2)
	<< "PNETCDF::_writeArrayData:"
	<< "offsets[0]=" << offsets[0] << ", "
	<< "counts[0]=" << counts[0] << ", "
	<< "strides[0]=" << strides[0];

    int variableId, ierr;
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    ierr = ncmpi_inq_varid(_ncid, varName.c_str(), &variableId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeArrayData::ncmpi_inq_varid("
            << variable.c_str()
            << "): failed to get id (error code:" << ierr << "):";
        return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeArrayData(" << variable.c_str()
	    << "): failed to get dataset type";
	return false;
    }

    switch (fqType) {
    case FQ::FQT_FLOAT:
	ierr = ncmpi_put_vars_float_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = ncmpi_put_vars_double_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (double*)data);
	break;
    case FQ::FQT_BYTE:
        ierr = ncmpi_put_vars_schar_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (signed char*)data);
	break;
    case FQ::FQT_UBYTE:
        ierr = ncmpi_put_vars_schar_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (unsigned char*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = ncmpi_put_vars_long_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
	ierr = ncmpi_put_vars_ulonglong_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (uint16_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = ncmpi_put_vars_int_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (int32_t*)data);
	break;
    case FQ::FQT_UINT:
	ierr = ncmpi_put_vars_uint_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
        ierr = ncmpi_put_vars_long_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (int64_t*)data);
	break;
    case FQ::FQT_ULONG:
	ierr = ncmpi_put_vars_ulonglong_all
	    (_ncid, variableId, (MPI_Offset*)(&offsets[0]),
	     (MPI_Offset*)(&counts[0]), (MPI_Offset*)(&strides[0]),
	     (uint64_t*)data);
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__writeArrayData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the PNETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- PNETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " failed to write data"
 	    << "(error code:" << ierr << "):";
	return false;
    }
    return true;
} // PNETCDF::__writeArrayData

bool PNETCDF::__getAllVariables(const std::string &path,
				std::vector<std::string> &variables)
{
    variables.clear();

    if (_ncid < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getAllVariables("
            << path.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false; 
    }
    // if the path is not valid file location, return false
    if (path.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::__getAllVariables(" 
	    << path.c_str() << "):"
 	    << " path is an empty string";
	return false;
    }

    if (path[0] != '/') {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "PNETCDF::__getAllVariables(" 
	    << path.c_str() << "):"
 	    << " not an absolute path";
	return false;
    }

    __traverseVariables("/", _ncid, variables);
    if (variables.size() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "PNETCDF::__getAllVariables "
            << " no variable is found";
    }
    return true;
} // PNETCDF::__getAllVariables

void PNETCDF::__traverseVariables(const std::string &path,
				  const int parentId,
				  std::vector<std::string> &variables)
{
    int nvars;
    ncmpi_inq_nvars(parentId, &nvars);
    if (nvars > 0) {
	for (int i=0; i<nvars; i++) {
	    char varname[NC_MAX_NAME+1];
	    ncmpi_inq_varname(parentId, i, varname);
	    std::string str = varname;
	    if (str.find(".bitmap") == str.npos &&
		str.find(".bitmapKeys") == str.npos &&
		str.find(".bitmapOffsets") == str.npos ) {
		std::string variable = path;
		if (path.compare("/") != 0) {
		    variable += '/';
		}
		variable += str;
		variables.push_back(variable);
	    }
	}
    }
#if 0
    int ngroups;
    ncmpi_inq_grps(parentId, &ngroups, NULL);
    if (ngroups > 0) {
        int groupIds[ngroups];
        ncmpi_inq_grps(parentId, NULL, groupIds);
        for (int i=0; i<ngroups; i++) {
            size_t len;
            ncmpi_inq_grpname_full(groupIds[i], &len, NULL);
            char groupName[len+1];
            ncmpi_inq_grpname_full(groupIds[i], NULL, groupName);
	    std::string groupNameStr = groupName;
            __traverseVariables(groupNameStr, groupIds[i], variables);
        }
    }
#endif
} // PNETCDF::traverseVariables

#ifndef FQ_NOMPI
bool PNETCDF::__getOffsets(const std::string &variable, const int column, 
			   uint64_t *start, uint64_t *end,
			   const uint64_t mpi_idx)
{
    bool berr = true;
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".MPIoffsetTable";
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    offsets.push_back(mpi_idx);
    counts.push_back(2);  // read the start and end offset
    strides.push_back(1);
    offsets.push_back(column); 
    counts.push_back(1);
    strides.push_back(1);
    uint64_t vals[2]; // TODO: indexes are stored in "int" instead of "ulong"
    //int vals[2];
    if (! __readArrayData(datasetName, offsets, counts, strides, vals)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getOffsets(" << variable.c_str()
	    << "): failed to get offset table";
	return false;
    }
    if (vals[1] < vals[0]) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getOffsets(" << variable.c_str()
	    << "): invalid offset: " << vals[1] << " should be greater than "
	    << vals[0];
	return false;
    }
    *start = vals[0];
    *end = vals[1];
    LOGGER(ibis::gVerbose > 2)
    	<< "PNETCDF::__getOffsets("
    	<< variable.c_str() << "):"
    	<< " successfully got offset at position "
	<< offsets[0] << " with values "
	<< "("<<*start<<"-"<<*end<<"):";
    return true;
} // PNETCDF::__getOffsets

bool PNETCDF::__getOffsetLength(const std::string &variable, const int column, 
				uint64_t* len, const uint64_t mpi_idx)
{
    uint64_t start;
    uint64_t end;
    if (! __getOffsets(variable, column, &start, &end, mpi_idx)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- PNETCDF::__getOffsetLength("
	    << variable.c_str() << "):"
	    << " failed to get offset table";
	return false;
    }
    *len = (end-start);
    return true;
} // PNETCDF::__getOffsetLength
#endif
