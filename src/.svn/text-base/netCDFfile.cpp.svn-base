#include "netCDFfile.h"
#include <stack>

/*******************************************
 * Constructor & De-Constructor
 ********************************************/
NETCDF::NETCDF(const std::string fileName, const bool readOnly,
	       const std::string indexPath) 
{
    _fileName = fileName;
    _indexPath = indexPath;
    _fileId = -1;

    if (! __openFile(readOnly)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF:"
	    << " failed to open file \""
	    << _fileName.c_str() << "\"";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF:"
	    << " successfully opened file \""
	    << _fileName.c_str() << "\"";
    }
}

NETCDF::~NETCDF() 
{
    if (! __closeFile()) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::~NETCDF:"
	    << " failed to close file \""
	    << _fileName.c_str() << "\"";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::~NETCDF:"
	    << " successfully closed file \""
	    << _fileName.c_str() << "\"";
    }
}

/*******************************************
 * Public Functions
 ********************************************/
bool NETCDF::getAttribute(const std::string &variable,
			  const std::string &attrName, void *values)
{
    return __getAttribute(variable, attrName, values);
} // NETCDF::getAttribute

bool NETCDF::setAttribute(const std::string &variable,
			  const std::string &attrName,
			  const void *values,
			  const uint64_t len,
			  const FQ::DataType fqType)
{
    return __setAttribute(variable, attrName, values, len, fqType);
} // NETCDF::setAttribute

bool NETCDF::getAttributeInfo(const std::string &variable,
			      const std::string &attrName,
                              uint64_t *length, FQ::DataType *type)
{
    if (! __getAttributeInfo(variable, attrName, length, type)) {
	return false;
    }
    return true;
} // NETCDF::getAttributeInfo

bool NETCDF::createDataset(const std::string& datasetName,
			   const std::vector<uint64_t> dims,
			   const FQ::DataType fqType)
{
    if (dims.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::createDataset("
            << datasetName.c_str() << "):"
            << " the number of dataset dimensions is 0";
        return false;
    }
    if (! __createDataset(datasetName, dims, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::createDataset("
            << datasetName.c_str() << "):"
            << " failed to create the dataset";
        return false;
    }
    return true;
} // NETCDF::createDataset

bool NETCDF::setData(const std::string &variable, const void *data)
{
    if (! __writeData(variable, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::setData("
	    << variable.c_str() << "):"
	    << " failed to set data";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::setData("
	    << variable.c_str() << "):"
	    << " successfully set data";
	return true;
    }
} // NETCDF::setData

bool NETCDF::getData(const std::string &variable, void *data)
{
    if (! __readData(variable, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getData("
	    << variable.c_str() << "):"
	    << " failed to get data";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getData("
	    << variable.c_str() << "):"
	    << " successfully got data";
	return true;
    }
} // NETCDF::getData

bool NETCDF::setArrayData(const std::string &variable,
			  const std::vector<uint64_t> &offsets,
			  const std::vector<uint64_t> &counts,
			  const std::vector<uint64_t> &strides, const void *data)
{
    if (! __writeArrayData(variable, offsets, counts, strides, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::setArrayData("
            << variable.c_str() << "):"
            << " failed to set data";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 3)
            << "NETCDF::setArrayData("
            << variable.c_str() << "):"
            << " successfully set data";
        return true;
    }
} // NETCDF::setDataSection

bool NETCDF::getArrayData(const std::string &variable, 
			  const std::vector<uint64_t> &offsets,
			  const std::vector<uint64_t> &counts, 
			  const std::vector<uint64_t> &strides, void *data)
{
    if (! __readArrayData(variable, offsets, counts, strides, data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getArrayData("
	    << variable.c_str() << "):"
	    << " failed to get data";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getArrayData("
	    << variable.c_str() << "):"
	    << " successfully got data";
	return true;
    }
} // NETCDF::getData

bool NETCDF::getPointData(const std::string &variable,
			  const std::vector<uint64_t> &coords, void *data)
{
    std::vector <uint64_t> dims;
    if (! __getDatasetDimension(variable, dims)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getPointData("
	    << variable.c_str() << "):"
	    << " failed to get dataset dimension";
	return false;
    }
    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
	    std::vector<char> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- NETCDF::getPointData("
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
	    std::vector<char> values;
	    values.resize(nElements);
	    if (__readData(variable, &(values[0]))) {
		for (unsigned int i=0; i<nSelects; i++) {
		    uint64_t pos = coords[i*dims.size()];
		    for (unsigned int dim=1; dim<dims.size(); dim++) {
			pos = pos*dims[dim] + coords[i*dims.size()+dim];
		    }
		    if (pos > values.size()) {
			LOGGER(ibis::gVerbose > 0)
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
			    << "Warning -- NETCDF::getPointData("
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
		<< "Warning -- NETCDF::getPointData("
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
		    << "Warning -- NETCDF::getPointData("
		    << variable.c_str() << "):"
		    << " unknown FQ data type " << fqType;
		return false;
	    }
	    if (! ierr) {
	        LOGGER(ibis::gVerbose > 0)
		    << "Warning -- NETCDF::getPointData("
		    << variable.c_str() << "):"
		    << " failed to read array data";
	        return false;
	    }
	}
    }

    LOGGER(ibis::gVerbose > 2)
    	<< "NETCDF::getPointData("
	<< variable.c_str() << "):"
	<< " successfully got data";
    return true;
} // NETCDF::getPointData

bool NETCDF::getBitmapKeys(const std::string &variable, void *keys,
			   const uint64_t mpi_idx)
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapKeys("
            << variable.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false;
    }
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __readData(datasetName, keys) ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapKeys("
	    << variable.c_str() << "):"
	    << " failed to get bitmap keys";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getBitmapKeys("
	    << variable.c_str() << "):"
	    << " successfully got bitmap keys";
	return true;
    }
} // NETCDF::getBitmapKeys

bool NETCDF::setBitmapKeys(const std::string &variable, const void *keys, 
			   const uint64_t nkeys, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";

    if (! __writeData(datasetName, keys)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::setBitmapKeys("
	    << variable.c_str() << "):"
	    << " failed to set bitmap keys";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::setBitmapKeys("
	    << variable.c_str() << "):"
	    << " successfully set bitmap keys with size " << nkeys;
    	return true;
    }
} // NETCDF::setBitmapKeys

bool NETCDF::getBitmapKeyLength(const std::string &variable,
				uint64_t *nkeys, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getDatasetLength(datasetName, nkeys)) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- NETCDF::getBitmapKeyLength("
	    << variable.c_str() << "):"
	    << " failed to get bitmap key length";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getBitmapKeyLength("
	    << variable.c_str() << "):"
	    << " successfully got bitmap key length " << *nkeys;
    	return true;
    }
} // NETCDF::getBitmapKeyLength

bool NETCDF::getBitmapOffsets(const std::string &variable,
			      void *offsets, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";
    if (! __readData(datasetName, offsets)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapOffsets("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offsets";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getBitmapOffsets("
	    << variable.c_str() << "):"
	    << " successfully got bitmap offsets";
    	return true;
    }
} // NETCDF::getBitmapOffsets

bool NETCDF::setBitmapOffsets(const std::string &variable, const void *offsets, 
			      const uint64_t noffsets, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";

    if (! __writeData(datasetName, offsets)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::setBitmapOffsets("
	    << variable.c_str() << "):"
	    << " failed to set bitmap offsets";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::setBitmapOffsets("
	    << variable.c_str() << "):"
	    << " successfully set bitmap offsets with size " << noffsets;
    	return true;
    }
} // NETCDF::setBitmapOffsets

bool NETCDF::getBitmapOffsetLength(const std::string &variable,
				   uint64_t *noffsets, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";
    if (! __getDatasetLength(datasetName, noffsets)) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- NETCDF::getBitmapOffsetLength("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offset length";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getBitmapOffsetLength("
	    << variable.c_str() << "):"
	    << " successfully got bitmap offset length " << *noffsets;
    	return true;
    }
} // NETCDF::getBitmapOffsetLength

FQ::DataType NETCDF::getBitmapOffsetType(const std::string& variable)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";
  
    FQ::DataType fqType;   
    if (!  __getDatasetType(datasetName, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapOffsetsType("
	    << variable.c_str() << "):"
	    << " failed to get bitmap offset type";
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getBitmapOffsetsType("
	    << variable.c_str() << "):"
	    << " successfully got bitmap offset FQ type " << fqType;
    }
    return fqType;
} // NETCDF::getBitmapOffsetsType

bool NETCDF::getBitmapLength(const std::string &variable,
			     uint64_t *size,
			     const uint64_t mpi_idx)
{
    *size = 0;
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmap";

    if (!  __getDatasetLength(datasetName, size)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapLength("
	    << variable.c_str() << "):"
	    << " failed to get bitmap length";
	return false;
    } 
    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapLength("
	    << variable.c_str() << "):"
	    << " failed to get dataset data type";
	return false;
    }

    switch(fqType){
    case FQ::FQT_FLOAT:
        *size *= sizeof(float);
    	break;
    case FQ::FQT_DOUBLE:
	*size *= sizeof(double);
        break;
    case FQ::FQT_BYTE:
    case FQ::FQT_UBYTE:
        *size *= sizeof(char);
        break;
    case FQ::FQT_SHORT:
    case FQ::FQT_USHORT:
	*size *= sizeof(int16_t);
        break;
    case FQ::FQT_INT:
    case FQ::FQT_UINT:
	*size *= sizeof(int32_t);
        break;
    case FQ::FQT_LONG:
    case FQ::FQT_ULONG:
        *size *= sizeof(int64_t);
        break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::getBitmapLength("
	    << variable.c_str() << "):"
	    << " unknown FQ data type " << fqType;
	*size =0;
	return false;
    }
    LOGGER(ibis::gVerbose > 2)
    	<< "NETCDF::getBitmapLength("
	<< variable.c_str() << "):"
	<< " successfully got bitmap length" << *size;
    return true;
} // NETCDF::getBitmapLength


bool NETCDF::readBitmap(const std::string& variable,
			const uint64_t startoffset, 
			const uint64_t endoffset,
			uint32_t *data, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmap";
    
    uint64_t count = endoffset - startoffset;
    if (count == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::readBitmap "
	    << "(" << variable.c_str() << "):"
	    << " offset size is 0";
	return false;
    }

    std::vector<uint64_t> offsets;
    offsets.clear();
    offsets.push_back(startoffset);
    std::vector<uint64_t> counts;
    counts.clear();
    counts.push_back(count);
    std::vector<uint64_t> strides;
    strides.clear();
    strides.push_back(1);

    if (! __readArrayData(datasetName, offsets, counts, strides,
			  (void*)data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::readBitmap("
	    << variable.c_str() << "):"
	    << " failed to read bitmap with size " << count
	    << "[" << startoffset << " - " << endoffset << "]";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::readBitmap("
	    << variable.c_str() << "):"
	    << " successfully read bitmap with size " << count
	    << "[" << startoffset << " - " << endoffset << "]";
    	return true;
    }
} // NETCDF::readBitmap

bool NETCDF::writeBitmap(const std::string& variable,
			 const uint64_t startoffset, 
			 const uint64_t endoffset,
			 const uint32_t *data,
			 const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmap";

    uint64_t count = endoffset - startoffset;
    if (count == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::readBitmap "
	    << "(" << variable.c_str() << "):"
	    << " offset size is 0";
	return false;
    }

    std::vector<uint64_t> offsets;
    offsets.clear();
    offsets.push_back(startoffset);
    std::vector<uint64_t> counts;
    counts.clear();
    counts.push_back(count);
    std::vector<uint64_t> strides;
    strides.clear();
    strides.push_back(1);

    if (! __writeArrayData(datasetName, offsets, counts, strides,
			   (void*)data)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::writeBitmap("
	    << variable.c_str() << "):"
	    << " failed to write bitmap with size " << count
	    << "[" << startoffset << " - " << endoffset << "]";
	return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::writeBitmap("
	    << variable.c_str() << "):"
	    << " successfully wrote bitmap with size " << count
	    << "[" << startoffset << " - " << endoffset << "]";
    	return true;
    }
} // NETCDF::writeBitmap

bool NETCDF::createBitmapKeys(const std::string &variable,
			      const uint64_t nkeys, 
			      const FQ::DataType fqType,
			      const uint64_t mpi_iter,
			      const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";

    std::vector<uint64_t> dims;
    dims.resize(1);
    dims[0] = nkeys;
    if (! __createDataset(datasetName, dims, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::createBitmapKeys("
            << variable.c_str() << "):"
            << " failed to create bitmap keys dataset";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "NETCDF::createBitmapKeys("
        << variable.c_str() << "):"
        << " successfully create bitmap keys dataset";
    return true;
} // NETCDF::createBitmapKeys

bool NETCDF::createBitmapOffsets(const std::string &variable,
				 const uint64_t noffsets,
				 const FQ::DataType fqType,
				 const uint64_t mpi_iter, const uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapOffsets";

    std::vector<uint64_t> dims;
    dims.resize(1);
    dims[0] = noffsets;
    if (! __createDataset(datasetName, dims, fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::createBitmapOffsets("
            << variable.c_str() << "):"
            << " failed to create bitmap offsets dataset";
        return false;
    }
    LOGGER(ibis::gVerbose > 2)
        << "NETCDF::createBitmapOffsets("
        << variable.c_str() << "):"
        << " successfully created the bitmap offsets dataset";
    return true;
} // NETCDF::createBitmapOffsets

bool NETCDF::createBitmap(const std::string& variable, uint64_t nElements,
			  const uint64_t mpi_iter, uint64_t mpi_idx)
{
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmap";

    std::vector<uint64_t> dims;
    dims.clear();
    dims.push_back(nElements);
    if (! __createDataset(datasetName, dims, FQ::FQT_UINT)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::createBitmap(" 
	    << variable.c_str() << "):"
 	    << " failed to create bitmap dataset with size " << nElements;
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::createBitmap(" << variable.c_str()
	    << "): successfully created bitmap dataset with size " << nElements;
	return true;
    }
} // NETCDF::createBitmap

//#ifdef FQ_NOMPI
bool NETCDF::getActualRange(const std::string &variable, void *range)
{
    std::string attrName = "actualRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getAttribute(datasetName, attrName, range)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::getActualRange(" 
	    << variable.c_str() << "):"
 	    << " failed to get actual range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getActualRange(" 
	    << variable.c_str() << "):"
 	    << " successfully got actual range";
    	return true;
    }
} // NETCDF::getActualRange

bool NETCDF::setActualRange(const std::string &variable, const void *range,
			    const FQ::DataType fqType)
{
    std::string attrName = "actualRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __setAttribute(datasetName, attrName, range, 2, fqType)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::setActualRange(" 
	    << variable.c_str() << "):"
 	    << " failed to set actual range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::setActualRange(" 
	    << variable.c_str() << "):"
 	    << " successfully set actual range";
    	return true;
    }
} // NETCDF::setActualRange

bool NETCDF::getExpectedRange(const std::string &variable, void *range)
{
    std::string attrName = "expectedRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __getAttribute(datasetName, attrName, range)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::getExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " failed to get expected range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " successfully got expected range";
    	return true;
    }
} // NETCDF::getExpectedRange

bool NETCDF::setExpectedRange(const std::string &variable,
			      const void *range, const FQ::DataType fqType)
{
    std::string attrName = "expectedRange";
    std::string datasetName = _indexPath;
    datasetName += variable;
    datasetName += ".bitmapKeys";
    if (! __setAttribute(datasetName, attrName, range, 2, fqType)) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::setExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " failed to set expected range";
	return false;
    } else {
	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::setExpectedRange(" 
	    << variable.c_str() << "):"
 	    << " successfully set expected range";
    	return true;
    }
} // NETCDF::setExpectedRange
//#endif

bool NETCDF::getVariableInfo(const std::string &variable,
			     std::vector <uint64_t> &dims, FQ::DataType *fqType)
{
    if (! __getDatasetType(variable, fqType)) {
    	LOGGER(ibis::gVerbose > 0)
    	    << "Warning -- " << "NETCDF::getVariableInfo("
	    << variable.c_str() << "): failed to get variable type";
	return false;
    } 
    if (! __getDatasetDimension(variable, dims)) {
    	LOGGER(ibis::gVerbose > 0)
    	    << "Warning -- " << "NETCDF::getVariableInfo"
	    << "(" << variable.c_str() << "):"
    	    << " failed to get variable dimension";
	return false;
    } 
    LOGGER(ibis::gVerbose > 2)
    	<< "NETCDF::getVariableInfo(" << variable.c_str() << "):"
    	<< " successfully get variable information dimension "
	<< dims.size() << " FQ type " << *fqType;
    return true;
} // NETCDF::getVariableInfo

bool NETCDF::getAllVariables(const std::string &path,
			     std::vector<std::string> &variables)
{
    if (! __getAllVariables(path, variables)) {
    	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getAllVariables(" 
	    << path.c_str() << "):"
	    << " failed to find variables"; 
	return false;
    } else {
    	LOGGER(ibis::gVerbose > 2)
            << "NETCDF::getAllVariables(" 
	    << path.c_str() << "): successfully found number of variables: "
	    << variables.size(); 
	return true;
    }
}

std::string NETCDF::getSortedFieldName()
{
    LOGGER(ibis::gVerbose > 2)
	<< "NETCDF::getSortedFieldName"
    	<< " it is not implemented yet";
    return "";
} // NETCDF::getSortedFieldName

#ifndef FQ_NOMPI
bool NETCDF::setBitmapKeyLength(const std::string &variable,
				uint64_t nkeys, const uint64_t mpi_iter)
{
    LOGGER(ibis::gVerbose > 2)
	<< "NETCDF::setBitmapKeyLength"
    	<< " it is not implemented yet";
    return "";
} // NETCDF::setBitmapKeyLength

bool NETCDF::setBitmapOffsetLength(const std::string &variable
				   , uint64_t noffsets, const uint64_t mpi_iter)
{
    LOGGER(ibis::gVerbose > 2)
	<< "NETCDF::setBitmapOffsetLength"
    	<< " it is not implemented yet";
    return "";
} // NETCDF::setBitmapOffsetLength

bool NETCDF::setBitmapLength(const std::string &variable,
			     const uint64_t nElements, const uint64_t mpi_iter)
{
    LOGGER(ibis::gVerbose > 2)
	<< "NETCDF::setBitmapLength"
    	<< " it is not implemented yet";
    return "";
} // NETCDF::setBitmapLength

bool NETCDF::createOffsetTable(const std::string &variable,
			       const uint64_t mpi_max_iter)
{
    LOGGER(ibis::gVerbose > 2)
	<< "NETCDF::createOffsetTable"
    	<< " it is not implemented yet";
    return "";
} // NETCDF::createOffsetTable
#endif


/*******************************************
 * Private Functions
 ********************************************/
bool NETCDF::__openFile(const bool readOnly) 
{
    //make sure the file is a valid NETCDF file and that is exists...
    int ierr;
    if ( readOnly ) {
    	ierr = nc_open(_fileName.c_str(), NC_NOWRITE, &_fileId);
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- " << "NETCDF::openFile" 
		<< " cannot open non-exist file " << _fileName.c_str() 
		<< " with read permission"
		<< " (error code:" << ierr << "):";
	    _fileId = -1;	
	    return false;
	}
    } else {
    	ierr = nc_open(_fileName.c_str(), NC_WRITE, &_fileId);
	if (ierr != NC_NOERR) {
    	    ierr = nc_create(_fileName.c_str(), NC_NETCDF4, &_fileId);
	    if (ierr != NC_NOERR) {
		LOGGER(ibis::gVerbose > 2)
		    << "NETCDF::openFile" 
		    << " cannot create file " << _fileName.c_str() 
		    << " with write permission"
		    << " (error code:" << ierr << "):";
		_fileId = -1;
		return false;
	    }
	}
    }
    nc_enddef(_fileId); // leave define mode
    return true;
} // NETCDF::__openFile

bool NETCDF::__closeFile(){
    int ierr = nc_close(_fileId);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::closeFile"
            << " fail to close file " << _fileName.c_str()
	    << " (error code:" << ierr << "):";
        return false;
    }
    _fileId = -1;
    return true;
} // NETCDF::closeFile

nc_type NETCDF::__getNETCDFDataType(const FQ::DataType fqType)
{
    switch(fqType) {
    case FQ::FQT_FLOAT:
	return NC_FLOAT;
    case FQ::FQT_DOUBLE:
	return NC_DOUBLE;
    case FQ::FQT_BYTE:
	return NC_BYTE;
    case FQ::FQT_UBYTE:
	return NC_UBYTE;
    case FQ::FQT_SHORT:
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
	return -1;
    }
} // NETCDF::__getNETCDFDataType

FQ::DataType NETCDF::__getFQDataType(const nc_type type)
{
    switch(type) {
    case NC_FLOAT:
	return FQ::FQT_FLOAT;
    case NC_DOUBLE:
	return FQ::FQT_DOUBLE;
    case NC_BYTE:
    case NC_CHAR:
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
} // NETCDF::__getFQDataType

int NETCDF::__getGroupId(const std::string &variable) 
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getGroupId("
            << variable.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return -1; 
    }
    int parentId = _fileId;
    int prePos = 1;
    int pos = variable.find('/', prePos);
    while(pos != variable.npos) {
    	std::string groupName = variable.substr(prePos, pos-prePos);
	int groupId = -1;
	int ierr = nc_inq_grp_ncid(parentId, groupName.c_str(), &groupId);
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 1)
		<< "Warning -- NETCDF::__getGroupId("
		<< variable.c_str() << "):"
		<< " failed to open the group " << groupName.c_str()
		<< " (error code:" << ierr << "):";
	    return -1;
	}
	parentId = groupId;
	prePos = pos+1;
    	pos = variable.find('/', prePos);
    }
    return parentId;
} // NETCDF::__getGroupId

bool NETCDF::__getDatasetId(const std::string &variable, int *datasetId,
			    int *groupId)
{
    *groupId = __getGroupId(variable);
    if (*groupId < 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- NETCDF::__getDatasetId("
            << variable.c_str() << "): failed to open group for variable \""
	    << variable.c_str() << "\"";
        return false; 
    }
    int pos = variable.find_last_of('/');
    std::string varName = variable.substr(pos+1);
    int ierr = nc_inq_varid(*groupId, varName.c_str(), datasetId);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 1)
	    << "Warning -- NETCDF::__getDatasetId("
	    << variable.c_str() << "):"
	    << " failed to open the dataset"
	    << " (error code:" << ierr << "):";
	return false;
    }
    return true;
} // NETCDF::__getDatasetId

bool NETCDF::__getDatasetType(const std::string &variable, FQ::DataType *fqType)
{
    *fqType = FQ::FQT_UNKNOWN;

    int datasetId;
    int groupId;
    if ( ! __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getDatasetType("
	    << variable.c_str() << "):"
	    << " cannot open the dataset";
	return false;
    }

    nc_type type;
    int ierr = nc_inq_vartype(groupId, datasetId, &type);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getDatasetType("
	    << variable.c_str() << "):"
	    << " failed to get the NETCDF data type"
	    << " (error code:" << ierr << "):";
	return false;
    }
    
    switch (type) {
    case NC_FLOAT:
        *fqType = FQ::FQT_FLOAT;
	break;
    case NC_DOUBLE:
        *fqType = FQ::FQT_DOUBLE;
	break;
    case NC_CHAR:
    case NC_BYTE:
        *fqType = FQ::FQT_BYTE;
	break;
    case NC_UBYTE:
        *fqType = FQ::FQT_UBYTE;
	break;
    case NC_INT:
        *fqType = FQ::FQT_INT;
	break;
    case NC_INT64:
        *fqType = FQ::FQT_LONG;
	break;
    case NC_UINT:
        *fqType = FQ::FQT_UINT;
	break;
    case NC_UINT64:
        *fqType = FQ::FQT_ULONG;
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getDatasetType("
	    << variable.c_str() << "):"
            << " FQ does not yet support the NETCDF data type";
	return false;
    }
    return true;
} // NETCDF::__gtDatasetType

bool NETCDF::__getDatasetLength(const std::string &variable, uint64_t *len)
{
    *len = 0;
    std::vector <uint64_t> dims;
    if (! __getDatasetDimension(variable, dims)) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- NETCDF::__getDatasetLength"
	    << "(" << variable.c_str() << "):"
	    << " failed to get dataset dimension";
	return false;
    }
    *len = dims[0];
    for (unsigned int i=1; i<dims.size(); i++) {
	*len *= dims[i];
    }
    return true;
} // NETCDF::__getDatasetSize

bool NETCDF::__getDatasetDimension(const std::string &variable,
				   std::vector <uint64_t> &dims)
{
    dims.clear();
    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
	LOGGER(ibis::gVerbose > 1)
            << "Warning -- " << "NETCDF::__getDatasetDimension" 
	    << "(" << variable.c_str() << "):"
	    << " failed to open the dataset";
	return false;
    }

    int ndims;
    int ierr = nc_inq_varndims(groupId, datasetId, &ndims);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getDatasetDimension("
            << variable.c_str() << "):"
            << " failed to get number of dimension"
	    << " (error code:" << ierr << "):";
        return false;
    }

    if (ndims == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getDatasetDimension("
            << variable.c_str() << "):"
            << " number of dimension is 0";
        return false;
    }

    int dimIds[ndims];
    ierr = nc_inq_vardimid(groupId, datasetId, dimIds);
    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getDatasetDimension("
            << variable.c_str() << "):"
            << " failed to get dimension ids"
	    << " (error code:" << ierr << "):";
        return false;
    }

    for (int i=0; i<ndims; i++) {
        size_t len;
        ierr = nc_inq_dimlen(groupId, dimIds[i], &len);
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
	    	<< "Warning -- NETCDF::__getDatasetDimension("
		<< variable.c_str() << "):"
		<< " failed to get dimension length"
	        << " (error code:" << ierr << "):";
	    return false;
	}
        dims.push_back(len);
    }
    return true;
} // NETCDF::__getDatasetDimension

bool NETCDF::__getAttribute(const std::string &variable,
			    const std::string &attrName, void *values)
{
    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getAttribute("
            << variable.c_str() << "," << attrName.c_str() << "):"
            << " cannot find the variable";
        return false;
    }

    int ierr = nc_get_att(groupId, datasetId, attrName.c_str(), values);
    if (ierr != NC_NOERR ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getAttribute("
            << variable.c_str() << "," << attrName.c_str() << "):"
            << " fail to get attribute values"
	    << " (error code:" << ierr << "):";
        return false;
    }
    return true;

} // NETCDF::__getAttribute

bool NETCDF::__setAttribute(const std::string &variable,
			    const std::string &attrName, 
			    const void *values, const size_t len,
			    const FQ::DataType fqType)
{
    if (len == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__setAttribute"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " cannot create an attribute with size 0";
	return false;
    }

    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__setAttribute("
            << variable.c_str() << "," << attrName.c_str() << "):"
            << " cannot find the variable";
        return false;
    }

    nc_redef(_fileId); // enter define mode

    int attId;
    if (nc_inq_attid(groupId, datasetId, attrName.c_str(), &attId) == NC_NOERR) {
        // delete the existing attribute
        int ierr = nc_del_att(groupId, datasetId, attrName.c_str());
	if (ierr != NC_NOERR) {
	    LOGGER(ibis::gVerbose > 0)
		<< "Warning -- NETCDF::__setAttribute("
		<< variable.c_str() << "," << attrName.c_str() << "):"
		<< " failed to remove existing attribute"
		<< " (error code:" << ierr << "):";
	    nc_enddef(_fileId); // leave define mode
	    return false;
        }
    }

    // create attribute
    nc_type attrType;
    attrType = __getNETCDFDataType(fqType);
    if (attrType < 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__setAttribute"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " not yet supported FQ data type " << fqType;
	nc_enddef(_fileId); // leave define mode
	return false;
    }

    int ierr = nc_put_att(groupId, datasetId, attrName.c_str(), attrType,
			  len, values);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__setAttribute"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " failed to set attribute"
 	    << " (error code:" << ierr << "):";
	nc_enddef(_fileId); // leave define mode
	return false;
    }
    nc_enddef(_fileId); // leave define mode
    return true;
} // NETCDF::__setAttribute

bool NETCDF::__getAttributeInfo(const std::string &variable,
				const std::string &attrName,
				uint64_t *length, FQ::DataType *fqType)
{
    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getAttributeInfo("
            << variable.c_str() << "," << attrName.c_str() << "):"
            << " cannot find the variable";
        return false;
    }
    nc_type type;
    int ierr = nc_inq_atttype(groupId, datasetId, attrName.c_str(), &type);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__getAttributeInfo"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " failed to get attribute type"
 	    << " (error code:" << ierr << "):";
	return false;
    }
    *fqType = __getFQDataType(type);
    if (*fqType == FQ::FQT_UNKNOWN) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getAttributeInfo("
            << "(" << variable.c_str() << "," << attrName.c_str() << "):"
            << " failed to get the attribute type";
        return false;
    }
    size_t ltmp;
    ierr = nc_inq_attlen(groupId, datasetId, attrName.c_str(), &ltmp);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__getAttributeInfo"
	    << "(" << variable.c_str() << "," << attrName.c_str() << "):"
	    << " failed to get attribute length"
 	    << " (error code:" << ierr << "):";
	return false;
    }
    *length = ltmp;
    return true;
} // NETCDF::__getAttributeInfo

bool NETCDF::__readData(const std::string &variable, void *data)
{
    if (data == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readData("
	    << variable.c_str() << "):"
	    << " no data needs to be read";
	return false;
    }

    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readData("
	    << variable.c_str() << "):"
	    << " failed to open the dataset";
	return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    int ierr;
    switch (fqType) {
    case FQ::FQT_FLOAT:
        ierr = nc_get_var_float(groupId, datasetId, (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = nc_get_var_double(groupId, datasetId, (double*)data);
	break;
    case  FQ::FQT_BYTE:
        ierr = nc_get_var_schar(groupId, datasetId, (signed char*)data);
	break;
    case  FQ::FQT_UBYTE:
        ierr = nc_get_var_uchar(groupId, datasetId, (unsigned char*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = nc_get_var_short(groupId, datasetId, (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
        ierr = nc_get_var_ushort(groupId, datasetId, (uint16_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = nc_get_var_int(groupId, datasetId, (int32_t*)data);
	break;
    case FQ::FQT_UINT:
        ierr = nc_get_var_uint(groupId, datasetId, (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
	if (sizeof(long) == 8)
	    ierr = nc_get_var_long(groupId, datasetId, (long*)data);
	else
	    ierr = nc_get_var_longlong(groupId, datasetId, (long long*)data);
	break;
    case FQ::FQT_ULONG:
	ierr = nc_get_var_ulonglong(groupId, datasetId,
				    (unsigned long long*)data);
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the NETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readData("
            << variable.c_str() << "):"
            << " failed to read data (error code:" << ierr << "):";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::__readData("
            << variable.c_str() << "):"
            << " successfully read data";
        return true;
    }
} // NETCDF::__readData

bool NETCDF::__writeData(const std::string &variable, const void *data)
{
    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__writeData("
	    << variable.c_str() << "):"
	    << " failed to open the dataset";
	return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__writeData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    int ierr;
    switch (fqType) {
    case FQ::FQT_FLOAT:
        ierr = nc_put_var_float(groupId, datasetId, (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = nc_put_var_double(groupId, datasetId, (double*)data);
	break;
    case FQ::FQT_BYTE:
        ierr = nc_put_var_schar(groupId, datasetId, (signed char*)data);
	break;
    case FQ::FQT_UBYTE:
        ierr = nc_put_var_uchar(groupId, datasetId, (unsigned char*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = nc_put_var_short(groupId, datasetId, (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
        ierr = nc_put_var_ushort(groupId, datasetId, (uint16_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = nc_put_var_int(groupId, datasetId, (int32_t*)data);
	break;
    case FQ::FQT_UINT:
        ierr = nc_put_var_uint(groupId, datasetId, (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
	if (sizeof(long) == 8)
	    ierr = nc_put_var_long(groupId, datasetId, (long*)data);
	else
	    ierr = nc_put_var_longlong(groupId, datasetId, (long long*)data);
	break;
    case FQ::FQT_ULONG:
	ierr = nc_put_var_ulonglong(groupId, datasetId,
				    (unsigned long long*)data);
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__writeData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the NETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__writeData("
            << variable.c_str() << "):"
            << " failed to write data (error code:" << ierr << "):";
        return false;
    } else {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::__writeData("
            << variable.c_str() << "):"
            << " successfully wrote data";
        return true;
    }
} // NETCDF::__writeData

bool NETCDF::__createDataset(const std::string& variable, 
			     const std::vector<uint64_t > dims,
			     const FQ::DataType fqType)
{
    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__createDataset("
            << variable.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false; 
    }

    if (dims.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__createDataset("
	    << variable.c_str() << "):"
	    << " number of dimension is 0";
	return false;
    }

    if (variable.size() == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::__createDataset("
            << variable.c_str() << "):"
            << " variable is an empty string";
        return false;
    }

    if (variable[0] != '/') {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::__createDataset("
            << variable.c_str() << "):"
            << " variable must be given as an absolute path starting with '/'";
        return false;
    }

    //check if dataset already exist
    int datasetId;
    int groupId;
    if ( __getDatasetId(variable, &datasetId, &groupId)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::__createDataset("
            << variable.c_str() << "):"
	    << " dataset already exists and NetCDF can not re-create it";
	return false;
    } 
    // create the intermediate groups to the dataset and get parent group id
    int parentId = _fileId;
    size_t curPos = 0;
    size_t prePos = 0;
    while (1) {
    	curPos = variable.find('/', prePos+1);
	if (curPos == variable.npos) break;
	std::string groupName = variable.substr(prePos+1, curPos-prePos-1);
	int groupId;
	if (nc_inq_ncid(parentId, groupName.c_str(), &groupId) != NC_NOERR) {
	    // create the group
	    int ierr = nc_def_grp(parentId, groupName.c_str(), &groupId);
	    if (ierr != NC_NOERR) {
		LOGGER(ibis::gVerbose > 0)
		    << "Warning -- NETCDF::__createDataset("
		    << variable.c_str() << "):"
		    << " failed to create group " << groupName.c_str()
		    << "(error code:" << ierr << "):";
		nc_enddef(_fileId); // leave define mode
		return false;
	    }
	}
        parentId = groupId;
        prePos = curPos;
    }

    // get NETCDF data type
    nc_type type = __getNETCDFDataType(fqType);
    if (type < 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__createDataset("
	    << variable.c_str() << "):"
	    << " not yet supported FQ data type " << fqType;
	return false;
    }

    nc_redef(_fileId); // enter define mode

    // create dimension
    int dimIds[dims.size()];
    for (unsigned int i=0; i<dims.size(); i++) {
	// create dimension name by replacing '/' with _ and add dimension
	// rank at the end
        std::string dimName = variable;
	int pos = dimName.find('/');
	while( pos != dimName.npos ) {
	    dimName[pos] = '_'; 	
	    pos = dimName.find('/');
	}
	std::ostringstream oss;
	oss << dimName << "__" << i;
	dimName = oss.str();
        int ierr = nc_def_dim(_fileId, dimName.c_str(), dims[i], &(dimIds[i]));
	if (ierr != NC_NOERR) {
            LOGGER(ibis::gVerbose > 0)
            	<< "Warning -- NETCDF::__createDataset("
                << variable.c_str() << "):"
		<< " failed to create dimension"
		<< "(error code:" << ierr << "):" << dimName.c_str();
            nc_enddef(_fileId); // leave define mode
            return false;
        }
    }

    std::string varName = variable.substr(prePos+1);
    int ierr = nc_def_var(parentId, varName.c_str(), type, dims.size(),
			  dimIds, &datasetId);
    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__createDataset("
	    << variable.c_str() << "):"
	    << " failed to create the dataset"
	    << " (error code:" << ierr << "):";
        nc_enddef(_fileId); // leave define mode
	return false;
    }
    nc_enddef(_fileId); // leave define mode
    return true;
} // NETCDF::__createDataset

bool NETCDF::__readArrayData(const std::string& variable, 
			     const std::vector<uint64_t> &offsets, 
			     const std::vector<uint64_t> &counts, 
			     const std::vector<uint64_t> &strides, void *data)
{
    if (offsets.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data offset is empty";
	return false;
    }

    if (counts.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data count is empty";
	return false;
    }

    if (counts.size() != offsets.size()) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " size of data count and offset does not match";
	return false;
    }

    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__readArrayData "
	    << "(" << variable.c_str() << "):"
	    << " cannot open the dataset";
	return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readArrayData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    int ierr;
    switch (fqType) {
    case FQ::FQT_FLOAT:
        ierr = nc_get_vars_float(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]), (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = nc_get_vars_double(groupId, datasetId, (size_t*)(&offsets[0]),
				  (size_t*)(&counts[0]),
				  (ptrdiff_t*)(&strides[0]), (double*)data);
	break;
    case  FQ::FQT_BYTE:
        ierr = nc_get_vars_schar(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]),
				 (signed char*)data);
	break;
    case  FQ::FQT_UBYTE:
        ierr = nc_get_vars_uchar(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]),
				 (unsigned char*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = nc_get_vars_short(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]), (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
        ierr = nc_get_vars_ushort(groupId, datasetId,
				  (size_t*)(&offsets[0]),
				  (size_t*)(&counts[0]),
				  (ptrdiff_t*)(&strides[0]),
				  (uint16_t*)data);
	break;
    case FQ::FQT_INT:
        ierr = nc_get_vars_int(groupId, datasetId, (size_t*)(&offsets[0]),
			       (size_t*)(&counts[0]),
			       (ptrdiff_t*)(&strides[0]), (int32_t*)data);
	break;
    case FQ::FQT_UINT:
        ierr = nc_get_vars_uint(groupId, datasetId, (size_t*)(&offsets[0]),
				(size_t*)(&counts[0]),
				(ptrdiff_t*)(&strides[0]), (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
	if (sizeof(long) == 8)
	    ierr = nc_get_vars_long(groupId, datasetId, (size_t*)(&offsets[0]),
				    (size_t*)(&counts[0]),
				    (ptrdiff_t*)(&strides[0]),
				    (long*)data);
	else
	    ierr = nc_get_vars_longlong(groupId, datasetId,
					(size_t*)(&offsets[0]),
					(size_t*)(&counts[0]),
					(ptrdiff_t*)(&strides[0]),
					(long long*)data);
	break;
    case FQ::FQT_ULONG:
	ierr = nc_get_vars_ulonglong(groupId, datasetId,
				     (size_t*)(&offsets[0]),
				     (size_t*)(&counts[0]),
				     (ptrdiff_t*)(&strides[0]),
				     (unsigned long long*)data);
	break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readArrayData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the NETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__readArrayData("
            << variable.c_str() << "):"
            << " failed to read data (error code:" << ierr << "):";
        return false;
    } 
    return true;
} // NETCDF::__readArrayData

bool NETCDF::__writeArrayData(const std::string& variable, 
			      const std::vector<uint64_t> &offsets, 
			      const std::vector<uint64_t> &counts, 
			      const std::vector<uint64_t> &strides, 
			      const void *data)
{
    if (offsets.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data offset is empty";
	return false;
    }

    if (counts.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " data count is empty";
	return false;
    }

    if (counts.size() != offsets.size()) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " size of data count and offset does not match";
	return false;
    }

    int datasetId;
    int groupId;
    if (! __getDatasetId(variable, &datasetId, &groupId)) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " cannot open the dataset";
	return false;
    }

    FQ::DataType fqType;
    if (! __getDatasetType(variable, &fqType)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__writeArrayData("
	    << variable.c_str() << "):"
	    << " failed to get dataset type";
	return false;
    }

    int ierr;
    switch (fqType) {
    case FQ::FQT_FLOAT:
        ierr = nc_put_vars_float(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]), (float*)data);
	break;
    case FQ::FQT_DOUBLE:
        ierr = nc_put_vars_double(groupId, datasetId, (size_t*)(&offsets[0]),
				  (size_t*)(&counts[0]),
				  (ptrdiff_t*)(&strides[0]), (double*)data);
	break;
    case FQ::FQT_BYTE:
        ierr = nc_put_vars_schar(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]),
				 (signed char*)data);
	break;
    case FQ::FQT_UBYTE:
        ierr = nc_put_vars_uchar(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]),
				 (unsigned char*)data);
	break;
    case FQ::FQT_SHORT:
        ierr = nc_put_vars_short(groupId, datasetId, (size_t*)(&offsets[0]),
				 (size_t*)(&counts[0]),
				 (ptrdiff_t*)(&strides[0]),
				 (int16_t*)data);
	break;
    case FQ::FQT_USHORT:
        ierr = nc_put_vars_ushort(groupId, datasetId,
				  (size_t*)(&offsets[0]),
				  (size_t*)(&counts[0]),
				  (ptrdiff_t*)(&strides[0]),
				  (uint16_t*)data);
    case FQ::FQT_INT:
        ierr = nc_put_vars_int(groupId, datasetId, (size_t*)(&offsets[0]),
			       (size_t*)(&counts[0]),
			       (ptrdiff_t*)(&strides[0]), (int32_t*)data);
	break;
    case FQ::FQT_UINT:
        ierr = nc_put_vars_uint(groupId, datasetId, (size_t*)(&offsets[0]),
				(size_t*)(&counts[0]),
				(ptrdiff_t*)(&strides[0]), (uint32_t*)data);
	break;
    case FQ::FQT_LONG:
	if (sizeof(long) == 8)
	    ierr = nc_put_vars_long(groupId, datasetId, (size_t*)(&offsets[0]),
				    (size_t*)(&counts[0]),
				    (ptrdiff_t*)(&strides[0]),
				    (long*)data);
	else
	    ierr = nc_put_vars_longlong(groupId, datasetId,
					(size_t*)(&offsets[0]),
					(size_t*)(&counts[0]),
					(ptrdiff_t*)(&strides[0]),
					(long long*)data);
	break;
    case FQ::FQT_ULONG:
	ierr = nc_put_vars_ulonglong(groupId, datasetId,
				     (size_t*)(&offsets[0]),
				     (size_t*)(&counts[0]),
				     (ptrdiff_t*)(&strides[0]),
				     (unsigned long long*)data);
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__writeArrayData("
	    << variable.c_str() << "):"
            << " FQ does not yet support the NETCDF data type";
	return false;
    }

    if (ierr != NC_NOERR) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- NETCDF::__writeArrayData "
	    << "(" << variable.c_str() << "):"
	    << " failed to write data"
 	    << "(error code:" << ierr << "):";
	return false;
    }
    return true;
} // NETCDF::__writeArrayData

bool NETCDF::__getAllVariables(const std::string &path,
			       std::vector<std::string> &variables)
{
    variables.clear();

    if (_fileId < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- NETCDF::__getAllVariables("
            << path.c_str() << "):"
            << " file is not opened " << _fileName.c_str();
        return false; 
    }
    // if the path is not valid file location, return false
    if (path.size() == 0) {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::__getAllVariables(" 
	    << path.c_str() << "):"
 	    << " path is an empty string";
	return false;
    }

    if (path[0] != '/') {
	LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "NETCDF::__getAllVariables(" 
	    << path.c_str() << "):"
 	    << " not an absolute path";
	return false;
    }

    __traverseVariables("/", _fileId, variables);
    if (variables.size() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "NETCDF::__getAllVariables "
            << " no variable is found";
    }
    return true;
} // NETCDF::__getAllVariables

void NETCDF::__traverseVariables(const std::string &path,
				 const int parentId,
				 std::vector<std::string> &variables)
{
    int nvars;
    nc_inq_varids(parentId, &nvars, NULL);
    if (nvars > 0) {
        int varIds[nvars];
        int ierr = nc_inq_varids(parentId, NULL, varIds);
	if (ierr != NC_NOERR) {
            LOGGER(ibis::gVerbose > 0)
	    	<< "Warning -- NETCDF::__traverseVariables "
	    	<< "(" << path.c_str() << "):"
	    	<< " failed to get variable id "
 	    	<< "(error code:" << ierr << "):";
	} else {
            for (int i=0; i<nvars; i++) {
            	char varname[NC_MAX_NAME+1];
            	nc_inq_varname(parentId, varIds[i], varname);
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
    }

    int ngroups;
    nc_inq_grps(parentId, &ngroups, NULL);
    if (ngroups > 0) {
        int groupIds[ngroups];
        nc_inq_grps(parentId, NULL, groupIds);
        for (int i=0; i<ngroups; i++) {
            size_t len;
            nc_inq_grpname_full(groupIds[i], &len, NULL);
            char groupName[len+1];
            nc_inq_grpname_full(groupIds[i], NULL, groupName);
	    std::string groupNameStr = groupName;
            __traverseVariables(groupNameStr, groupIds[i], variables);
        }
    }
} // NETCDF::__traverseVariables
