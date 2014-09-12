// $Id$
//  Author: Jinoh Kim <jinohkim at lbl.gov>
//          John Wu <john.wu@nersc.gov>
//              Lawrence Berkeley National Laboratory
//  Copyright 2011-2013 the Regents of the University of California

#include "BPArrayIODriver.h"
#include <memory>       // std::auto_ptr
#include <ctype.h>      // isdigit

/**
 * Variable naming conventions (NO LONGER accurate as of ADIOS 1.4)
 * For data variables:
 *      /time#/group/variable
 *              ==> group name = group
 *              ==> variable name = variable
 *              ==> timestep = time#
 *
 * For index variables:
 *      /time#/group/variable/index_field
 *              ==> group name = "/time#/group/variable"
 *              ==> variable name = index_field
 *              ==> timestep = time#
 **/

ibis::horometer BPArrayIODriver::readDataTimer;
ibis::horometer BPArrayIODriver::writeIndexTimer;

BPArrayIODriver::BPArrayIODriver(const std::string &fileName,
                                 const std::string &indexPrefix,
                                 MPI_Comm comm,
                                 ADIOS_READ_METHOD rm,
                                 float timeout,
                                 bool streaming)
    : _adiosFile(new ADIOS_File(fileName, comm, rm, timeout, streaming)),
      _indexPrefix(indexPrefix), _streaming(streaming), _read_method(rm),
      _timeout(timeout), _comm(comm)
{
    _fileName = fileName; // member of the base class
    // can not use because adios_finalize will be called after MPI_Finalize
    // BPCommon::init();
}

BPArrayIODriver::~BPArrayIODriver()
{
    delete _adiosFile;

    for (std::map<std::string,BPBitmapInfo*>::iterator it=_bitmapInfos.begin();
         it != _bitmapInfos.end(); ++it) {
        delete (*it).second;
    }
    _bitmapInfos.clear();

    for (std::map<std::string,BPDataInfo*>::iterator it=_dataInfos.begin();
         it != _dataInfos.end(); ++it) {
        delete (*it).second;
    }
    _dataInfos.clear();
}

static int tokenize(char* s, char** tok)
{
    int tokIdx = 0;
    char* p = strtok(s, "/");
    while (p != 0) {
        tok[tokIdx++] = p;
        p = strtok(0, "/");
    }

    return tokIdx;
}

bool BPArrayIODriver::getAllVariables
(const std::string& path, std::vector<std::string>& variables)
{
    if (_adiosFile == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getAllVariables("
            << path << ") can not proceed without an open file";
        return false;
    }

    const ADIOS_FILE &afile(*_adiosFile->getHandle());
    variables.resize(afile.nvars);
    for (int j = 0; j < afile.nvars; ++ j) {
        variables[j] = afile.var_namelist[j];
    }

    LOGGER(ibis::gVerbose > 8)
        << "BPArrayIODriver::getAllVariables(" << path << ") found "
        << variables.size() << " variable" << (variables.size()>1?"s":"");
    return true;
}

// std::string BPArrayIODriver::makeVariableName
// (const std::string& grpName, const std::string& varName, int64_t time)
// {
//     std::ostringstream oss;
//     oss << "/time" << time << "/" << grpName << varName;
//     return oss.str();
// }

/// Create ADIOS_Var from a variable name.
ADIOS_Var* BPArrayIODriver::getADIOSVar(const std::string& variable)
{
    if (_adiosFile == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getADIOSVar("
            << variable << ") can not proceed with an open file";
        return 0;
    }

    std::string grpName;
    int64_t time;
    getGroupInfo(variable, grpName, &time);

    return new ADIOS_Var(_adiosFile, variable.c_str());
}

// std::string BPArrayIODriver::getVarName(const std::string& variable)
// {
//     char* tok[1024];
//     char* s = strdup(variable.c_str());
//     int tokIdx = tokenize(s, tok);
//     if (tokIdx < 2) {
//      LOGGER(ibis::gVerbose >= 0)
//          << "Warning -- BPArrayIODriver::getVarName("
//          << variable << ") tokIdx < 2";
//      exit(1);
//     }

//     std::string str;
//     for (int i=2; i<tokIdx; i++) {
//      str += "/";
//      str += tok[i];
//     }
//     free(s);
//     return str;
// }

/// Extract group and time index from variable.  It assumes everyone
/// before the last delimiter '/' to be the group name.  The first segement
/// of the group starting with 'time' and followed by decimal digits is
/// assumed to contain the time index.
int BPArrayIODriver::getGroupInfo
(const std::string& variable, std::string& grpName, int64_t* time)
{
    char* tok[1024];
    char* s = strdup(variable.c_str());
    int ntok = tokenize(s, tok);
    grpName.clear();
    *time = -1;
    if (ntok > 1) { // assume timeddd to encode the time index
        grpName = tok[0];
        if (tok[0][0] == 't' && tok[0][1] == 'i' &&
            tok[0][2] == 'm' && tok[0][3] == 'e' && isdigit(tok[0][4]))
            *time = atol(tok[0]+4);

        for (int i = 1; i < ntok-1; i++) {
            if (*time < 0 && tok[i][i] == 't' && tok[i][1] == 'i' &&
                tok[i][2] == 'm' && tok[i][3] == 'e' && isdigit(tok[0][4]))
                *time = atol(tok[i]+4);

            grpName += "/";
            grpName += tok[i];
        }
    }

    free(s);
    return 0;
}

bool BPArrayIODriver::getVariableInfo
(const std::string &variable, std::vector<uint64_t> &dims, FQ::DataType *type)
{
    BPDataInfo* dinfo = dataInfos()[variable];
    if (dinfo != 0) {
        dims = dinfo->getDim();
        *type = dinfo->getType();

        LOGGER(ibis::gVerbose > 5)
            << "BPArrayIODriver::getVariableInfo(" << variable
            << ") completed (from temporal info) with  dimension = "
            << dims.size() << ", FQ type = " << *type;
        return true;
    }

    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(variable));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getVariableInfo("
            << variable << ") found no such variable";
        return false;
    }

    getDimensionArray(adiosVar.get(), &dims);
    *type = getType(adiosVar.get());
    LOGGER(ibis::gVerbose > 5)
        << "BPArrayIODriver::getVariableInfo(" << variable
        << "): completed with dimension = " << dims.size()
        << ", FQ type = " << *type;
    return true;
}

bool
BPArrayIODriver::createDataset(const std::string& variable,
                               const std::vector<uint64_t> dims,
                               const FQ::DataType fqType)
{
    BPDataInfo* dinfo = dataInfos()[variable];
    if (dinfo != 0) {
        delete dinfo;
        dataInfos().erase(variable);
    }
    dinfo = new BPDataInfo(dims, fqType);
    dataInfos()[variable] = dinfo;
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::createDataset(" << variable
        << ") complete";
    return true;
}

bool
BPArrayIODriver::setData(const std::string &variable, const void *data)
{
    BPDataInfo* dinfo = dataInfos()[variable];
    if (dinfo == 0) {
        LOGGER(ibis::gVerbose > 6)
            << "Warning -- BPArrayIODriver::setData(" << variable
            << ") no such variable to set";

        return false;
    }


    if (dinfo->valid()) {
        int64_t nbytes = writeData(variable, dinfo, data);

        LOGGER(ibis::gVerbose > 6)
            << "BPArrayIODriver::setData(" << variable
            << "): # bytes written is " << nbytes;

        delete dinfo;
        dataInfos().erase(variable);

        return nbytes > 0;
    }

    LOGGER(ibis::gVerbose > 2)
        << "Warning -- BPArrayIODriver::setData(" << variable
        << ") can not proceed due to incomplete information";
    return false;
}

bool
BPArrayIODriver::getData(const std::string &variable, void *data)
{
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(variable));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getData("
            << variable << "): no such variable";
        return false;
    }

    int64_t time = 0;

    readDataTimer.resume();
    int64_t nelm = readData(adiosVar.get(), time, data, -1);
    readDataTimer.stop();

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getData(" << variable << ""
        << "): bytes read = " << nelm;
    return nelm > 0;
}

bool
BPArrayIODriver::setArrayData(const std::string &variable,
                              const std::vector<uint64_t> &offsets,
                              const std::vector<uint64_t> &counts,
                              const std::vector<uint64_t> &strides,
                              const void *data)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::setArrayData"
        << "(" << variable << "):"
        << " not supported --- call setData instead!";

    return false;
}

bool
BPArrayIODriver::getArrayData(const std::string &variable,
                              const std::vector<uint64_t> &offsets,
                              const std::vector<uint64_t> &counts,
                              const std::vector<uint64_t> &strides,
                              void *data)
{
    for (int i=0; i<strides.size(); i++) {
        if (strides[i] > 1) {
            LOGGER(ibis::gVerbose > 1)
                << "Warning -- BPArrayIODriver::getArrayData"
                << "(" << variable << "):"
                << " not support for stride>1";
            return false;
        }
    }

    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(variable));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getArrayData(" << variable
            << ") found no such variable";
        return false;
    }

    int64_t time = 0;

    readDataTimer.resume();
    int64_t nelm = readData(adiosVar.get(), time, data, -1, offsets, counts);
    readDataTimer.stop();

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getArrayData(" << variable
        << "): # bytes read = " << nelm;
    return nelm > 0;
}

bool
BPArrayIODriver::getPointData(const std::string &variable,
                              const std::vector<uint64_t> &coords,
                              void *data)
{
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(variable));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getArrayData(" << variable
            << ") found no such variable";
        return false;
    }

    std::vector<uint64_t> dims;
    int64_t time = 0;
    if (getDimensionArray(adiosVar.get(), &dims) <= 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getPointData("
            << variable << ") failed to get dataset dimension";
        return false;
    }
    FQ::DataType fqType;
    if ((fqType = getType(adiosVar->getType())) == FQ::FQT_UNKNOWN) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getPointData("
            << variable << ") failed to get dataset type";
        return false;
    }

    unsigned int nSelects = coords.size()/dims.size();
    unsigned int nElements = dims[0];
    for(unsigned int i=1; i<dims.size(); i++) nElements *= dims[i];

    if (nSelects > (nElements >> 6)) {
        // read all then return selected points
        switch(fqType) {
        case FQ::FQT_FLOAT:{
            float* ptr = (float*)data;
            float* buf = new float[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i = 0; i < nSelects; ++ i) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim = 1; dim < dims.size(); ++dim) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_DOUBLE:{
            double* ptr = (double*)data;
            double* buf = new double[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_BYTE:{
            signed char* ptr = (signed char*)data;
            signed char* buf = new signed char[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_UBYTE:{
            unsigned char* ptr = (unsigned char*)data;
            unsigned char* buf = new unsigned char[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_SHORT: {
            int16_t* ptr = (int16_t*)data;
            int16_t* buf = new int16_t[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_USHORT:{
            uint16_t* ptr = (uint16_t*)data;
            uint16_t* buf = new uint16_t[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_INT:{
            int32_t* ptr = (int32_t*)data;
            int32_t* buf = new int32_t[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_UINT:{
            uint32_t* ptr = (uint32_t*)data;
            uint32_t* buf = new uint32_t[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_LONG:{
            int64_t* ptr = (int64_t*)data;
            int64_t* buf = new int64_t[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        case FQ::FQT_ULONG:{
            uint64_t* ptr = (uint64_t*)data;
            uint64_t* buf = new uint64_t[nElements];
            if (readData(adiosVar.get(), time, buf, -1) > 0) {
                for(unsigned int i=0; i<nSelects; i++) {
                    uint64_t pos = coords[i*dims.size()];
                    for (unsigned int dim=1; dim<dims.size(); dim++) {
                        pos = pos*dims[dim] + coords[i*dims.size()+dim];
                    }
                    if (pos > nElements) {
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning -- BPArrayIODriver::getPointData("
                            << variable << "): coords out of range";
                        return false;
                    }
                    ptr[i] = buf[pos];
                }
            }
            delete[] buf;
            break;}
        default:{
            LOGGER(ibis::gVerbose > 1)
                << "Warning -- BPArrayIODriver::getPointData("
                << variable << ") can not handle FQ data type " << fqType;
            return false;}
        }
    }
    else {
        // read point by point
        for(unsigned int i=0; i<nSelects; i++) {
            std::vector<uint64_t> offsets;
            std::vector<uint64_t> counts;
            std::vector<uint64_t> strides;
            offsets.clear();
            counts.clear();
            strides.clear();
            for (unsigned int dim=0; dim<dims.size(); dim++) {
                offsets.push_back(coords[i*dims.size()+dim]);
                counts.push_back(1);
                strides.push_back(1);
            }
            int64_t nelm=0;
            switch(fqType){
            case FQ::FQT_BYTE:{
                nelm = readData(adiosVar.get(), time,
                                &(((signed char*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_UBYTE:{
                nelm = readData(adiosVar.get(), time,
                                &(((unsigned char*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_SHORT:{
                nelm = readData(adiosVar.get(), time,
                                &(((int16_t*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_USHORT:{
                nelm = readData(adiosVar.get(), time,
                                &(((uint16_t*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_INT:{
                nelm = readData(adiosVar.get(), time,
                                &(((int32_t*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_UINT:{
                nelm = readData(adiosVar.get(), time,
                                &(((uint32_t*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_LONG:{
                nelm = readData(adiosVar.get(), time,
                                &(((int64_t*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_ULONG:{
                nelm = readData(adiosVar.get(), time,
                                &(((uint64_t*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_FLOAT:{
                nelm = readData(adiosVar.get(), time,
                                &(((float*)data)[i]), -1,
                                offsets, counts);
                break;}
            case FQ::FQT_DOUBLE:{
                nelm = readData(adiosVar.get(), time,
                                &(((double*)data)[i]), -1,
                                offsets, counts);
                break;}
            default:
                LOGGER(ibis::gVerbose > 1)
                    << "Warning -- BPArrayIODriver::getPointData("
                    << variable << "):"
                    << " unknown FQ data type " << fqType;
                return false;
            }
            if (nelm <= 0) {
                LOGGER(ibis::gVerbose > 1)
                    << "Warning -- BPArrayIODriver::getPointData("
                    << variable << "):"
                    << " failed to read array data";
                return false;
            }
        }
    }

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getPointData(" << variable
        << "): successfully got data";
    return true;
}

bool
BPArrayIODriver::getAttributeInfo(const std::string &variable,
                                  const std::string &attrName,
                                  uint64_t *length,
                                  FQ::DataType *type)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::getAttributeInfo(" << variable
        << ") is not implemented yet";
    return false;
}

bool
BPArrayIODriver::getAttribute(const std::string &variable,
                              const std::string &attrName,
                              void *values)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::getAttribute(" << variable
        << ") is no implemented yet";
    return false;
}

bool
BPArrayIODriver::setAttribute(const std::string &variable,
                              const std::string &attrName,
                              const void *values,
                              const uint64_t length,
                              const FQ::DataType fqType)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::setAttribute"
        << "(" << variable << "):"
        << " not supported!";
    return false;
}

bool
BPArrayIODriver::createBitmapKeys(const std::string &variable,
                                  const uint64_t nkeys,
                                  const FQ::DataType fqType,
                                  const uint64_t mpi_iter,
                                  const uint64_t mpi_idx)
{
    BPBitmapInfo* binfo = bitmapInfos()[variable];
    if (binfo == 0) {
        binfo = new BPBitmapInfo();
        bitmapInfos()[variable] = binfo;
    }
    binfo->setKeySize(nkeys);
    binfo->setKeyType(fqType);
    binfo->setMPIIter(mpi_iter);
    binfo->setMPIIdx(mpi_idx);

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::createBitmapKeys(" << variable
        << ") initialized the bitmap keys";

    return true;
}

bool
BPArrayIODriver::createBitmapOffsets(const std::string &variable,
                                     const uint64_t noffsets,
                                     const FQ::DataType fqType,
                                     const uint64_t mpi_iter,
                                     const uint64_t mpi_idx)
{
    BPBitmapInfo* binfo = bitmapInfos()[variable];
    if (binfo == 0) {
        binfo = new BPBitmapInfo();
        bitmapInfos()[variable] = binfo;
    }
    binfo->setOffsetSize(noffsets);
    binfo->setOffsetType(fqType);
    binfo->setMPIIter(mpi_iter);
    binfo->setMPIIdx(mpi_idx);

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::createBitmapOffsets(" << variable
        << ") initialized the bitmap offsets";

    return true;
}

bool
BPArrayIODriver::getBitmapKeys(const std::string& variable,
                               void *keys,
                               const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::getBitmapKeys(" << variable
        << ") is called...";

    std::string var = indexFieldName("keyArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::getBitmapKeys("
            << variable << ", " << var << ") found no such variable";
        return false;
    }

    int64_t nelm = readData(adiosVar.get(), -1, keys, -1);
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getBitmapKeys(" << variable << ", "
        << var << ") read " << nelm << " word" << (nelm>1?"s":"");

    return nelm > 0;
}

bool
BPArrayIODriver::getBitmapOffsets(const std::string& variable,
                                  void *offsets,
                                  const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::getBitmapOffsets(" << variable
        << ") is called...";

    std::string var = indexFieldName("offsetArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getBitmapOffsets("
            << variable << ", " << var << ") found no such variable";
        return false;
    }

    int64_t nelm = readData(adiosVar.get(), -1, offsets, -1);
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getBitmapOffsets(" << variable << ", "
        << var << "): # bytes read is " << nelm;
    return nelm > 0;
}

bool
BPArrayIODriver::setBitmapKeys(const std::string& variable,
                               const void *keys,
                               const uint64_t nkeys,
                               const uint64_t mpi_idx)
{
    BPBitmapInfo* binfo = bitmapInfos()[variable];
    if (binfo == 0) {
        binfo = new BPBitmapInfo();
        bitmapInfos()[variable] = binfo;
    }
    binfo->setKeyArr(keys);
    binfo->setKeySize(nkeys);
    binfo->setMPIIdx(mpi_idx);

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::setBitmapKeys(" << variable
        << ") recorded the bitmap keys";

    return true;
}

bool
BPArrayIODriver::setBitmapOffsets(const std::string& variable,
                                  const void *offsets,
                                  const uint64_t noffsets,
                                  const uint64_t mpi_idx)
{
    BPBitmapInfo* binfo = bitmapInfos()[variable];
    if (binfo == 0) {
        binfo = new BPBitmapInfo();
        bitmapInfos()[variable] = binfo;
    }
    binfo->setOffsetArr(offsets);
    binfo->setOffsetSize(noffsets);
    binfo->setMPIIdx(mpi_idx);

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::setBitmapOffsets(" << variable
        << ") recorded the bitmap offsets";

    return true;
}

bool
BPArrayIODriver::getBitmapKeyLength(const std::string &variable,
                                    uint64_t *nkeys,
                                    const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::getBitmapKeyLength(" << variable
        << ") is called...";

    std::string var = indexFieldName("keyArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::getBitmapKeyLength("
            << variable << ", " << var << ") found no such variable";
        return false;
    }

    *nkeys = adiosVar->getDimensionSize(0);
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getBitmapKeyLength(" << variable << ", "
        << var << ") key size is " << *nkeys;
    return (*nkeys > 0);
}

bool
BPArrayIODriver::getBitmapOffsetLength(const std::string &variable,
                                       uint64_t *noffsets,
                                       const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::getBitmapOffsetLength"
        << "(" << variable << ")"
        << " is called...";

    std::string var = indexFieldName("offsetArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::getBitmapOffsetLength("
            << variable << ", " << var << ") found no such variable";
        return false;
    }

    *noffsets = adiosVar->getDimensionSize(0);
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getBitmapOffsetLength(" << variable << ", "
        << var << "): size = " << *noffsets;
    return (*noffsets > 0);
}

FQ::DataType
BPArrayIODriver::getBitmapOffsetType(const std::string& variable)
{
    std::string var = indexFieldName("offsetArr", variable);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getBitmapOffsetType("
            << variable << ") found no such variable";
        return FQ::FQT_UNKNOWN;
    }

    FQ::DataType type = getType(adiosVar.get());
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getBitmapOffsetType(" << variable
        << ") type is " << type;
    return type;
}

bool
BPArrayIODriver::createBitmap(const std::string &variable,
                              const uint64_t nElements,
                              const uint64_t mpi_iter,
                              const uint64_t mpi_idx)
{
    BPBitmapInfo* binfo = bitmapInfos()[variable];
    if (binfo == 0) {
        binfo = new BPBitmapInfo();
        bitmapInfos()[variable] = binfo;
    }
    binfo->setNBitmaps(nElements);
    binfo->setMPIIter(mpi_iter);
    binfo->setMPIIdx(mpi_idx);

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::createBitmap(" << variable
        << ") complete";

    return true;
}

bool
BPArrayIODriver::getBitmapLength(const std::string &variable,
                                 uint64_t *len,
                                 const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::getBitmapLength(" << variable
        << ") is called...";

    std::string var = indexFieldName("bitmapArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getBitmapLength("
            << variable << ", " << var << ") found no such variable";
        return false;
    }

    *len = adiosVar->getDimensionSize(0);
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getBitmapLength(" << variable << ", "
        << var << "): size = " << *len;
    return (*len > 0);
}

bool
BPArrayIODriver::readBitmap(const std::string& variable,
                            const uint64_t startoffset,
                            const uint64_t endoffset,
                            uint32_t *data,
                            const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::readBitmap(" << variable
        << ") is called...";

    std::string var = indexFieldName("bitmapArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::readBitmap("
            << variable << ", " << var << ") found no such variable";
        return false;
    }

    const int64_t bufLen = (endoffset - startoffset) * sizeof(uint32_t);
    const int64_t nelm = readData(adiosVar.get(), data, bufLen, startoffset,
                                  endoffset);
    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::readBitmap(" << variable << ", " << var
        << ") read " << nelm << " word" << (nelm>1?"s":"")
        << ", expected " << bufLen;
    return (nelm == (endoffset - startoffset));
}

bool
BPArrayIODriver::writeBitmap(const std::string& variable,
                             const uint64_t startoffset,
                             const uint64_t endoffset,
                             const uint32_t *data,
                             const uint64_t mpi_idx)
{
    BPBitmapInfo* binfo = bitmapInfos()[variable];
    if (binfo == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "BPArrayIODriver::writeBitmap(" << variable
            << ") found no such variable to write";

        return false;
    }

    binfo->setBitmapArr(data, startoffset, endoffset);
    binfo->setMPIIdx(mpi_idx);

    if (binfo->setBitmapComplete()) {
        writeIndexTimer.resume();
        int64_t nbytes = writeFastBitIndex(variable, binfo, startoffset,
                                           endoffset, mpi_idx);
        writeIndexTimer.stop();

        LOGGER(ibis::gVerbose > 6)
            << (nbytes<0?"Warning -- ":"") << "BPArrayIODriver::writeBitmap("
            << variable << "): writeFastBitIndex returned " << nbytes;

        delete binfo;
        bitmapInfos().erase(variable);

        return nbytes > 0;
    }
    else {
        LOGGER(ibis::gVerbose > 6)
            << "Warning -- BPArrayIODriver::writeBitmap(" << variable
            << ") will not write an incomplete index";
    }

    return true;
}

#ifndef FQ_NOMPI
bool
BPArrayIODriver::setBitmapKeyLength(const std::string &variable,
                                    const uint64_t nkeys,
                                    const uint64_t mpi_iter)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::setBitmapKeyLength"
        << "(" << variable << "):"
        << " not supported -- ignore!";
    return true;
}

bool
BPArrayIODriver::setBitmapOffsetLength(const std::string &variable,
                                       const uint64_t noffsets,
                                       const uint64_t mpi_iter)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::setBitmapOffsetLength"
        << "(" << variable << "):"
        << " not supported -- ignore!";
    return true;
}

bool
BPArrayIODriver::setBitmapLength(const std::string &variable,
                                 const uint64_t len,
                                 const uint64_t mpi_iter)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::setBitmapLength"
        << "(" << variable << "):"
        << " not supported -- ignore!";
    return true;
}

bool
BPArrayIODriver::createOffsetTable(const std::string &variable,
                                   const uint64_t nElements)
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::createOffsetTable"
        << "(" << variable << "):"
        << " not supported -- ignore!";
    return true;
}
#endif

bool
BPArrayIODriver::getActualRange(const std::string &variable,
                                void *range)
{
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(variable));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::getActualRange("
            << variable << ") found no such variable";
        return false;
    }

    double v[2];
    int retVal = readFastBitKeyRange(variable, v);

    if (retVal <= 0) {
        LOGGER(ibis::gVerbose > 6)
            << "BPArrayIODriver::getActualRange(" << variable
            << "): failed to read the keys for min and max values";
        return false;
    }

    FQ::DataType type = getType(adiosVar.get());
    if (type == FQ::FQT_BYTE) {
        ((char*)range)[0] = (signed char)v[0];
        ((char*)range)[1] = (signed char)v[1];
    }
    else if (type == FQ::FQT_UBYTE) {
        ((uint64_t*)range)[0] = (unsigned char)v[0];
        ((uint64_t*)range)[1] = (unsigned char)v[1];
    }
    else if (type == FQ::FQT_SHORT) {
        ((int64_t*)range)[0] = (int16_t)v[0];
        ((int64_t*)range)[1] = (int16_t)v[1];
    }
    else if (type == FQ::FQT_USHORT) {
        ((uint64_t*)range)[0] = (uint16_t)v[0];
        ((uint64_t*)range)[1] = (uint16_t)v[1];
    }
    else if (type == FQ::FQT_INT) {
        ((int32_t*)range)[0] = (int32_t)v[0];
        ((int32_t*)range)[1] = (int32_t)v[1];
    }
    else if (type == FQ::FQT_UINT) {
        ((uint32_t*)range)[0] = (uint32_t)v[0];
        ((uint32_t*)range)[1] = (uint32_t)v[1];
    }
    else if (type == FQ::FQT_LONG) {
        ((int64_t*)range)[0] = (int64_t)v[0];
        ((int64_t*)range)[1] = (int64_t)v[1];
    }
    else if (type == FQ::FQT_ULONG) {
        ((uint64_t*)range)[0] = (uint64_t)v[0];
        ((uint64_t*)range)[1] = (uint64_t)v[1];
    }
    else if (type == FQ::FQT_FLOAT) {
        ((float*)range)[0] = (float)v[0];
        ((float*)range)[1] = (float)v[1];
    }
    else if (type == FQ::FQT_DOUBLE) {
        ((double*)range)[0] = (double)v[0];
        ((double*)range)[1] = (double)v[1];
    }
    else {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::getActualRange(" << variable
            << ") does not support this data type: " << type;
    }

    LOGGER(ibis::gVerbose > 6)
        << "BPArrayIODriver::getActualRange(" << variable << ""
        << ") complete";

    return true;
}

bool
BPArrayIODriver::setActualRange(const std::string& variable,
                                const void *range,
                                const FQ::DataType fqType)
{
    LOGGER(ibis::gVerbose > 2)
        << "Warning -- BPArrayIODriver::setActualRange(" << variable
        << ") not yet supported -- ignore!";
    return true;
}

bool
BPArrayIODriver::getExpectedRange(const std::string &variable,
                                  void *range)
{
    LOGGER(ibis::gVerbose > 2)
        << "Warning -- BPArrayIODriver::getExpectedRange("
        << variable << ") is forwarded to getActualRange";
    return getActualRange(variable, range);
}

bool
BPArrayIODriver::setExpectedRange(const std::string& variable,
                                  const void *range,
                                  const FQ::DataType fqType)
{
    LOGGER(ibis::gVerbose > 2)
        << "Warning -- BPArrayIODriver::setExpectedRange("
        << variable << ") not supported -- ignore!";
    return true;
}

std::string BPArrayIODriver::getSortedFieldName()
{
    LOGGER(ibis::gVerbose > 1)
        << "Warning -- BPArrayIODriver::getSortedFieldName:"
        << " not supported -- ignore!";
    return "";
}


/*
 * BP specific
 */

/// Test if the variable exists in the file.  This function uses
/// adios_inq_var to determine whether the named file is in the file.
bool
BPArrayIODriver::exists(const std::string& variable)
{
    if (_adiosFile == 0 || _adiosFile->getHandle() == 0) return false;

    ADIOS_VARINFO *info =
        adios_inq_var(_adiosFile->getHandle(), variable.c_str());
    bool yes = (info != 0);
    adios_free_varinfo(info);
    return yes;
}

FQ::DataType
BPArrayIODriver::getType(ADIOS_DATATYPES adios_type)
{
    if (adios_type == adios_byte)
        return FQ::FQT_BYTE;
    else if (adios_type == adios_integer)
        return FQ::FQT_INT;
    else if (adios_type == adios_unsigned_integer)
        return FQ::FQT_UINT;
    else if (adios_type == adios_long)
        return FQ::FQT_LONG;
    else if (adios_type == adios_unsigned_long)
        return FQ::FQT_ULONG;
    else if (adios_type == adios_real)
        return FQ::FQT_FLOAT;
    else if (adios_type == adios_double)
        return FQ::FQT_DOUBLE;

    LOGGER(ibis::gVerbose >= 0)
        << "Warning -- BPArrayIODriver::getType does not support ADIOS type: "
        << (int)adios_type;
    return FQ::FQT_UNKNOWN;
}


void
BPArrayIODriver::offset2index(int64_t offset, int ndim, int64_t* dims,
                              int64_t* index)
{
    int64_t p = offset;
    for (int i=ndim-1; i>=0; i--) {
        index[i] = p % dims[i];
        p /= dims[i];
    }
}

int64_t
BPArrayIODriver::readData(ADIOS_Var* adiosVar, int64_t time, void* buf,
                          int64_t bufLen)
{
    uint64_t start[FQ_ADIOS_MAX_DIM];
    uint64_t count[FQ_ADIOS_MAX_DIM];

    for (int i=0; i<FQ_ADIOS_MAX_DIM; i++) {
        start[i] = count[i] = 0;
    }

    for (int i=0; i<adiosVar->getNumDimension(); i++) {
        count[i] = adiosVar->getDimensionSize(i);
    }

    return adiosVar->readData(buf, bufLen, start, count);
}

int64_t
BPArrayIODriver::readData(ADIOS_Var* adiosVar, int64_t time, void* buf,
                          int64_t bufLen, const std::vector<uint64_t>& offsets,
                          const std::vector<uint64_t>& counts)
{
    if (time < 0) {
        uint64_t start[offsets.size()];
        uint64_t count[offsets.size()];

        for (int i=0; i<offsets.size(); i++) {
            start[i] = offsets[i];
            count[i] = counts[i];
        }
        return adiosVar->readData(buf, bufLen, start, count);
    }
    else {
        uint64_t start[offsets.size()+1];
        uint64_t count[offsets.size()+1];

        int index=0;
        for (int i=0; i<offsets.size()+1; i++) {
            start[i] = offsets[index];
            count[i] = counts[index];
            ++index;
        }

        return adiosVar->readData(buf, bufLen, start, count);
    }

}

int64_t
BPArrayIODriver::readData(ADIOS_Var* adiosVar, void* buf, int64_t bufLen,
                          int64_t startOffset, int64_t endOffset)
{
    if (startOffset > endOffset) return -1;
    // 1D only!
    int64_t ndim = adiosVar->getNumDimension();
    if (ndim > 1) {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- BPArrayIODriver::readData can only handle array "
            "with 1 dimension, ndim = " << ndim;
        return -2;
    }

    uint64_t start[1];
    uint64_t count[1];
    start[0] = startOffset;
    count[0] = endOffset - startOffset;
    return adiosVar->readData(buf, bufLen, start, count);
}

/// Writes the data to this file.  The file name is given at the time of
/// construction of this object.  This function opens the named file in
/// write mode, and closes the file prior to the return from this function.
int64_t BPArrayIODriver::writeData
(const std::string& variable, const BPDataInfo* dinfo, const void* data)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::writeData(" << variable
        << ") is called...";

    if (dinfo == 0 || !dinfo->valid()) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::writeData(" << variable
            << ") found no data to write";
        return -2;
    }

    MPI_Comm mycomm = _adiosFile->getComm();
    ADIOS_Var* adiosVar(getADIOSVar(variable));
    if (adiosVar != 0) {
        if (dinfo->getDim().size() != adiosVar->getNumDimension()) {
            LOGGER(ibis::gVerbose > 2)
                << "Warning -- BPArrayIODriver::writeData"
                << "(" << variable << "):"
                << " different # dimensions with existing one:"
                << " existing=" << adiosVar->getNumDimension()
                << " new=" << dinfo->getDim().size();
            return -3;
        }
        if (dinfo->getType() != getType(adiosVar->getType())) {
            LOGGER(ibis::gVerbose > 2)
                << "Warning -- BPArrayIODriver::writeData"
                << "(" << variable << "):"
                << " different data type with existing one:"
                << " existing=" << getType(adiosVar->getType())
                << " new=" << dinfo->getType();
            return -4;
        }
        for (int i=0; i<dinfo->getDim().size(); i++) {
            if (dinfo->getDim()[i] != adiosVar->getDimensionSize(i)) {
                LOGGER(ibis::gVerbose > 2)
                    << "Warning -- BPArrayIODriver::writeData"
                    << "(" << variable << "):"
                    << " different dimension size with existing one for "
                    << i << "th dimension:"
                    << " existing=" << adiosVar->getDimensionSize(i)
                    << " new=" << dinfo->getDim()[i];
                return -5;
            }
        }

        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::writeData"
            << "(" << variable << "):"
            << " data already exists --- update!";
    }

    // explicit file close before writing
    _adiosFile->close();

    //int64_t adios_handle;
    uint64_t adios_groupsize=0;
    uint64_t adios_totalsize=0;
    //int adios_err;

    int ierr;
    int64_t adios_file;
    int64_t adios_group;
    std::string grpName;
    int64_t time;
    getGroupInfo(variable, grpName, &time);
    ierr = adios_declare_group(&adios_group, grpName.c_str(), 0, adios_flag_yes);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeData(" << variable
        << ") -- adios_declare_group(" << grpName << ") returned "
        << ierr << ", produced group " << adios_group;
    if (adios_errno != 0) {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- adios_declare_group encounters adios_errno "
            << adios_errno << " -- " << adios_errmsg();
        return -10;
    }
#ifdef FQ_NOMPI
    const char *amethod = "POSIX";
#else
    const char *amethod = "MPI";
#endif
    ierr = adios_select_method(adios_group, amethod, "", "");
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeData(" << variable
        << ") -- adios_select_method(" << adios_group << ", " << amethod
        << ") returned " << ierr;
    if (adios_errno != 0) {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- adios_select_method encounters adios_errno "
            << adios_errno << " -- " << adios_errmsg();
        return -11;
    }

    const char* path = "";
    uint64_t array_nelem = 1;
    std::ostringstream oss, s_dim;
    oss << variable << "_D1";
    s_dim << oss.str().c_str()+1;

    adios_define_var(adios_group, oss.str().c_str()+1, path,
                     adios_unsigned_integer, 0, 0, 0);
    adios_groupsize += sizeof(uint32_t);
    array_nelem *= dinfo->getDim()[0];
    for (int i=1; i<dinfo->getDim().size(); i++) {
        std::ostringstream oss2;
        oss2 << variable << "_D" << (i+1);
        s_dim << "," << oss2.str().c_str()+1;

        adios_define_var(adios_group, oss2.str().c_str()+1, path,
                         adios_unsigned_integer, 0, 0, 0);
        adios_groupsize += sizeof(uint32_t);
        array_nelem *= dinfo->getDim()[i];
    }

    switch(dinfo->getType()){
    case FQ::FQT_BYTE:{
        adios_define_var(adios_group, variable.c_str()+1, path, adios_byte,
                         s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(char);
        break;}
    case FQ::FQT_INT:{
        adios_define_var(adios_group, variable.c_str()+1, path, adios_integer,
                         s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(int32_t);
        break;}
    case FQ::FQT_UINT:{
        adios_define_var(adios_group, variable.c_str()+1, path,
                         adios_unsigned_integer, s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(uint32_t);
        break;}
    case FQ::FQT_LONG:{
        adios_define_var(adios_group, variable.c_str()+1, path, adios_long,
                         s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(int64_t);
        break;}
    case FQ::FQT_ULONG:{
        adios_define_var(adios_group, variable.c_str()+1, path,
                         adios_unsigned_long, s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(uint64_t);
        break;}
    case FQ::FQT_FLOAT:{
        adios_define_var(adios_group, variable.c_str()+1, path, adios_real,
                         s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(float);
        break;}
    case FQ::FQT_DOUBLE:{
        adios_define_var(adios_group, variable.c_str()+1, path, adios_double,
                         s_dim.str().c_str(), 0, 0);
        adios_groupsize += array_nelem*sizeof(double);
        break;}
    default:
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::writeData("
            << variable << "):"
            << " unknown FQ data type " << dinfo->getType();
        return 0;
    }

    // "a" for append, "w" for write
    const char* option = (ADIOS_File::exists(_fileName.c_str()) ? "a" : "w");
    // open the named file for writing
    ierr = adios_open(&adios_file, grpName.c_str(), getFileName().c_str(),
                      option, &mycomm);
    if (ibis::gVerbose > 3 || ierr < 0 || adios_errno != 0) {
        ibis::util::logger lg;
        if (adios_errno != 0)
            lg() << "Warning -- ";
        lg() << "BPArrayIODriver::writeData(" << variable
             << ") -- adios_open(" << getFileName() << ", " << option
             << ") returned " << ierr << ", file descriptor = " << adios_file
             << ", adios_errno = " << adios_errno;
        if (adios_errno != 0) {
            lg() << " -- \"" << adios_errmsg() << '"';
            return -18;
        }
    }

    ierr = adios_group_size(adios_file, adios_groupsize, &adios_totalsize);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeData(" << variable
        << ") -- adios_group_size(" << adios_file << ", "
        << adios_groupsize << ") returned " << ierr;

    for (int i=0; i<dinfo->getDim().size(); i++) {
        std::ostringstream oss2;
        oss2 << variable << "_D" << (i+1);
        int index = dinfo->getDim().size()-(i+1);
        adios_write(adios_file, oss2.str().c_str(),
                    (void*)&(dinfo->getDim()[index]));
    }

    adios_write(adios_file, variable.c_str(), (void*)data);
    ierr = adios_close(adios_file);
    if (ierr < 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- BPArrayIODriver::writeData(" << variable
            << ") failed due to (" << ierr << ") " << adios_errmsg();
        return ierr;
    }
    LOGGER(ibis::gVerbose > 3)
        << "BPArrayIODriver::writeData(" << variable
        << ") completed successfully with adios_close(" << getFileName()
        << ") returning " << ierr;

    // open file for subsequent reading
    _adiosFile->open(_adiosFile->getName(), _read_method, _timeout,
                     _streaming);

    return adios_groupsize;
} // BPArrayIODriver::writeData

/// Create
std::string BPArrayIODriver::indexGroupName
(const std::string& variable, const int64_t mpi_idx)
{
    std::ostringstream oss;
    oss << "FastBitIndex";
    if (variable[0] != '/')
        oss << '/';
#ifndef FQ_NOMPI
    oss << variable << "-" << mpi_idx;
#else
    oss << variable;
#endif
    return oss.str();
}

std::string BPArrayIODriver::indexFieldName
(const std::string& indexField, const std::string& variable,
 const int64_t mpi_idx)
{
    std::ostringstream oss;
    oss << "FastBitIndex";
#ifndef FQ_NOMPI
    oss << variable << "-" << mpi_idx << "/" << indexField;
#else
    oss << variable << "/" << indexField;
#endif
    return oss.str();
}

// std::string BPArrayIODriver::getFastBitGrpName(const std::string& indexVar)
// {
//     char* tok[1024];
//     char* s = strdup(indexVar.c_str());
//     int ntok = tokenize(s, tok);
//     if (ntok < 2) {
//      LOGGER(ibis::gVerbose >= 0)
//          << "Warning -- BPArrayIODriver::getFastBitGrpName("
//          << indexVar << ") ntok < 2";
//      exit(1);
//     }

//     std::string str;
//     for (int i=0; i<ntok-1; i++) {
//      str += "/";
//      str += tok[i];
//     }
//     free(s);
//     return str;
// }

// std::string BPArrayIODriver::getFastBitVarName(const std::string& indexVar)
// {
//     char* tok[1024];
//     char* s = strdup(indexVar.c_str());
//     int ntok = tokenize(s, tok);
//     if (ntok < 2) {
//      LOGGER(ibis::gVerbose >= 0)
//          << "Warning -- BPArrayIODriver::getFastBitVarName("
//          << indexVar << ") encounters a name with ntok < 2";
//      exit(1);
//     }

//     std::string str = "/";
//     str += tok[ntok-1];
//     free(s);
//     return str;
// }

// int BPArrayIODriver::getFastBitVariableInfo
// (const std::string& indexVar, std::string& grpName,
//  std::string& indexField, int64_t* time)
// {
//     char* tok[1024];
//     char* s = strdup(indexVar.c_str());
//     int ntok = tokenize(s, tok);
//     if (ntok < 2) {
//      LOGGER(ibis::gVerbose >= 0)
//          << "Warning -- BPArrayIODriver::getFastBitVariableInfo("
//          << indexVar << ") encounters a name with ntok < 2";
//      exit(1);
//     }

//     for (int i=0; i<ntok-1; i++) {
//      grpName += "/";
//      grpName += tok[i];
//     }
//     indexField = "/";
//     indexField += tok[ntok-1];
//     *time = atol(tok[0]+4);

//     free(s);
//     return 0;
// }

// int64_t
// BPArrayIODriver::readFastBitKeySize
// (const std::string& variable, const int64_t mpi_idx)
// {
//     LOGGER(ibis::gVerbose > 7)
//      << "BPArrayIODriver::readFastBitKeySize(" << variable
//      << ") is called...";

//     std::string var = indexFieldName("keyArr", variable, mpi_idx);
//     ADIOS_Var* adiosVar = getADIOSVar(var);
//     if (adiosVar == 0) {
//      LOGGER(ibis::gVerbose > 2)
//          << "Warning -- BPArrayIODriver::readFastBitKeySize("
//          << variable << ", " << var << ") found no such variable";
//      return -1;
//     }

//     LOGGER(ibis::gVerbose > 5)
//      << "BPArrayIODriver::readFastBitKeySize("
//      << variable << ", " << var << ") key size is "
//      << adiosVar->getDimensionSize(0);

//     return adiosVar->getDimensionSize(0);
// }

// int64_t
// BPArrayIODriver::readFastBitKeyArr(const std::string& variable, void* keyArr,
//                                 const int64_t mpi_idx)
// {
//     LOGGER(ibis::gVerbose > 7)
//      << "BPArrayIODriver::readFastBitKeyArr(" << variable
//      << ") is called...";

//     std::string var = indexFieldName("keyArr", variable, mpi_idx);
//     ADIOS_Var* adiosVar = getADIOSVar(var);
//     if (adiosVar == 0) {
//      LOGGER(ibis::gVerbose > 2)
//          << "Warning -- BPArrayIODriver::readFastBitKeyArr("
//          << variable << ", " << var << ") found no such variable";
//      return -1;
//     }

//     int64_t nbytes = readData(adiosVar, -1, keyArr, -1);

//     LOGGER(ibis::gVerbose > 6)
//      << "BPArrayIODriver::readFastBitKeyArr("
//      << "(" << variable << ", " << var << ") # bytes read is " << nbytes;

//     return nbytes;
// }

// int64_t
// BPArrayIODriver::readFastBitOffsetSize(const std::string& variable,
//                                     const int64_t mpi_idx)
// {
//     LOGGER(ibis::gVerbose > 7)
//      << "BPArrayIODriver::readFastBitOffsetSize"
//      << "(" << variable << ")"
//      << " is called...";

//     std::string var = indexFieldName("offsetArr", variable, mpi_idx);
//     ADIOS_Var* adiosVar = getADIOSVar(var);
//     if (adiosVar == 0) {
//      LOGGER(ibis::gVerbose > 2)
//          << "Warning -- BPArrayIODriver::readFastBitOffsetSize("
//          << variable << ", " << var << ") found no such variable";
//      return -1;
//     }

//     LOGGER(ibis::gVerbose > 6)
//      << "BPArrayIODriver::readFastBitOffsetSize(" << variable << ", "
//      << var << "): size = " << adiosVar->getDimensionSize(0);

//     return adiosVar->getDimensionSize(0);
// }

// int64_t BPArrayIODriver::readFastBitOffsetArr
// (const std::string& variable, void* offsetArr, const int64_t mpi_idx)
// {
//     LOGGER(ibis::gVerbose > 7)
//      << "BPArrayIODriver::readFastBitOffsetArr(" << variable
//      << ") is called...";

//     std::string var = indexFieldName("offsetArr", variable, mpi_idx);
//     ADIOS_Var* adiosVar = getADIOSVar(var);
//     if (adiosVar == 0) {
//      LOGGER(ibis::gVerbose > 1)
//          << "Warning -- BPArrayIODriver::readFastBitOffsetArr("
//          << variable << ", " << var << ") found no such variable";
//      return -1;
//     }

//     int64_t nbytes = readData(adiosVar, -1, offsetArr, -1);

//     LOGGER(ibis::gVerbose > 6)
//      << "BPArrayIODriver::readFastBitOffsetArr(" << variable << ", "
//      << var << "): # bytes read is " << nbytes;

//     return nbytes;
// }

// int64_t
// BPArrayIODriver::readFastBitBitmapSize(const std::string& variable,
//                                     const int64_t mpi_idx)
// {
//     LOGGER(ibis::gVerbose > 7)
//      << "BPArrayIODriver::readFastBitBitmapSize(" << variable
//      << ") is called...";

//     std::string var = indexFieldName("bitmapArr", variable, mpi_idx);
//     ADIOS_Var* adiosVar = getADIOSVar(var);
//     if (adiosVar == 0) {
//      LOGGER(ibis::gVerbose > 1)
//          << "Warning -- BPArrayIODriver::readFastBitBitmapSize("
//          << variable << ", " << var << ") found no such variable";
//      return -1;
//     }

//     return adiosVar->getDimensionSize(0);
// }

// int64_t BPArrayIODriver::readFastBitBitmapArr
// (const std::string& variable, uint32_t* bitmapArr, uint32_t startOffset,
//  uint32_t endOffset, const int64_t mpi_idx)
// {
//     LOGGER(ibis::gVerbose > 7)
//      << "BPArrayIODriver::readFastBitBitmapArr(" << variable
//      << ") is called...";

//     std::string var = indexFieldName("bitmapArr", variable, mpi_idx);
//     ADIOS_Var* adiosVar = getADIOSVar(var);
//     if (adiosVar == 0) {
//      LOGGER(ibis::gVerbose > 1)
//          << "Warning -- BPArrayIODriver::readFastBitBitmapArr("
//          << variable << ", " << var << ") found no such variable";
//      return -1;
//     }
//     int64_t bufLen = (endOffset-startOffset+1)*sizeof(uint32_t);
//     int64_t nbytes = readData(adiosVar, bitmapArr, bufLen, startOffset,
//                            endOffset);

//     LOGGER(ibis::gVerbose > 1)
//      << "Warning -- BPArrayIODriver::readFastBitBitmapArr"
//      << "(" << variable << ", " << var << "):"
//      << " # bytes read is " << nbytes;

//     if (nbytes <= 0) {
//      return -1;
//     }

//     return nbytes;
// }

/// Read FastBit key range.  The range is expressed at two doubles, where
/// the first one is the minimum and the second one is the maximum.
int64_t
BPArrayIODriver::readFastBitKeyRange(const std::string& variable,
                                     double* range,
                                     const int64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::readFastBitKeyRange(" << variable
        << ") is called...";


    std::string var = indexFieldName("keyArr", variable, mpi_idx);
    std::auto_ptr<ADIOS_Var> adiosVar(getADIOSVar(var));
    if (adiosVar.get() == 0) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::readFastBitKeyRange("
            << variable << ", " << var << ") found no such variable";
        return -1;
    }

    range[0] = adiosVar->getGlobalMin();
    range[1] = adiosVar->getGlobalMax();

    LOGGER(ibis::gVerbose > 5)
        << "BPArrayIODriver::readFastBitKeyRange(" << variable << ", "
        << var << ") range = [" << range[0] << ", " << range[1] << "]";
    return 0;
}

/// Write a FastBit index to file.  This function actually opens the named
/// file and output the content stored in BPBitmapInfo object.  Note that
/// the BPBitmapInfo object must be properly initialized.
int64_t BPArrayIODriver::writeFastBitIndex
(const std::string& variable, const BPBitmapInfo* binfo,
 uint64_t startOffset, uint64_t endOffset, const uint64_t mpi_idx)
{
    LOGGER(ibis::gVerbose > 7)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") is called...";

    if (binfo == 0 || !binfo->valid()) {
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPArrayIODriver::writeFastBitIndex"
            << "(" << variable << "):"
            << " index is not valid!";
        return -2;
    }

    // std::string var = indexFieldName("bitmapArr", variable, mpi_idx);
    // ADIOS_Var* adiosVar = getADIOSVar(var);

    // close the read-only file pointer
    _adiosFile->close();

    //int64_t adios_handle;
    uint64_t adios_groupsize=0;
    uint64_t adios_totalsize=0;
    //int adios_err;

    int64_t adios_file;
    int64_t adios_group;

    std::string grpName = indexGroupName(variable, mpi_idx);
    int ierr = adios_declare_group(&adios_group,
                                   grpName.c_str()+(grpName[0]=='/'),
                                   0, adios_flag_no);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_declare_group(" << grpName << ") returned "
        << ierr << ", produced group " << adios_group;
    if (adios_errno != 0) {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- adios_errno = " << adios_errno << " -- "
            << adios_errmsg();
        return -10;
    }
#ifdef FQ_NOMPI
    const char *amethod = "POSIX";
#else
    const char *amethod = "MPI";
#endif
    ierr = adios_select_method(adios_group, amethod, "", "");
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_select_method(" << adios_group << ", " << amethod
        << ") returned " << ierr;
    if (adios_errno != 0) {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- adios_select_method encounters adios_errno "
            << adios_errno << " -- " << adios_errmsg();
        return -11;
    }

    const char *str1, *str2;
    const char *path = "";
    std::string s_bitmapSize =
        indexFieldName("bitmapSize", variable, mpi_idx);
    std::string s_bitmapArr =
        indexFieldName("bitmapArr", variable, mpi_idx);
    std::string s_offsetSize =
        indexFieldName("offsetSize", variable, mpi_idx);
    std::string s_offsetArr =
        indexFieldName("offsetArr", variable, mpi_idx);
    std::string s_keySize =
        indexFieldName("keySize", variable, mpi_idx);
    std::string s_keyArr =
        indexFieldName("keyArr", variable, mpi_idx);
    MPI_Comm mycomm = _adiosFile->getComm();

    // skip the leading '/' in the variable names
    str1 = s_bitmapSize.c_str()+(s_bitmapSize[0]=='/');
    ierr = adios_define_var
        (adios_group, str1, path, adios_unsigned_integer, 0, 0, 0);
    adios_groupsize += sizeof(uint32_t);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_define_var(" << str1 << ") returned " << ierr;
    if (ierr < 0) return -12;

    str2 = s_bitmapArr.c_str()+(s_bitmapArr[0]=='/');
    ierr = adios_define_var
        (adios_group, str2, path, adios_unsigned_integer, str1, 0, 0);
    adios_groupsize += binfo->getBitmapSize()*sizeof(uint32_t);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_define_var(" << str2 << ") returned "
        << ierr;
    if (ierr < 0) return -13;

    str1 = s_offsetSize.c_str()+(s_offsetSize[0]=='/');
    ierr = adios_define_var
        (adios_group, str1, path, adios_unsigned_integer, 0, 0, 0);
    adios_groupsize += sizeof(uint32_t);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_define_var(" << str1 << ") returned " << ierr;
    if (ierr < 0) return -14;

    str2 = s_offsetArr.c_str()+(s_offsetArr[0]=='/');
    if (binfo->getOffsetType() == FQ::FQT_LONG) {
        ierr = adios_define_var
            (adios_group, str2, path, adios_unsigned_long, str1, 0, 0);
        adios_groupsize += binfo->getOffsetSize()*sizeof(uint64_t);
    }
    else {
        ierr = adios_define_var
            (adios_group, str2, path, adios_unsigned_integer, str1, 0, 0);
        adios_groupsize += binfo->getOffsetSize()*sizeof(uint32_t);
    }
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_define_var(" << str2 << ") returned " << ierr;
    if (ierr < 0) return -15;

    str1 = s_keySize.c_str()+(s_keySize[0]=='/');
    ierr = adios_define_var
        (adios_group, str1, path, adios_unsigned_integer, 0, 0, 0);
    adios_groupsize += sizeof(uint32_t);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_define_var(" << str1 << ") returned "
        << ierr;
    if (ierr < 0) return -16;

    str2 = s_keyArr.c_str()+(s_keyArr[0]=='/');
    ierr = adios_define_var
        (adios_group, str2,
         path, adios_double, str1, 0, 0);
    adios_groupsize += binfo->getKeySize()*sizeof(double);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_define_var(" << str2 << ") returned " << ierr;
    if (ierr < 0) return -17;

    // "a" for append, "w" for write
    const char* option = (ADIOS_File::exists(_fileName.c_str()) ? "a" : "w");
    // need to initialize for the write operation
    adios_init_noxml();
    // open the named file for writing
    ierr = adios_open(&adios_file, grpName.c_str(), _fileName.c_str(),
                      option, &mycomm);
    if (ibis::gVerbose > 3 || ierr < 0 || adios_errno != 0) {
        ibis::util::logger lg;
        lg() << "BPArrayIODriver::writeFastBitIndex(" << variable
             << ") -- adios_open(" << grpName << ", " << _fileName << ", "
             << option << ") returned " << ierr << ", file descriptor = "
             << adios_file << ", adios_errno = " << adios_errno;
        if (adios_errno != 0) {
            lg() << " -- \"" << adios_errmsg() << '"';
            //return -18;
        }
    }
    ierr = adios_group_size(adios_file, adios_groupsize, &adios_totalsize);
    LOGGER(ibis::gVerbose > 4 || ierr < 0)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") -- adios_group_size(" << adios_file << ", "
        << adios_groupsize << ") returned " << ierr;

    uint32_t bitmapSize = binfo->getBitmapSize();
    const uint32_t* bitmapArr = binfo->getBitmapArr();
    uint32_t offsetSize = binfo->getOffsetSize();
    const uint32_t* offsetArr = (uint32_t*) binfo->getOffsetArr();
    uint32_t keySize = binfo->getKeySize();
    const double* keyArr = (double*) binfo->getKeyArr();

    if (bitmapSize > 0 && offsetSize > 0 && keySize > 0) {
        ierr = adios_write(adios_file, s_bitmapSize.c_str(), &bitmapSize);
        LOGGER(ibis::gVerbose > 5 || adios_errno != 0)
            << "BPArrayIODriver::writeFastBitIndex(" << variable
            << ") -- adios_write(s_bitmapSize) returned " << adios_errno;
        //if (adios_errno != 0) return -19;
        adios_write(adios_file, s_bitmapArr.c_str(), (void*)bitmapArr);
        adios_write(adios_file, s_offsetSize.c_str(), &offsetSize);
        adios_write(adios_file, s_offsetArr.c_str(), (void*)offsetArr);
        adios_write(adios_file, s_keySize.c_str(), &keySize);
        adios_write(adios_file, s_keyArr.c_str(), (void*)keyArr);
    }
    ierr = adios_close(adios_file);
    if (adios_errno != 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- BPArrayIODriver::writeFastBitIndex(" << variable
            << ") failed due to (" << adios_errno << ") " << adios_errmsg();
        //return -20;
    }

    // need to finalize to avoid the group definition to interfere with
    // future writing operations
    int mpi_rank;
    MPI_Comm_rank(_comm, &mpi_rank);
    ierr = adios_finalize(mpi_rank);
    LOGGER(ibis::gVerbose > 3)
        << "BPArrayIODriver::writeFastBitIndex(" << variable
        << ") completed successfully with adios_close(" << getFileName()
        << ") returning " << ierr;

    // open the file for subsequent reading
    _adiosFile->open(_adiosFile->getName(), _read_method, _timeout,
                     _streaming);

    return adios_groupsize;
} // BPArrayIODriver::writeFastBitIndex

/// Returns true only if the function adios_advance_step completes without
/// any error.  Note that the error checking is based on adios_errno, not
/// the return value of the function adios_advance_step.
bool BPArrayIODriver::nextStep() const {
    bool res = _adiosFile->nextStep(_timeout);
    if (adios_errno == err_no_error) {
        LOGGER(ibis::gVerbose > 3)
            << "BPArrayIODriver::nextStep advances to step "
            << _adiosFile->currentStep() << " in " << _adiosFile->getName();
    }
    else {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- BPArrayIODriver::nextStep encountered "
            << adios_errmsg() << " when attempting to move to the next step";
    }
    return res;
} // BPArrayIODriver::nextStep

void BPBitmapInfo::free(void* arr, FQ::DataType type) {
    switch (type) {
    case FQ::FQT_BYTE:
        delete[] (signed char*)arr;
        break;
    case FQ::FQT_UBYTE:
        delete[] (unsigned char*)arr;
        break;
    case FQ::FQT_SHORT:
        delete[] (int16_t*)arr;
        break;
    case FQ::FQT_USHORT:
        delete[] (uint16_t*)arr;
        break;
    case FQ::FQT_INT:
        delete[] (int32_t*)arr;
        break;
    case FQ::FQT_UINT:
        delete[] (uint32_t*)arr;
        break;
    case FQ::FQT_LONG:
        delete[] (int64_t*)arr;
        break;
    case FQ::FQT_ULONG:
        delete[] (uint64_t*)arr;
        break;
    case FQ::FQT_FLOAT:
        delete[] (float*)arr;
        break;
    case FQ::FQT_DOUBLE:
        delete[] (double*)arr;
        break;
    default:
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPBitmapInfo::free can not handle type " << type;
    }
} // BPBitmapInfo::free


int BPBitmapInfo::setKeyArr(const void* in) {
    size_t sz = -1;
    switch (keyType) {
    case FQ::FQT_BYTE:
        sz = keySize * sizeof(signed char);
        keyArr = new signed char[keySize];
        break;
    case FQ::FQT_UBYTE:
        sz = keySize * sizeof(unsigned char);
        keyArr = new unsigned char[keySize];
        break;
    case FQ::FQT_SHORT:
        sz = keySize * sizeof(int16_t);
        keyArr = new int16_t[keySize];
        break;
    case FQ::FQT_USHORT:
        sz = keySize * sizeof(uint16_t);
        keyArr = new uint16_t[keySize];
        break;
    case FQ::FQT_INT:
        sz = keySize * sizeof(int32_t);
        keyArr = new int32_t[keySize];
        break;
    case FQ::FQT_UINT:
        sz = keySize * sizeof(uint32_t);
        keyArr = new uint32_t[keySize];
        break;
    case FQ::FQT_LONG:
        sz = keySize * sizeof(int64_t);
        keyArr = new int64_t[keySize];
        break;
    case FQ::FQT_ULONG:
        sz = keySize * sizeof(uint64_t);
        keyArr = new uint64_t[keySize];
        break;
    case FQ::FQT_FLOAT:
        sz = keySize * sizeof(float);
        keyArr = new float[keySize];
        break;
    case FQ::FQT_DOUBLE:
        sz = keySize * sizeof(double);
        keyArr = new double[keySize];
        break;
    default:
        LOGGER(ibis::gVerbose > 2)
            << "Warning -- BPBitmapInfo::setKeyArr can not handle type "
            << keyType;
        return -1;
    }

    memcpy(keyArr, in, sz);
    return 0;
}
