// $Id$
//  Author: Jinoh Kim <jinohkim at lbl.gov>
//          John Wu <john.wu at nersc.gov>
//              Lawrence Berkeley National Laboratory
//  Copyright 2011-2013 the Regents of the University of California

#include "ADIOS_Wrapper.h"
#include <array_t.h>
#include <limits>

//
// ADIOS_File class
//

ADIOS_File::ADIOS_File(const std::string &fn, MPI_Comm comm,
                       ADIOS_READ_METHOD rm, float to, bool strm)
    : m_comm(comm), m_timeout(to), m_handle(0)
{
    LOGGER(ibis::gVerbose > 2)
        << "ADIOS_File::ctor for " << fn;
    open(fn, rm, to, strm);
}

ADIOS_File::~ADIOS_File()
{
    LOGGER(ibis::gVerbose > 2)
        << "ADIOS_File::dtor for " << m_fname;
    close();
}

/// Open the named file in read-only mode.  If a write operation is
/// desired, the caller is expected to close this file and call adios_open
/// explicitly.  This object is used to carry the file name when any write
/// operations are to be performed.  Furthermore, the file is always open
/// in streaming mode.
void ADIOS_File::open(const std::string &fn, ADIOS_READ_METHOD rm,
                      float to, bool strm)
{
    const char *cptr = (fn.empty() ? m_fname.c_str() : fn.c_str());
    if (*cptr == 0) // can not proceed with a valid file name
        return;

    if (m_handle)
        (void) adios_read_close(m_handle);

    if (strm) { // variables are visible one step at a time
        m_handle = adios_read_open_stream
            (cptr, rm, m_comm, ADIOS_LOCKMODE_CURRENT, to);
        if (m_handle != 0) { // successfully opened the stream for reading
            if (m_fname.compare(cptr) != 0)
                m_fname = cptr;
            LOGGER(ibis::gVerbose > 4)
                << "ADIOS_File::open(" << m_fname
                << ") successfully opened the name file in streaming mode";
        }
        else {
            LOGGER(ibis::gVerbose >= 0)
                << "Warning -- ADIOS_File::open(" << cptr
                << ") failed to locate an active stream";
        }
    }
    else { // open_file makes all variables visible at once
        if (m_fname.compare(cptr) != 0)
            m_fname = cptr;
        m_handle = adios_read_open_file(cptr, rm, m_comm);
        if (m_handle != 0) {
            LOGGER(ibis::gVerbose > 4)
                << "ADIOS_File::open(" << m_fname
                << ") opened the named file for reading";
        }
        else {
            LOGGER(ibis::gVerbose >= 0)
                << "Warning -- ADIOS_File::open(" << cptr
                << ") failed to open the named file in read-only mode";
        }
    }
} // ADIOS_File::open

/// Close the open file.  This function frees the content of the ADIOS_FILE
/// object under the control of this object.
void ADIOS_File::close()
{
    if (m_handle == 0) return;

    int ierr = adios_read_close(m_handle);
    if (ierr == 0) {
        LOGGER(ibis::gVerbose > 4)
            << "ADIOS_File::close successfully closed file \"" << m_fname
            << '"';
    }
    else {
        LOGGER(ibis::gVerbose >= 0)
            << "Warning -- ADIOS_File::close failed to close file \""
            << m_fname << "\" because of " << adios_errmsg();
    }
    m_handle = 0;
} // ADIOS_File::close

/// Advance to the next available step.  Returns true if successful,
/// false otherwise.
bool ADIOS_File::nextStep(float to) {
    if (m_handle == 0) return false;

    adios_release_step(m_handle);
    return (0 == adios_advance_step(m_handle, 0, to));
} // ADIOS_File::nextStep

/// Test if the file is present.  This function accomplishes its task by
/// actually attempting to open it.  It returns true if the open attempt
/// was successful, returns false otherwise.
bool ADIOS_File::exists(const std::string& fileName)
{
    ADIOS_FILE* f = adios_read_open_file
        (fileName.c_str(), ADIOS_READ_METHOD_BP, MPI_COMM_WORLD);
    if (f != 0) {
        adios_read_close(f);
        return true;
    }
    return false;
}

/// Retrieve the named variable.
ADIOS_Var* ADIOS_File::getVariable(const char* vn) {
    return new ADIOS_Var(this, vn);
} // ADIOS_File::getVariable

/// Retrieve the ith variable.
ADIOS_Var* ADIOS_File::getVariable(int i) {
    return new ADIOS_Var(this, i);
} // ADIOS_File::getVariable


//
// ADIOS_Var class
//

/// Constructor.
ADIOS_Var::ADIOS_Var(ADIOS_File* afile, int idx)
    : m_file(afile), m_handle(0)
{
    if (afile != 0 && afile->getHandle() != 0)
        m_handle = adios_inq_var_byid(afile->getHandle(), idx);
    if (m_handle !=0)
        m_name = afile->getHandle()->var_namelist[idx];
    LOGGER(ibis::gVerbose > 8)
        << "adios_inq_var_byid is called for " << m_name << ": m_handle="
        << m_handle;
}

/// Constructor.
ADIOS_Var::ADIOS_Var(ADIOS_File* afile, const char* vn)
    : m_file(afile), m_name(vn), m_handle(0)
{
    if (afile != 0 && afile->getHandle() != 0)
        m_handle = adios_inq_var(afile->getHandle(), vn);
    LOGGER(ibis::gVerbose > 8)
        << "adios_inq_var is called for " << m_name << ": m_handle="
        << m_handle;
}

/// Destructor.
ADIOS_Var::~ADIOS_Var()
{
    LOGGER(ibis::gVerbose > 8)
        << "adios_free_varinfo is called for " << m_name << ": m_handle="
        << m_handle;
    adios_free_varinfo(m_handle);
}

double
ADIOS_Var::getGlobalMin() {
    return ADIOS_Misc::convert_to_double(getType(), m_handle->statistics->min);
}

double
ADIOS_Var::getGlobalMax() {
    return ADIOS_Misc::convert_to_double(getType(), m_handle->statistics->max);
}

double
ADIOS_Var::getGlobalAvg() {
    return ADIOS_Misc::convert_to_double(getType(), m_handle->statistics->avg);
}

double
ADIOS_Var::getGlobalStd() {
    return ADIOS_Misc::convert_to_double(getType(),
                                         m_handle->statistics->std_dev);
}

int64_t
ADIOS_Var::getDimensionArray(std::vector<uint64_t>* dims)
{
    dims->resize(getNumDimension());
    for (int i=0; i<getNumDimension(); i++) {
        (*dims)[i] = getDimensionSize(i);
    }
    return dims->size();
}

int
ADIOS_Var::getNumDataElem()
{
    int getNumElem = 1;
    for (int i=0; i<getNumDimension(); i++) {
        getNumElem *= getDimensionSize(i);
    }
    return getNumElem;
}

int64_t
ADIOS_Var::getDataLen()
{
    uint64_t total_bytes = adios_type_size(getType(), readValue());
    for (int i=0; i<getNumDimension(); i++) {
        total_bytes *= getDimensionSize(i);
    }
    return total_bytes;
}

int64_t
ADIOS_Var::readData(void* buf, int64_t bufLen)
{
    ibis::array_t<uint64_t> start(m_handle->ndim, 0);
    return readData(buf, bufLen, start.begin(), m_handle->dims);
}

int64_t
ADIOS_Var::readData(void* buf, int64_t bufLen,
                    const std::vector<uint64_t>& start,
                    const std::vector<uint64_t>& count)
{
    ibis::array_t<uint64_t> s(start.size());
    ibis::array_t<uint64_t> c(start.size());

    for (int i=0; i<start.size(); i++) {
        s[i] = start[i];
        c[i] = count[i];
    }

    return readData(buf, bufLen, s.begin(), c.begin());
}

/// Read a sub-array.
///
/// @Note The length of the buffer is assumed to be sufficiently large
/// if the argument bufLen is less or equal to 0.  The caller is
/// responsible for allocating the buffer (@c buf) of the currect size.
int64_t
ADIOS_Var::readData(void* buf, int64_t bufLen,
                    const uint64_t* start,
                    const uint64_t* count)
{
    uint64_t total_elm = 1;
    for (int i=0; i<getNumDimension(); i++) {
        total_elm *= count[i];
#ifdef DEBUG
        LOGGER(ibis::gVerbose > 5)
            << "ADIOS_Var::readData: [" << i << "] start=" << start[i]
            << ", count=" << count[i];
#endif
    }
    {
        uint64_t total_bytes = total_elm *
            adios_type_size(getType(), readValue());
        if (bufLen > 0 && (uint64_t)bufLen < total_bytes) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- ADIOS_Var::readData: bufLen (" << bufLen
                << ") < total_bytes (" << total_bytes << ")";
            return -1;
        }
    }
    ADIOS_SELECTION *sel = adios_selection_boundingbox
        (m_handle->ndim, start, count);
    if (sel == 0) {
        LOGGER(ibis::gVerbose > 1)
            << "Warning -- ADIOS_Var::readData failed to create a selection";
        return -2;
    }

    IBIS_BLOCK_GUARD(adios_selection_delete, sel);
    int ierr = adios_schedule_read_byid
        (getFile(), sel, index(), getFile()->current_step, 1, buf);
    if (ierr != 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- ADIOS_Var::readData call to "
            "adios_schedule_read_byid failed due to " << adios_errmsg();
        return -3;
    }

    ierr = adios_perform_reads(getFile(), 1); // 1 == blocking
    if (ierr != 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- ADIOS_Var::readData call to adios_perform_reads on "
            << getFileName() << " failed due to " << adios_errmsg();
        return -4;
    }
    LOGGER(ibis::gVerbose > 5)
        << "ADIOS_Var::readData: competeled reading " << total_elm
        << " element" << (total_elm>1?"s":"") << " for " << getName()
        << " from " << getFileName();
    return total_elm;
} // ADIOS_Var::readData

/// Cast the incoming value into a double.
double
ADIOS_Misc::convert_to_double(enum ADIOS_DATATYPES type, void * data) {
    double ret = std::numeric_limits<double>::quiet_NaN();
    switch (type) {
    case adios_unsigned_byte:
        ret = (*((uint8_t *) data));
        break;
    case adios_byte:
        ret = (double)(*((int8_t *) data));
        break;
    case adios_short:
        ret = (double)(* ((int16_t *) data));
        break;
    case adios_unsigned_short:
        ret = (double)(* ((uint16_t *) data));
        break;
    case adios_integer:
        ret = (double)(*((int32_t *) data));
        break;
    case adios_unsigned_integer:
        ret = (double)(*((uint32_t *) data));
        break;
    case adios_long:
        ret = (double)(*((int64_t *) data));
        break;
    case adios_unsigned_long:
        ret = (double)(*((uint64_t *) data));
        break;
    case adios_real:
        ret = (double)(*((float *) data));
        break;
    case adios_double:
        ret = (double)(*((double *) data));
        break;
    case adios_long_double:
        ret = (double)(*((long double *) data));
        break;
    case adios_string:
    case adios_complex:
    case adios_double_complex:
    default:
        LOGGER(ibis::gVerbose > 0)
            << "ADIOS_Misc::convert_to_double does not support data type: "
            << (int)type;
        break;
    }
    return ret;
}


double ADIOS_Misc::convert_to_double(enum ADIOS_DATATYPES type,
                                     void * data, int idx) {
    double ret = std::numeric_limits<double>::quiet_NaN();
    switch (type) {
    case adios_unsigned_byte:
        ret = (double)(((uint8_t *) data)[idx]);
        break;

    case adios_byte:
        ret = (double)(((int8_t *) data)[idx]);
        break;

    case adios_short:
        ret = (double)( ((int16_t *) data)[idx]);
        break;

    case adios_unsigned_short:
        ret = (double)( ((uint16_t *) data)[idx]);
        break;

    case adios_integer:
        ret = (double)(((int32_t *) data)[idx]);
        break;

    case adios_unsigned_integer:
        ret = (double)(((uint32_t *) data)[idx]);
        break;

    case adios_long:
        ret = (double)(((int64_t *) data)[idx]);
        break;

    case adios_unsigned_long:
        ret = (double)(((uint64_t *) data)[idx]);
        break;

    case adios_real:
        ret = (double)(((float *) data)[idx]);
        break;

    case adios_double:
        ret = (double)(((double *) data)[idx]);
        break;

    case adios_long_double:
        ret = (double)(((long double *) data)[idx]);
        break;

    case adios_string:
    case adios_complex:
    case adios_double_complex:
    default:
        LOGGER(ibis::gVerbose > 0)
            << "ADIOS_Misc::convert_to_double does not support data type "
            << (int)type;
        break;
    }
    return ret;
}

//
// ADIOS_Misc class
//

// const char*
// ADIOS_Misc::value_to_string(enum ADIOS_DATATYPES type, void * data, int idx)
// {
//     static char s [100];
//     s [0] = 0;


//     switch (type) {
//     case adios_unsigned_byte:
//         sprintf (s, "%u", ((uint8_t *) data)[idx]);
//         break;

//     case adios_byte:
//         sprintf (s, "%d", ((int8_t *) data)[idx]);
//         break;

//     case adios_short:
//         sprintf (s, "%hd", ((int16_t *) data)[idx]);
//         break;

//     case adios_unsigned_short:
//         sprintf (s, "%hu", ((uint16_t *) data)[idx]);
//         break;

//     case adios_integer:
//         sprintf (s, "%d", ((int32_t *) data)[idx]);
//         break;

//     case adios_unsigned_integer:
//         sprintf (s, "%u", ((uint32_t *) data)[idx]);
//         break;

//     case adios_long:
//         sprintf (s, "%lld", ((long long int *) data)[idx]);
//         break;

//     case adios_unsigned_long:
//         sprintf (s, "%llu", ((long long unsigned int *) data)[idx]);
//         break;

//     case adios_real:
//         sprintf (s, "%g", ((float *) data)[idx]);
//         break;

//     case adios_double:
//         sprintf (s, "%lg", ((double *) data)[idx]);
//         break;

//     case adios_long_double:
//         sprintf (s, "%Lg", ((long double *) data)[idx]);
//         break;

//     case adios_string:
//         return (char*) ((char *)data+idx);
//         break;

//     case adios_complex:
//         sprintf (s, "(%g, %g)",
//                  ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
//         break;

//     case adios_double_complex:
//         sprintf (s, "(%lg, %lg)",
//                  ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
//         break;

//     default:
//         fprintf(stderr, "No such adios type: %d\n", (int)type);
//         exit(1);
//     }

//     return s;
// }

/// Convert ADIOS data type to fastbit data type
ibis::TYPE_T
ADIOS_Misc::convert_adios_to_fastbit_type(ADIOS_DATATYPES adios_type)
{
    ibis::TYPE_T ret = ibis::UNKNOWN_TYPE;
    switch (adios_type) {
    case adios_byte:
        ret = ibis::BYTE; break;
    case adios_unsigned_byte:
        ret = ibis::UBYTE; break;
    case adios_integer:
        ret = ibis::INT; break;
    case adios_unsigned_integer:
        ret = ibis::UINT; break;
    case adios_long:
        ret = ibis::LONG; break;
    case adios_unsigned_long:
        ret = ibis::ULONG; break;
    case adios_real:
        ret = ibis::FLOAT; break;
    case adios_double:
        ret = ibis::DOUBLE; break;
    default:
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- unsupported ADIOS data type: "
            << (int)adios_type;
        break;
    }
    return ret;
}

/// The default constructor.  It is private and can not be directly called
/// by any user code.
BPCommon::BPCommon() {
    int ierr = adios_init_noxml();
    if (ierr != 0) {
        std::cerr << "BPCommon::ctor failed to initialize adios, "
                  << adios_errmsg() << std::endl;
        throw "BPCommon failed to initialize adios";
    }

    ierr = adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW,
                                 FQ_ADIOS_DEFAULT_BUFFER_MB);
    if (ierr != 0) {
        std::cerr << "BPCommon::ctor failed to allocated "
                  << FQ_ADIOS_DEFAULT_BUFFER_MB << " MB for ADIOS buffer, "
                  << adios_errmsg() << std::endl;
    }
}

/// Destructor.
BPCommon::~BPCommon() {
    int rank=0;
    int ierr;
#ifndef FQ_NOMPI
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (ierr != MPI_SUCCESS) {
        std::cerr << "BPCommon::dtor failed to determine the MPI rank"
                  << std::endl;
        exit(-1);
    }
#endif
    ierr = adios_finalize(rank);
    if (ierr != 0) {
        std::cerr << "BPCommon::dtor failed to invoke adios_finalize, "
                  << adios_errmsg() << std::endl;
    }
}

/// Initialize the single object.  Return a refernece to the singleton
/// BPCommon object.
///
/// @note When using MPI, all processes need to call this function in order
/// for ADIOS to initialize correctly.
const BPCommon& BPCommon::init() {
    static BPCommon _ins;
    return _ins;
} // BPCommon::instance

