// $Id$
//  Author: Jinoh Kim <jinohkim at lbl.gov>
//          John Wu <john.wu at nersc.gov>
//              Lawrence Berkeley National Laboratory
//  Copyright 2011-2013 the Regents of the University of California

#ifndef _ADIOS_Wrapper_h
#define _ADIOS_Wrapper_h

/// @file ADIOS_Wrapper.h
///
/// This header file provides definitions of wrapping classes for ADIOS
/// functions.
///

#include <vector>
#include <ibis.h>      // ibis::TYPE_T

#include <adios.h>
#include <adios_read.h>

#ifndef FQ_ADIOS_DEFAULT_BUFFER_MB
#define FQ_ADIOS_DEFAULT_BUFFER_MB 32
#endif

#ifndef FQ_ADIOS_MAX_DIM
#define FQ_ADIOS_MAX_DIM   10
#endif

/// The default time out value is set for stream operations.
#ifndef FQ_ADIOS_STREAM_TIMEOUT
#define FQ_ADIOS_STREAM_TIMEOUT 60.0
#endif

/// A singleton for invoking adios_init and adios_finalize.  It assumes the
/// MPI communicator is MPI_COMM_WORLD.  The constructor also allocates a
/// buffer of FQ_ADIOS_DEFAULT_BUFFER_MB megabytes.
class BPCommon
{
public:
    static const BPCommon& init();
    ~BPCommon();

private:
    BPCommon();
    BPCommon(const BPCommon&);
    BPCommon operator=(const BPCommon&);
};


// Forward declaration
class ADIOS_Var;

///
/// The class ADIOS_Misc provides common functions used by ADIOS-related classes.
///
class ADIOS_Misc
{
public:
    static const char* type_to_string(enum ADIOS_DATATYPES type) {
        return adios_type_to_string(type);
    }
    // static const char* value_to_string(enum ADIOS_DATATYPES type, void* data,
    //                                    int idx);
    static double convert_to_double(enum ADIOS_DATATYPES type, void* data);
    static double convert_to_double(enum ADIOS_DATATYPES type, void* data,
                                    int idx);
    static ibis::TYPE_T
    convert_adios_to_fastbit_type(ADIOS_DATATYPES adiosType);
};

///
/// The class ADIOS_File is a wrapping class for ADIOS file functions
///
class ADIOS_File
{
public:
    /// Constructor.
    ADIOS_File(const std::string &fileName,
               MPI_Comm comm        = MPI_COMM_WORLD,
               ADIOS_READ_METHOD rm = ADIOS_READ_METHOD_BP,
               float to             = 0.0,
               bool strm            = true);
    /// Destructor.
    virtual ~ADIOS_File();

    void open(const std::string &fn,
              ADIOS_READ_METHOD rm = ADIOS_READ_METHOD_BP,
              float to             = FQ_ADIOS_STREAM_TIMEOUT,
              bool strm            = false);
    void close();

    /// The MPI communicator used to open this file.
    MPI_Comm getComm() const {return m_comm;}
    /// The name of the file.
    const std::string& getName() const {return m_fname;}
    /// The ADIOS file handle.
    ADIOS_FILE*  getHandle() {return m_handle;}
    /// The current time step.  If the file handle is nil, return -1.
    int currentStep() const {return (m_handle ? m_handle->current_step : -1);}
    bool nextStep(float to = 0.0);

    ADIOS_Var* getVariable(const char *);
    ADIOS_Var* getVariable(int);

    static bool exists(const std::string& fileName);

private:
    std::string m_fname;
    MPI_Comm    m_comm;
    float       m_timeout;
    ADIOS_FILE *m_handle;

    // to prevent copying
    ADIOS_File(const ADIOS_File&);
    ADIOS_File& operator=(const ADIOS_File&);
}; // ADIOS_File

///
/// The class ADIOS_Var is a wrapping class for ADIOS variable functions
///
class ADIOS_Var
{
public:
    ADIOS_Var(ADIOS_File *, int idx);
    ADIOS_Var(ADIOS_File *, const char*);
    virtual ~ADIOS_Var();

    // Is this object valid?  Return true if yes, false otherwise.
    bool valid() {return(m_handle!=0);}

    int getIndex() const {return (m_handle ? m_handle->varid : -1);}
    const std::string& getName() const {return m_name;}
    ADIOS_VARINFO* getHandle() {return m_handle;}

    ADIOS_FILE* getFile() const {
        return m_file ? m_file->getHandle() : 0;}
    const char* getFileName() const {
        return m_file ? m_file->getHandle()->path : 0;}

    /// Get ADIOS data type
    enum ADIOS_DATATYPES getType() {
        return m_handle ? m_handle->type : adios_unknown;}
    /// Get ADIOS data type as a string
    std::string getTypeString() {return ADIOS_Misc::type_to_string(getType());}

    /// Get the number of dimension
    int64_t getNumDimension() {return m_handle ? m_handle->ndim : 0;}
    /// Get the (dim)-th dimension size
    int64_t getDimensionSize(int dim) {
        return m_handle ? m_handle->dims[dim] : 0;}
    /// Get all dimension sizes (0..ndim-1) in a vector
    int64_t getDimensionArray(std::vector<uint64_t>* vec);

    // statistics for this variable
    double getGlobalMin();
    double getGlobalMax();
    double getGlobalAvg();
    double getGlobalStd();

    /// Get data length in bytes
    int64_t getDataLen();
    /// Get the number of elements in the data
    int getNumDataElem();

    /// Read a scalar value
    void*   readValue() {return m_handle->value;}
    int64_t readData(void* buf, int64_t bufLen);
    int64_t readData(void* buf, int64_t bufLen,
                     const uint64_t* start, const uint64_t* count);
    int64_t readData(void* buf, int64_t bufLen,
                     const std::vector<uint64_t>& start,
                     const std::vector<uint64_t>& count);

protected:
    int index() {return m_handle->varid;}

private:
    std::string    m_name;
    ADIOS_File    *m_file;
    ADIOS_VARINFO *m_handle;
}; // ADIOS_Var

#endif  // _ADIOS_Wrapper_h
