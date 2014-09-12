// $Id$
//  Author: Jinoh Kim <jinohkim at lbl.gov>
//          John Wu <john.wu at nersc.gov>
//              Lawrence Berkeley National Laboratory
//  Copyright 2011-2013 the Regents of the University of California

#ifndef _BPArrayIODriver_h
#define _BPArrayIODriver_h

#include "arrayIODriver.h"
#include "ADIOS_Wrapper.h"
#include <map>	// std::map

// Forward declaration
class BPDataInfo;
class BPBitmapInfo;

/**
   The ADIOS BP implementation of Array I/O Driver.

   @note This version require ADIOS 1.4 and will not work with earlier
   versions of ADIOS.
*/
class BPArrayIODriver : public ArrayIODriver {
public:
    /*!
      Constructor.
    */
    BPArrayIODriver(const std::string& fileName,
		    const std::string& indexPrefix,
                    MPI_Comm comm        = MPI_COMM_WORLD,
                    ADIOS_READ_METHOD rm = ADIOS_READ_METHOD_BP,
                    float timeout        = 0.0,
                    bool streaming       = true);
    /*!
      Destructor.
    */
    virtual ~BPArrayIODriver();

    /*!
      Check if the file is opened successfully.
    */
    virtual bool isValid() const {
        return _adiosFile != 0;
    }

    //*********************
    // API for variables
    //*********************/
    /*!  Get the location(i.e. full-path) of all variables available
      in the file with some prefix path.

      \param path IN: A prefix path string.

      \param variables OUT: vector containing the locations(i.e. full-path)
             of each variable.

      \return True if success, otherwise return False.
    */
    virtual bool getAllVariables(const std::string &path,
                                 std::vector<std::string> &variables);

    /*!
      Retrieve information of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param dims OUT: Vector containing the dimensions for a given variable.

      \param type OUT: Enumerated type of variable

      \return True if success, otherwise return False.
    */
    virtual bool getVariableInfo(const std::string &variable,
                                 std::vector<uint64_t> &dims,
                                 FQ::DataType *type);

    /*!
      Create new dataset for a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param dims IN: Vector containing the dimensions for a given variable.

      \param type IN: Enumerated type of variable.

      \return True if success, otherwise return False.
    */
    virtual bool createDataset(const std::string& variable,
                               const std::vector<uint64_t> dims,
                               const FQ::DataType fqType);

    /*!
      Set the data values of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param data IN: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool setData(const std::string &variable,
                         const void *data);

    /*!
      Get all the data values of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param data OUT: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool getData(const std::string &variable,
                         void *data);

    /*!
      Set the data values of a subarray for a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param offsets IN: The offsets of the subarray at the each of the
             data dimensions.

      \param counts IN: The counts position of the subarray at the each of
             the data dimensions.

      \param strides IN: The strides of the subarray at the each of the
             data dimensions.

      \param data IN: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool setArrayData(const std::string &variable,
                              const std::vector<uint64_t> &offsets,
                              const std::vector<uint64_t> &counts,
                              const std::vector<uint64_t> &strides,
                              const void *data);

    /*!
      Get the data values of a subarray for a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
      the file structure.

      \param offsets IN: The offsets of the subarray at the each of the
             data dimensions.

      \param counts IN: The counts position of the subarray at the each of
             the data dimensions.

      \param strides IN: The strides of the subarray at the each of the
             data dimensions.

      \param data OUT: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool getArrayData(const std::string &variable,
                              const std::vector<uint64_t> &offsets,
                              const std::vector<uint64_t> &counts,
                              const std::vector<uint64_t> &strides,
                              void *data);
    /*!
      Get the values of a set of data selected by points from a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param coords IN: The coordinations for the data values to be
             retrieved.  All coordinates have the same dimensionality
             (rank) as the dataspace they are located within.

      \param data OUT: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool getPointData(const std::string &variable,
                              const std::vector<uint64_t> &coords,
                              void *data);

    //*********************
    // API for variable attribute
    //*********************/
    /*!
      Retrieve information of a variable attribute.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param attrNameStr IN: A attribute name of the given variable.

      \param length OUT: The number of values for the variable attribute.

      \param type OUT: Enumerated type of variable

      \return True if success, otherwise return False.
    */
    virtual bool getAttributeInfo(const std::string &variable,
                                  const std::string &attrName,
                                  uint64_t *length,
                                  FQ::DataType *type);

    /*!
      Get all the data values of a variable attribute.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param attrNameStr IN: A attribute name of the given variable.

      \param values OUT: Attribute data values. The memory size of the data
             must be reserved.

      \return True if success, otherwise return False.
    */
    virtual bool getAttribute(const std::string &variable,
                              const std::string &attrName,
                              void *values);

    /*!
      Create a vraiable attribute and set its values.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param attrNameStr IN: A attribute name of the given variable.

      \param values IN: Attribute data values. The memory size of the data
             must be reserved.

      \param length IN: The number of values for the attribute.

      \param fqType IN:  Enumerated type of variable.

      \return True if success, otherwise return False.
    */
    virtual bool setAttribute(const std::string &variable,
                              const std::string &attrName,
                              const void *values,
                              const uint64_t length,
                              const FQ::DataType fqType);
    //*********************
    // API for bitmap metadata
    //*********************/
    /*!
      Create the dataset to store bitmap keys of a variable

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param nkeys IN: The number of bitmap keys.

      \param fqType IN:  Enumerated type of variable.

      \param mpi_iter IN: The MPI iteration during the index building
             process.  It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool createBitmapKeys(const std::string &variable,
                                  const uint64_t nkeys,
                                  const FQ::DataType,
                                  const uint64_t mpi_iter,
                                  const uint64_t mpi_idx);

    /*!
      Create the dataset to store bitmap offsets of a variable

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param noffsets IN: The number of bitmap offsets.

      \param fqType IN:  Enumerated type of variable.

      \param mpi_iter IN: The MPI iteration during the index building
             process.  It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool createBitmapOffsets(const std::string &variable,
                                     const uint64_t noffsets,
                                     const FQ::DataType fqType,
                                     const uint64_t mpi_iter,
                                     const uint64_t mpi_idx);

    /*!
      Set bitmap keys of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param keys IN: The values of bitmap keys.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapKeys(const std::string& variable,
                               void *keys,
                               const uint64_t mpi_idx);

    /*!
      Retrieve bitmap offsets of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param offsets OUT: The values of bitmap offsets.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapOffsets(const std::string& variable,
                                  void *offsets,
                                  const uint64_t mpi_idx);
    /*!
      Set bitmap keys of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param keys IN: The values of bitmap keys.

      \param nkeys IN: The number of bitmap keys.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapKeys(const std::string& variable,
                               const void *keys,
                               const uint64_t nkeys,
                               const uint64_t mpi_idx);

    /*!
      Set bitmap offsets of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param offsets IN: The values of bitmap keys.

      \param noffsets IN: The number of bitmap offsets.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapOffsets(const std::string& variable,
                                  const void *offsets,
                                  const uint64_t noffsets,
                                  const uint64_t mpi_idx);

    /*!
      Get the number of bitmap keys of a variable from file.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param nkeys OUT: The number of bitmap keys.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapKeyLength(const std::string &variable,
                                    uint64_t *nkeys,
                                    const uint64_t mpi_idx);
    /*!
      Get the number of bitmap offsets of a variable from file.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param noffsets OUT: The number of bitmap offsets.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapOffsetLength(const std::string &variable,
                                       uint64_t *noffsets,
                                       const uint64_t mpi_idx);

    /*!
      Get the data type of bitmap offset of a variable from file.
      This function is needed because offsets could be stored in the type
      of either uint32_t(older version) or uint64_t(current version).

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \return the data type of bitmap offset of a variable from file.
    */
    virtual FQ::DataType getBitmapOffsetType(const std::string& variable);

    //*********************
    // API for bitmap
    //*********************/

    /*!
      Create the dataset to store bitmaps of a variable.
      The dataset must be created with the type uint32_t.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param nElements IN: The number of bitmaps (in the uint of uint32_t).

      \param mpi_iter IN: The MPI iteration during the index building process.
             It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool createBitmap(const std::string &variable,
                              const uint64_t nElements,
                              const uint64_t mpi_iter,
                              const uint64_t mpi_idx);

    /*!
      Get the length(i.e. uint32_t) of bitmaps of a variable from file.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param nkeys OUT: The number of bitmap keys.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapLength(const std::string &variable,
                                 uint64_t *len,
                                 const uint64_t mpi_idx);

    /*!
      Read a subsection of the bitmap indexes.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param startoffset IN: The start offset of the subsection.

      \param endoffset IN: The end offset of the subsection.  The
             startoffset is always *inclusive* and the endoffset is always
             *exclusive*.  The numbering of elements for the offsets starts
             from zero (as expected).

      \param data OUT: The values of bitmap indexes.  The bitmaps are
             stored in the type of uint32_t.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool readBitmap(const std::string& variable,
                            const uint64_t startoffset,
                            const uint64_t endoffset,
                            uint32_t *data,
                            const uint64_t mpi_idx);

    /*!
      Write a subsection of the bitmap indexes.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param startoffset IN: The start offset of the subsection.

      \param endoffset IN: The end offset of the subsection.  The
             startoffset is always *inclusive* and the endoffset is always
             *exclusive* The numbering of elements for the offsets starts
             from zero (as expected).

      \param data IN: The values of bitmap indexes.

      \param mpi_idx IN: The MPI index in the current iteration during the
             index building process.  It is only used under parallel MPI
             mode.

      \return True if success, otherwise return False.
    */
    virtual bool writeBitmap(const std::string& variable,
                             const uint64_t startoffset,
                             const uint64_t endoffset,
                             const uint32_t *data,
                             const uint64_t mpi_idx);

#ifndef FQ_NOMPI
    //*********************
    // API for offset table
    //*********************/
    /*!
      Set the number of bitmap keys of a variable to the offset table.
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param nkeys OUT: The number of bitmap keys.

      \param mpi_iter IN: The MPI iteration during the index building
             process.  It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapKeyLength(const std::string &variable,
                                    const uint64_t nkeys,
                                    const uint64_t mpi_iter);

    /*!
      Set the number of bitmap offsets of a variable to the offset table.
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param noffsets IN: The number of bitmap offsets.

      \param mpi_iter IN: The MPI iteration during the index building
             process.  It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapOffsetLength(const std::string &variable,
                                       const uint64_t noffsets,
                                       const uint64_t mpi_iter);

    /*!
      Set the number of bitmaps of a variable to the offset table.
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param len IN: The number of bitmaps.

      \param mpi_iter IN: The MPI iteration during the index building
             process.  It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapLength(const std::string &variable,
                                 const uint64_t len,
                                 const uint64_t mpi_iter);

    /*!
      Create the offset table
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param nElements OUT: The number of MPI iterations during the index
             building process.

      \return True if success, otherwise return False.
    */
    virtual bool createOffsetTable(const std::string &variable,
                                   const uint64_t nElements);
#endif

    //*********************
    // API for data value range
    //*********************/
    /*!
      Get the max/min value of a variable.

      We always know ActualRange are always two elements with the
      datatype the same as the indexed dataset.  So we do not need to
      query the length, nor do we need to supply the length when we
      read/write the Range data.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param range IN: The max and min values of a varaible.

      \return the max/min value.
    */
    virtual bool getActualRange(const std::string &variable,
                                void *range);

    /*!
      Get the max/min value of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param range IN: The max and min values of a varaible.

      \param fqType IN:  Enumerated type of variable.

      \return the max/min value.
    */
    virtual bool setActualRange(const std::string& variable,
                                const void *range,
                                const FQ::DataType fqType);

    /*!
      Get the max/min value of a variable.

      We always know ExpectedRange are always two elements with the
      datatype the same as the indexed dataset.  So we do not need to
      query the length, nor do we need to supply the length when we
      read/write the Range data.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param range IN: The max and min values of a varaible.

      \return the max/min value.
    */
    virtual bool getExpectedRange(const std::string &variable,
                                  void *range);

    /*!
      Get the max/min value of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in
             the file structure.

      \param range IN: The max and min values of a varaible.

      \param fqType IN:  Enumerated type of variable.

      \return the max/min value.
    */
    virtual bool setExpectedRange(const std::string& variable,
                                  const void *range,
                                  const FQ::DataType fqType);

    virtual std::string getSortedFieldName();


    /*
      Followings are ADIOS BP specific functions
    */
    /*! Is data available as a stream? */
    virtual bool isStreaming() const {return _streaming;}
    /*! Move to the next step in a stream. */
    virtual bool nextStep() const;
    /*! The current step. */
    virtual int currentStep() const {
        return (_adiosFile ? _adiosFile->currentStep() : -2);}

private:
    // member data

    /// The read-only file pointer.  The write operations will close this
    /// file before open the named file in write mode.
    ADIOS_File* _adiosFile;
    /// Index prefix (for superclass)
    std::string _indexPrefix;
    /// Is ADIOS read in streaming mode?  Ture if yes, false otherwise.
    const bool _streaming;
    /// The ADIOS read method to use.
    const ADIOS_READ_METHOD _read_method;
    /// The time out value.
    float _timeout;
    /// The MPI communicator used.
    MPI_Comm _comm;
    /// Map to maintain newly created data.  Key is the variable name.
    std::map<std::string, BPDataInfo*> _dataInfos;
    /// Map to maintain newly created FastBit index.  Key is the variable name.
    std::map<std::string, BPBitmapInfo*> _bitmapInfos;

protected:
    // access references
    std::map<std::string,BPDataInfo*>& dataInfos() { return _dataInfos; }
    std::map<std::string,BPBitmapInfo*>& bitmapInfos() { return _bitmapInfos; }

    static ibis::horometer readDataTimer;
    static ibis::horometer writeIndexTimer;

    ADIOS_Var* getADIOSVar(const std::string& variable);

    /// Convert an offset (a single dimension view) to count/start pair (a
    /// multi-dimension view used to call ADIOS read function)
    void offset2index(int64_t offset, int ndim, int64_t* dims,
                      int64_t* index);

public:
    static ibis::horometer getReadDataTimer() { return readDataTimer; }
    static ibis::horometer getWriteIndexTimer() { return writeIndexTimer; }

    // Naming conventions
    // /// Extract local variable name from variable
    // virtual std::string getVarName(const std::string& variable);
    int getGroupInfo(const std::string& variable,
                     std::string& grpName, int64_t* time);
    // /// Coin variable name from group/variable/timestep
    // virtual std::string makeVariableName
    // (const std::string& grpName, const std::string& varName, int64_t time=-1);

    bool exists(const std::string& variable);
    /// Access data type from ADIOS_Var
    FQ::DataType getType(ADIOS_Var* adiosVar) {
        return getType(adiosVar->getType());
    }
    /// Convert ADIOS_DATATYPES to FQ::DataType
    static  FQ::DataType getType(ADIOS_DATATYPES type);

    /// Access the number of dimensions
    int64_t getNumDimension(ADIOS_Var* adiosVar) {
        return adiosVar->getNumDimension();
    }
    /// Access size for a specific dimension
    int64_t getDimensionSize(ADIOS_Var* adiosVar, int dim) {
        return adiosVar->getDimensionSize(dim);
    }
    /// Access an array of sizes
    int64_t getDimensionArray(ADIOS_Var* adiosVar,
                              std::vector<uint64_t>* dims) {
        return adiosVar->getDimensionArray(dims);
    }

    /// Get data length in bytes
    int64_t getDataLen(ADIOS_Var* adiosVar) {
        return adiosVar->getDataLen();
    }
    /// Get the number of elements in the data
    int64_t getNumDataElem(ADIOS_Var* adiosVar) {
        return adiosVar->getNumDataElem();
    }

    /// Read a single scalar value
    void*   readValue(ADIOS_Var* adiosVar) {
        return adiosVar->readValue();
    }
    /// Read an array -- use bufLen=-1 for unknown buffer size
    int64_t readData(ADIOS_Var* adiosVar, int64_t time,
                     void* buf, int64_t bufLen);
    /// Read an array
    int64_t readData(ADIOS_Var* adiosVar, int64_t time,
                     void* buf, int64_t bufLen,
                     const std::vector<uint64_t>& offsets,
                     const std::vector<uint64_t>& counts);
    /// Read an array with offset (start) to offset (end)
    int64_t readData(ADIOS_Var* adiosVar,
                     void* buf, int64_t bufLen,
                     int64_t start, int64_t end);

    int64_t writeData(const std::string& variable,
                      const BPDataInfo* dinfo, const void* data);

    std::string indexGroupName
    (const std::string& variable, const int64_t mpi_idx=0);
    std::string indexFieldName
    (const std::string &indexField, const std::string &variable,
     const int64_t mpi_idx=0);
    // virtual std::string getFastBitGrpName(const std::string& indexVar);
    // virtual std::string getFastBitVarName(const std::string& indexVar);
    // //virtual int64_t getFastBitTimestep(const std::string& indexVar);
    // virtual int     getFastBitVariableInfo
    // (const std::string& indexVar, std::string& grpName,
    //  std::string& indexField, int64_t* time);

    int64_t readFastBitKeyRange(const std::string& variable,
                                double* range, const int64_t mpi_idx=0);
    // /// Read fastbit key size information
    // virtual int64_t readFastBitKeySize(const std::string& variable,
    //                                 const int64_t mpi_idx=0);
    // /// Read fastbit key array information
    // virtual int64_t readFastBitKeyArr(const std::string& variable,
    //                                void* keyArr, const int64_t mpi_idx=0);
    // /// Read fastbit offset size information
    // virtual int64_t readFastBitOffsetSize(const std::string& variable,
    //                                    const int64_t mpi_idx=0);
    // /// Read fastbit offset array information
    // virtual int64_t readFastBitOffsetArr(const std::string& variable,
    //                                   void* offsetArr,
    //                                   const int64_t mpi_idx=0);
    // /// Read fastbit bitmap size information
    // virtual int64_t readFastBitBitmapSize(const std::string& variable,
    //                                    const int64_t mpi_idx=0);
    // /// Read fastbit bitmap array
    // virtual int64_t readFastBitBitmapArr(const std::string& variable,
    //                                   uint32_t* bitmapArr,
    //                                   uint32_t startOffset,
    //                                   uint32_t endOffset,
    //                                   const int64_t mpi_idx=0);

    /// Write a fastbit index (for a column)
    int64_t writeFastBitIndex(const std::string& variable,
                              const BPBitmapInfo* binfo,
                              uint64_t startOffset,
                              uint64_t endOffset,
                              const uint64_t mpi_idx=0);
};      // BPArrayIODriver

/// Extra parameters passed through the constructor of FastQuery.  The
/// constructor of FastQuery has an optional argument for passing extra
/// arguments likes these.
struct BPExtras {
    enum ADIOS_READ_METHOD read_method;
    float timeout;
    bool streaming;

    /// The default constructor.  Suitable for reading a BP file on disk.
    BPExtras() : read_method(ADIOS_READ_METHOD_BP), timeout(0.0),
                 streaming(true) {}
}; // BPExtras

class BPDataInfo {
private:
    std::vector<uint64_t> dims;
    FQ::DataType type;

public:
    BPDataInfo(const std::vector<uint64_t>& d, FQ::DataType t)
        : dims(d), type(t) { }
    virtual ~BPDataInfo() { }

    bool valid() const {
        return dims.size() > 0;
    }

    const std::vector<uint64_t>& getDim() const { return dims; }
    FQ::DataType getType() const { return type; }

    void setDims(const std::vector<uint64_t>& d) { dims = d; }
    void setType(FQ::DataType t) { type = t; }
};    // BPDataInfo

class BPBitmapInfo {
private:
    uint64_t  keySize;
    void*     keyArr;
    FQ::DataType keyType;
    uint64_t  offsetSize;
    void*     offsetArr;
    FQ::DataType offsetType;
    uint64_t  bitmapSize;
    uint32_t* bitmapArr;
    uint64_t  bitmapEndOffset;
    uint64_t  mpi_iter;
    uint64_t  mpi_idx;

public:
    BPBitmapInfo()
        : keySize(0), keyArr(0),
          offsetSize(0), offsetArr(0),
          bitmapSize(0), bitmapArr(0),
          bitmapEndOffset(0),
          keyType(FQ::FQT_UNKNOWN), offsetType(FQ::FQT_UNKNOWN),
          mpi_iter(0), mpi_idx(0) {}
    virtual ~BPBitmapInfo() {
        if (keyArr != 0) free(keyArr, keyType);
        if (offsetArr != 0) free(offsetArr, offsetType);
        if (bitmapArr != 0) delete[] bitmapArr;
    }

    static void free(void* arr, FQ::DataType type);

    bool valid() const {
        return (keySize > 0 && offsetSize > 0 && bitmapSize > 0 &&
                keyArr != 0 && offsetArr != 0 && bitmapArr != 0);
    }

    uint64_t getKeySize() const { return keySize; }
    const void* getKeyArr() const { return keyArr; }
    FQ::DataType getKeyType() const { return keyType; }
    uint64_t getOffsetSize() const { return offsetSize; }
    const void* getOffsetArr() const { return offsetArr; }
    FQ::DataType getOffsetType() const { return offsetType; }
    uint64_t getBitmapSize() const { return bitmapSize; }
    const uint32_t* getBitmapArr() const { return bitmapArr; }
    uint64_t getMPIIter() const { return mpi_iter; }
    uint64_t getMPIIdx() const { return mpi_idx; }

    void setKeySize(uint64_t in) { keySize = in; }
    void setKeyType(FQ::DataType in) { keyType = in; }
    void setOffsetSize(uint64_t in) { offsetSize = in; }
    void setOffsetType(FQ::DataType in) { offsetType = in; }
    void setNBitmaps(uint64_t in) { bitmapSize = in; }
    void setMPIIter(uint64_t in) { mpi_iter = in; }
    void setMPIIdx(uint64_t in) { mpi_idx = in; }

    int setKeyArr(const void* in);

    int setOffsetArr(const void* in) {
        size_t sz = 0;
        if (offsetType == FQ::FQT_INT) {
            sz = offsetSize * sizeof(int32_t);
            offsetArr = new int32_t[offsetSize];
        }
        else if (offsetType == FQ::FQT_LONG) {
            sz = offsetSize * sizeof(int64_t);
            offsetArr = new int64_t[offsetSize];
        }
        else {
            LOGGER(ibis::gVerbose > 2)
                << "BPBitmapInfo::setOffsetArr only support type "
                << "FQT_INT and FQT_LONG, but got " << offsetType;
            return -1;
        }

        memcpy(offsetArr, in, sz);
        return 0;
    }

    int setBitmapArr(const uint32_t* in) {
        if (bitmapArr == 0)
            bitmapArr = new uint32_t[bitmapSize];
        memcpy(bitmapArr, in, bitmapSize*sizeof(uint32_t));
        bitmapEndOffset = bitmapSize;
        return 0;
    }

    int setBitmapArr(const uint32_t* in,
                     uint64_t startOffset, uint64_t endOffset) {
        if (bitmapArr == 0)
            bitmapArr = new uint32_t[bitmapSize];
        uint32_t* ptr = &bitmapArr[startOffset];
        memcpy(ptr, in, (endOffset-startOffset)*sizeof(uint32_t));
        bitmapEndOffset = endOffset;
        return 0;
    }

    bool setBitmapComplete() {
        if (! valid()) return false;
        return bitmapEndOffset == bitmapSize;
    }
};      // BPBitmapInfo

#endif
