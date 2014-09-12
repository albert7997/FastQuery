#ifndef _HDF5FILE_H
#define _HDF5FILE_H

#define HDF5_MPI_POSIX
#define HDF5_DISABLE_METADATA_FLUSH
#define HDF5_DATA_ALIGNMENT

#define HDF5_COLLECTIVE_WRITE false

#include <vector>
#include <string>
#include "arrayIODriver.h"
#include <hdf5.h>

/*
 * Define operator data structure type for H5Literate callback.
 * During recursive iteration, these structures will form a
 * linked list that can be searched for duplicate groups,
 * preventing infinite recursion.
 */
class Node {
public:
    Node(Node* node, const haddr_t addr, const char *name, bool dataset) {
        parent = node;
	this->addr = addr;
	this->name = name;
	this->dataset = dataset;
	children.clear();
    }
    ~Node() {
	for (unsigned int i = 0; i < getNumChildren(); ++ i) {
	    Node* childNode = getChild(i);
	    delete(childNode);
	}
    }
    void addChild(Node* node) {
	children.push_back(node);
    }
    unsigned int getNumChildren() {return children.size();}
    Node* getChild(const unsigned int idx) {return children[idx];}
    Node* getParent() {return parent;}
    haddr_t getAddr() {return addr;}
    std::string getName() {return name;}
    bool isDataset() {return dataset;}
private:
    Node *parent;                 /* Pointer to previous opdata */
    std::vector<Node*> children;  /* Pointer to previous opdata */
    haddr_t addr;	          /* Group address */
    std::string name;
    bool dataset;
};
/*
 * Operator function to be called by H5Literate.
 */
herr_t op_func(hid_t loc_id, const char *name, const H5L_info_t *info,
	       void *operator_data);

/*
 * Function to check for duplicate groups in a path.
 */
int group_check(Node* node, haddr_t target_addr);

class HDF5 : public ArrayIODriver {
public: 
#ifdef FQ_NOMPI
    HDF5(const std::string fileName,
	 bool readOnly=true,
	 const std::string indexPath="");
#else
    HDF5(const std::string fileName, 
	 MPI_Comm comm,
	 bool readOnly=true, 
	 const std::string indexPath="");
#endif
    virtual ~HDF5();

    /*!
      \brief Check if the file is opened successfully.
    */
    virtual bool isValid() const {return (_fileId<0?false:true);}

    //*********************
    // API for variables
    //*********************/
    /*!
      \brief Get the location(i.e. full-path) of all variables available in the file.

      \param path IN: A prefix path string.
      \param variables OUT: vector containing the locations(i.e. full-path) of each variable.

      \return True if success, otherwise return False.
    */
    virtual bool getAllVariables(const std::string &path,
                                 std::vector<std::string> &variables);

    /*!
      \brief Retrieve information of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param dims OUT: Vector containing the dimensions for a given variable.
      \param type OUT: Enumerated type of variable

      \return True if success, otherwise return False.
    */
    virtual bool getVariableInfo(const std::string &variable,
                                 std::vector <uint64_t> &dims,
                                 FQ::DataType *type);

    /*!
      \brief Create new dataset for a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param dims IN: Vector containing the dimensions for a given variable.
      \param type IN: Enumerated type of variable.

      \return True if success, otherwise return False.
    */
    virtual bool createDataset(const std::string& variable,
                               const std::vector<uint64_t> dims,
                               const FQ::DataType fqType);

    /*!
      \brief Set the data values of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param data IN: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool setData(const std::string &variable,
                         const void *data);

    /*!
      \brief Get all the data values of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param data OUT: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool getData (const std::string &variable,
			  void *data);

    /*!
      \brief Set the data values of a subarray for a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param offsets IN: The offsets of the subarray at the each of the data dimensions.
      \param counts IN: The counts position of the subarray at the each of the data dimensions.
      \param strides IN: The strides of the subarray at the each of the data dimensions.
      \param data IN: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool setArrayData(const std::string &variable,
                              const std::vector<uint64_t> &offsets,
                              const std::vector<uint64_t> &counts,
                              const std::vector<uint64_t> &strides,
                              const void *data);
    /*!
      \brief Get the data values of a subarray for a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param offsets IN: The offsets of the subarray at the each of the data dimensions.
      \param counts IN: The counts position of the subarray at the each of the data dimensions.
      \param strides IN: The strides of the subarray at the each of the data dimensions.
      \param data OUT: Data values.

      \return True if success, otherwise return False.
    */
    virtual bool getArrayData(const std::string &variable,
                              const std::vector<uint64_t> &offsets,
                              const std::vector<uint64_t> &counts,
                              const std::vector<uint64_t> &strides,
                              void *data);
    /*!
      \brief Get the values of a set of data selected by points from a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param coords IN: The coordinations for the data values to be retrieved.
      All coordinates have the same dimensionality (rank) 
      as the dataspace they are located within. 
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
      \brief Retrieve information of a variable attribute.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
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
      \brief Get all the data values of a variable attribute.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param attrNameStr IN: A attribute name of the given variable.
      \param values OUT: Attribute data values. The memory size of the data must be reserved.

      \return True if success, otherwise return False.
    */
    virtual bool getAttribute(const std::string &variable,
                              const std::string &attrName,
                              void *values);

    /*!
      \brief Create a vraiable attribute and set its values.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param attrNameStr IN: A attribute name of the given variable.
      \param values IN: Attribute data values. The memory size of the data must be reserved.
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
      \brief Create the dataset to store bitmap keys of a variable 

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param nkeys IN: The number of bitmap keys.
      \param fqType IN:  Enumerated type of variable.
      \param mpi_iter IN: The MPI iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool createBitmapKeys(const std::string &variable,
                                  const uint64_t nkeys,
                                  const FQ::DataType,
                                  const uint64_t mpi_iter,
				  const uint64_t mpi_idx);

    /*!
      \brief Create the dataset to store bitmap offsets of a variable 

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param noffsets IN: The number of bitmap offsets.
      \param fqType IN:  Enumerated type of variable.
      \param mpi_iter IN: The MPI iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool createBitmapOffsets(const std::string &variable,
                                     const uint64_t noffsets,
                                     const FQ::DataType fqType,
                                     const uint64_t mpi_iter,
				     const uint64_t mpi_idx);

    /*!
      \brief Retrieve bitmap keys of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param keys OUT: The values of bitmap keys.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapKeys(const std::string& variable,
                               void *keys,
                               const uint64_t mpi_idx);

    /*!
      \brief Retrieve bitmap offsets of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param offsets OUT: The values of bitmap offsets.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapOffsets(const std::string& variable,
                                  void *offsets,
                                  const uint64_t mpi_idx);

    /*!
      \brief Set bitmap keys of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param keys IN: The values of bitmap keys.
      \param nkeys IN: The number of bitmap keys.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapKeys(const std::string& variable,
                               const void *keys,
                               const uint64_t nkeys,
                               const uint64_t mpi_idx);
    /*!
      \brief Set bitmap offsets of a variable to file

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param keys IN: The values of bitmap keys.
      \param noffsets IN: The number of bitmap offsets.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapOffsets(const std::string& variable,
                                  const void *offsets,
                                  const uint64_t noffsets,
                                  const uint64_t mpi_idx);
    /*!
      \brief Get the number of bitmap keys of a variable from file.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param nkeys OUT: The number of bitmap keys.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapKeyLength(const std::string &variable,
                                    uint64_t *nkeys,
                                    const uint64_t mpi_idx);

    /*!
      \brief Get the number of bitmap offsets of a variable from file.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param noffsets OUT: The number of bitmap offsets.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapOffsetLength(const std::string &variable,
                                       uint64_t *noffsets,
                                       const uint64_t mpi_idx);

    /*!
      \brief Get the data type of bitmap offset of a variable from file.
      This function is needed because offsets could be stored 
      in the type of either uint32_t(older version) or uint64_t(current version).

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.

      \return the data type of bitmap offset of a variable from file.
    */
    virtual FQ::DataType getBitmapOffsetType(const std::string& variable);

    //*********************
    // API for bitmap
    //*********************/

    /*!
      \brief Create the dataset to store bitmaps of a variable.
      The dataset must be created with the type uint32_t.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
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
      \brief Get the length(i.e. uint32_t) of bitmaps of a variable from file.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param nElements IN: The number of bitmaps.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode.

      \return True if success, otherwise return False.
    */
    virtual bool getBitmapLength(const std::string &variable,
                                 uint64_t *len,
                                 const uint64_t mpi_idx);

    /*!
      \brief Read a subsection of the bitmap indexes.  

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param startoffset IN: The start offset of the subsection.
      \param endoffset IN: The end offset of the subsection.
      The startoffset is always *inclusive* and the endoffset is always *exclusive*.
      The numbering of elements for the offsets starts from zero (as expected).
      \param data OUT: The values of bitmap indexes.
      The bitmaps are stored in the type of uint32_t.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode.

      \return True if success, otherwise return False.
    */
    virtual bool readBitmap(const std::string& variable,
                            const uint64_t startoffset,
                            const uint64_t endoffset,
                            uint32_t *data,
                            const uint64_t mpi_idx);

    /*!
      \brief Write a subsection of the bitmap indexes.  

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param startoffset IN: The start offset of the subsection.
      \param endoffset IN: The end offset of the subsection.
      The startoffset is always *inclusive* and the endoffset
      is always *exclusive* The numbering of elements for the offsets
      starts from zero (as expected).
      \param data IN: The values of bitmap indexes.
      \param mpi_idx IN: The MPI index in the current iteration during the index building process.
      It is only used under parallel MPI mode.

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
      \brief Set the number of bitmap keys of a variable to the offset table.
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param nkeys IN: The number of bitmap keys.
      \param mpi_iter IN: The MPI iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapKeyLength(const std::string &variable,
                                    const uint64_t nkeys,
                                    const uint64_t mpi_iter);

    /*!
      \brief Set the number of bitmap offsets of a variable to the offset table.
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param noffsets IN: The number of bitmap offsets.
      \param mpi_iter IN: The MPI iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapOffsetLength(const std::string &variable,
                                       const uint64_t noffsets,
                                       const uint64_t mpi_iter);

    /*!
      \brief Set the number of bitmaps of a variable to the offset table.
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param len IN: The number of bitmaps.
      \param mpi_iter IN: The MPI iteration during the index building process.
      It is only used under parallel MPI mode

      \return True if success, otherwise return False.
    */
    virtual bool setBitmapLength(const std::string &variable,
                                 const uint64_t len,
                                 const uint64_t mpi_iter);

    /*!
      \brief Create the offset table
      This function is only needed if offset table is used.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param nElements OUT: The number of MPI iterations during the index building process.

      \return True if success, otherwise return False.
    */
    virtual bool createOffsetTable(const std::string &variable,
                                   const uint64_t nElements);
#endif

    //*********************
    // API for data value range
    //*********************/
    /*!
      \brief Get the max/min value of a variable.
      we always know ActualRange are always two elements
      with the datatype the same as the indexed dataset.  So we do not
      need to query the length, nor do we need to supply the length when
      we read/write the Range data.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param range IN: The max and min values of a varaible.

      \return the max/min value.
    */
    virtual bool getActualRange(const std::string &variable,
                                void *range);

    /*!
      \brief Get the max/min value of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param range IN: The max and min values of a varaible.
      \param fqType IN:  Enumerated type of variable.

      \return the max/min value.
    */
    virtual bool setActualRange(const std::string& variable,
                                const void *range,
                                const FQ::DataType fqType);

    /*!
      \brief Get the max/min value of a variable.
      we always know ExpectedRange are always two elements
      with the datatype the same as the indexed dataset.  So we do not
      need to query the length, nor do we need to supply the length when
      we read/write the Range data.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param range IN: The max and min values of a varaible.

      \return the max/min value.
    */
    virtual bool getExpectedRange(const std::string &variable,
                                  void *range);

    /*!
      \brief Get the max/min value of a variable.

      \param variable IN: The location (i.e. full-path) of a variable in the file structure.
      \param range IN: The max and min values of a varaible.
      \param fqType IN:  Enumerated type of variable.

      \return the max/min value.
    */
    virtual bool setExpectedRange(const std::string& variable,
                                  const void *range,
                                  const FQ::DataType fqType);

    virtual std::string getSortedFieldName();

private:
#ifndef FQ_NOMPI
    int mpi_rank;
    int mpi_size;
    MPI_Comm mpi_comm;
#endif
    hid_t _fileId;

    std::string _curDatasetName;
    hid_t _curDatasetId;

    std::string _indexPath;


    bool __openFile(bool readOnly);
    bool __closeFile();

    hid_t __getDatasetId(const std::string &datasetName);

    // return the HDF5 object type at path location
    H5O_type_t __getHDF5ObjType(const std::string &path);

    // return -1 if the type is not supported by HDF5
    hid_t __getHDF5DataType(const FQ::DataType fqType);

    // return the FQ data type
    FQ::DataType __getFQDataType(const hid_t type);

    bool __getAllVariables(const std::string &path,
			   std::vector<std::string> &variables);
    Node* __getVariableTree(const std::string &path);

    bool __getDatasetType(const std::string &datasetName, FQ::DataType *fqType);
    bool __getDatasetSize(const std::string &datasetName, uint64_t *size);
    bool __getDatasetDimension(const std::string &datasetName,
			       std::vector <uint64_t> &dims);

    bool __getAttribute(const std::string &path, const std::string &attrName,
			void *values);
    bool __setAttribute(const std::string &path, const std::string &attrName, 
			const void *values, const hsize_t len,
			const FQ::DataType fqType);
    bool __getAttributeInfo(const std::string &variable,
	 	            const std::string &attrName, 
                            uint64_t *size,
                            FQ::DataType *type);

    bool __readData(const std::string &datasetName, void *values);
    bool __writeData(const std::string &datasetName, const void *values);

    bool __readArrayData(const std::string& datasetName,
                         const std::vector<uint64_t> &offsets,
                         const std::vector<uint64_t> &counts,
                         const std::vector<uint64_t> &strides,
                         void *data);
    bool __writeArrayData(const std::string& datasetName,
                          const std::vector<uint64_t> &offsets,
                          const std::vector<uint64_t> &counts,
                          const std::vector<uint64_t> &strides,
                          const void *data,
			  const bool collective=false);

    bool __readPointData(const std::string &datasetName,
			 const std::vector<uint64_t> &coords,
			 void *data);

    /*    bool __readBlockData(const std::string &variable,
	  const std::vector<uint64_t> &coords,
	  void *data);*/

    bool __createDataset(const std::string& datasetName, 
			 const std::vector<uint64_t > dims, 
			 const FQ::DataType fqType,
			 const bool UNLIMITED=false);

    bool __removeDataset(const std::string &variable);

    bool __extendDataset(const std::string& datasetName, 
			 const std::vector<uint64_t> &curDims,
			 const uint64_t totLen);

#ifndef FQ_NOMPI
    bool __getOffsets(const std::string &variable, 
		      const int column, 
		      uint64_t *start, 
		      uint64_t *end, 
		      const uint64_t mpi_idx);

    bool __getOffsetLength(const std::string &variable, 
			   const int column, 
			   uint64_t* len, 
			   const uint64_t mpi_idx);
#endif
};
#endif
