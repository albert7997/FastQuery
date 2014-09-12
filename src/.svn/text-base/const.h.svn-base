#ifndef _CONST_H
#define _CONST_H

#include "fastquery-config.h"
#include <stdint.h>

// define it to enable packing bitmaps into fixed chunk before writing them
// to file
#define FQ_PACK_BITMAP

// the chunk size of packed bitmaps
#ifndef FQ_CHUNK_SIZE
#define FQ_CHUNK_SIZE 1000000
#endif

// #define FQ_PERFORMANCE_TEST  // insert I/O blocking to synchronize
//                              // processors for performance measurement
// #define FQ_ALWAYS_INDEX      // create temporary indexing in memory
//                              // for query processing

// default MPI subarray size in parallel mode
#define FQ_DEFAULT_MPI_LEN 1000000
// default dimension rank to split data into subarrays in parallel mode
#define FQ_DEFAULT_MPI_DIM 0
// additional field name for reporting performance measurement
#define FQ_REPORT_STATISTIC "FQ_REPORT_STATISTIC"

namespace FQ {
    /// Supported data types.
    enum DataType {
        /// Unknown type, a place holder.  Can not process data of this type!
        FQT_UNKNOWN,
        /// Four-byte IEEE floating-point numbers, internally float.
        FQT_FLOAT,
        /// Eight-byte IEEE floating-point numbers, internally double.
        FQT_DOUBLE,
        /// One-byte signed integers.  Internally signed char.
        FQT_BYTE,
        /// One-byte signed integers.  Internally unsigned char.
        FQT_UBYTE,
        /// Two-byte signed integers.  Internally int16_t.
        FQT_SHORT,
        /// Two-byte unsigned integers. Internally uint16_t.
        FQT_USHORT,
        /// Four-byte signed integers.  Internally int32_t.
        FQT_INT,
        /// Four-byte unsigned integers. Internally uint32_t.
        FQT_UINT,
        /// Eight-byte signed integers.  Internally int64_t.
        FQT_LONG,
        /// Eight-byte unsigned integers.  Internally uint64_t.
        FQT_ULONG,
    };

    /// Format of the underlying files.  Current supported file formats
    /// are: HDF5, H5Part, NetCDF, pnetCDF and ADIOS BP.
    enum FileFormat {FQ_HDF5, FQ_H5Part, FQ_NetCDF, FQ_pnetCDF, FQ_BP};
    /// How the selection are represented.
    enum SelectForm {POINTS_SELECTION, LINES_SELECTION, BOXES_SELECTION};
};
#endif
