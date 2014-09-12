#ifndef _FQ_VAR_H
#define _FQ_VAR_H

#include <vector>
#include <string>
#include <sstream>
#include "const.h"

#ifndef FQ_NOMPI
#include <mpi.h>
#endif

class VarSpace {
public:
    VarSpace(const std::vector<uint64_t> &o,
             const std::vector<uint64_t> &c,
             const std::vector<uint64_t> &s)
	: offsets(o), counts(c), strides(s){
        setText();
#ifndef FQ_NOMPI
        mpi_idx = 0;
        mpi_iter = 0;
        mpi_max_iter = 0;
#endif
    }

    const std::vector<uint64_t>& getOffsets() const {return offsets;}
    const std::vector<uint64_t>& getCounts() const {return counts;}
    const std::vector<uint64_t>& getStrides() const {return strides;}
    void setOffsets(const std::vector<uint64_t> &in) {
        offsets.reserve(in.size());
        offsets.clear();
        offsets.insert(offsets.end(), in.begin(), in.end());}
    void setCounts(const std::vector<uint64_t> &in) {
        counts.reserve(in.size());
        counts.clear();
        counts.insert(counts.end(), in.begin(), in.end());}
    void setStrides(const std::vector<uint64_t> &in) {
        strides.reserve(in.size());
        strides.clear();
        strides.insert(strides.end(), in.begin(), in.end());}

    uint64_t getSize() const {
        if (counts.size() == 0) return 0;
        uint64_t size = counts[0];
        for(unsigned int i=1; i<counts.size(); ++i) {
            size *= counts[i];
        }
        return size;
    }

    bool isEqual (const VarSpace &var) const {
        std::vector<uint64_t> varOffsets = var.getOffsets();
        if (offsets.size() != varOffsets.size()) return false;
        for(unsigned int i=0; i<offsets.size(); i++) {
            if (offsets[i] != varOffsets[i]) return false;
        }
        std::vector<uint64_t> varCounts = var.getCounts();
        if (counts.size() != varCounts.size()) return false;
        for(unsigned int i=0; i<counts.size(); i++) {
            if (counts[i] != varCounts[i]) return false;
        }
        std::vector<uint64_t> varStrides = var.getStrides();
        if (strides.size() != varStrides.size()) return false;
        for(unsigned int i=0; i<strides.size(); i++) {
            if (strides[i] != varStrides[i]) return false;
        }
        return true;
    }

    /// Extract the text version of the hyperslab.
    std::string getText() const {return name;}
    /// Generate the text version of the hyperslab.
    void setText() {
        std::ostringstream text("");
        text << "[";
        if (! offsets.empty())
            text << offsets[0];
        text << ":";
        if (! counts.empty())
            text << counts[0];
        text << ":";
        if (! strides.empty())
            text << strides[0];
        for(unsigned int i=1; i<counts.size(); i++) {
            text << ",";
            if (i < offsets.size())
                text << offsets[i];
            text << ":";
            if (i < counts.size())
                text << counts[i];
            text << ":";
            if (i < strides.size())
                text << strides[i];
        }
        text << "]";

        name = text.str();
    }

#ifndef FQ_NOMPI
    uint64_t getMpiDim() {return mpi_dim;}
    uint64_t getMpiLen() {return mpi_len;}
    uint64_t getMpiIdx() {return mpi_idx;}
    uint64_t getMpiIter() {return mpi_iter;}
    uint64_t getMpiMaxIter() {return mpi_max_iter;}
    MPI_Comm getMpiComm() {return mpi_comm;}
    void setMpiInfo(std::string text, uint64_t dim, uint64_t len,
                    uint64_t idx, uint64_t iter, uint64_t maxIter,
                    MPI_Comm comm) {
        mpi_dim = dim;
        mpi_len = len;
        mpi_idx = idx;
        mpi_iter = iter;
        mpi_max_iter = maxIter;
        mpi_comm = comm;
        std::ostringstream tmp("");
        tmp << text << "-" << dim << "-" << len;
        name = tmp.str();
    }
#endif

private:
    std::vector<uint64_t> offsets;
    std::vector<uint64_t> counts;
    std::vector<uint64_t> strides;
    std::string name;

#ifndef FQ_NOMPI
    std::string mpi_text;
    uint64_t mpi_dim;
    uint64_t mpi_len;
    uint64_t mpi_idx;
    uint64_t mpi_iter;
    uint64_t mpi_max_iter;
    MPI_Comm mpi_comm;
#endif
};

class VarInfo {
public:
    VarInfo(const std::string filePath,
            std::vector<uint64_t> dims,
            const FQ::DataType type) {
        this->filePath = filePath;
        this->type = type;
        this->dims = dims;
    }

    VarInfo(const VarInfo& var) {
        this->filePath = var.getPath();
        this->type = var.getType();
        this->dims = var.getDims();
    }

    std::string getPath() const {return filePath;}
    FQ::DataType getType() const {return type;}
    const std::vector<uint64_t>& getDims() const {return dims;}
    unsigned int getNDims() const {return dims.size();}
    uint64_t getSize() const {
        if (dims.size() == 0) return 0;
        uint64_t size = dims[0];
        for(unsigned int i=1; i<dims.size(); i++) {
            size*= dims[i];
        }
        return size;
    }

    void copy(const VarInfo &rhs) {
        filePath = rhs.filePath;
        type = rhs.type;
        dims.clear();
        dims.reserve(rhs.dims.size());
        dims.insert(dims.end(), rhs.dims.begin(), rhs.dims.end());}

private:
    std::string filePath;
    FQ::DataType type;
    std::vector<uint64_t> dims;
};
#endif
