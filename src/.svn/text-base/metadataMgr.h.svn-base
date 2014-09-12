#ifndef _METADATA_MGR_H
#define _METADATA_MGR_H

#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <ibis.h>
#include "const.h"
#include "arrayIODriver.h"
#include "fqParser.h"
#include "fqVar.h"

/// A data structure for storing different variables.
typedef std::map<std::string, VarInfo> VarMap;

/// A data structure to hold information about different variables.
class MetadataMgr {
public:
    MetadataMgr(const ArrayIODriver &fileMgr);

    ~MetadataMgr();

    bool addVariable(const std::string &variable);
    unsigned int populateAllVariables();

    unsigned int getAllVariables(std::vector<VarInfo> &varInfoList,
                                 std::vector<VarSpace> &varSpaceList,
                                 const std::string &varPathStr="",
                                 const std::string &varNameStr="");
    int getShortMatch(VarInfo &vi, VarSpace &vs,
                      const char * =0, const char * =0) const;

    bool parseQuery(const std::string &query,
                    const std::string &varPathStr,
                    std::string &queryFB,
                    std::vector<VarInfo> &varInfoList,
                    std::vector<VarSpace> &varSpaceList,
                    int mpi_dim = FQ_DEFAULT_MPI_LEN);

    const std::string& getFileName() const {return fileMgr.getFileName();}

private:
    bool _findPathStr(const std::string &variable,
                      const std::string &varPathStr) const;
    bool _findPathStr(const std::string &variable,
                      const char *vp) const;

    bool _findNameStr(const std::string &variable,
                      const std::string &varNameStr) const;

    bool _findNameStr(const std::string &variable,
                      const char *) const;

    bool _parseSpaceStr(const std::string &desc,
                        const std::vector<uint64_t> &dims,
                        std::vector<uint64_t> &offsets,
                        std::vector<uint64_t> &counts,
                        std::vector<uint64_t> &strides) const;
    bool _parseSpaceRange(const std::string &desc,
                          const uint64_t maxLen,
                          uint64_t *offset,
                          uint64_t *count,
                          uint64_t *stride) const;

    /// All known variables.
    std::map<std::string, VarInfo> varMap;
    /// Reference to the file.
    ArrayIODriver &fileMgr;
};
#endif
