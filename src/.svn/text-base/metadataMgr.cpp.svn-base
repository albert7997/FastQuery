#include "metadataMgr.h"
#include <stack>


/*******************************************
 * Constructor and Destructor
 ********************************************/
MetadataMgr::MetadataMgr(const ArrayIODriver &dataFile) :
    fileMgr(const_cast<ArrayIODriver&>(dataFile)) {
    populateAllVariables();
} // MetadataMgr::MetadataMgr

MetadataMgr::~MetadataMgr() {
    varMap.clear();
} // MetadataMgr::~MetadataMgr

/*******************************************
 * Public APIs
 ********************************************/

/*!
  \breif Add a information entry in the repository for a variable.
  The information of the varaible is read from file.

  \param variablePath IN: the full-path of a variable in the file.

  \return false if the information of the variable cannot be retrieved
  from file.
  Also return true if the variable already exists.
*/
bool MetadataMgr::addVariable(const std::string &variable) {
    VarMap::iterator it = varMap.find(variable);
    if (it != varMap.end()) {
        LOGGER(ibis::gVerbose > 5)
            << "MetadataMgr::addVariable(" << variable.c_str()
            << "): variable already exits";
        return true; // the variable has been added.
    }

    std::vector<uint64_t> dims;
    FQ::DataType type;
    if (! fileMgr.getVariableInfo(variable, dims, &type)) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "MetadataMgr::addVariable("
            << variable.c_str() << "): invalid variable";
        return false;
    }
    varMap.insert(std::pair<std::string, VarInfo>
                  (variable, VarInfo(variable, dims, type)));
    return true;
} // MetadataMgr::addVariable

/*!
  \breif Reset and populate information for all variables in the file.
*/
unsigned int MetadataMgr::populateAllVariables() {
    // clear the current variable information
    varMap.clear();
    // get all variables
    std::vector<std::string> variables;
    fileMgr.getAllVariables("/", variables);
    for(unsigned int i=0; i<variables.size(); i++) {
        bool berr = addVariable(variables[i]);
        if (! berr) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "MetadataMgr::populateAllVariables:"
                << " failed to add variable "
                << variables[i].c_str();
        }
    }
    return varMap.size();
} // MetadataMgr::populateAllVariables

/*!
  \brief Get the full-path of all variables having a specific prefix
  and/or name.
  Variable index is not considered by this function.
  So variable index should not be included in varNameStr.

  \param varInfoList OUT: A vector containing the basic information of
  variables that match to the given prefix and name.
  \param varSpaceList OUT: A vector containing the subarray space
  information of variables that match to the given prefix and name.
  \param varPathStr IN: A prefix-path in the file structure. (optional)
  \param varNameStr IN: A variable name. (optional)

  \return the number of variables found.
*/
unsigned int MetadataMgr::getAllVariables
(std::vector<VarInfo> &varInfoList, std::vector<VarSpace> &varSpaceList,
 const std::string &varPathStr, const std::string &varNameStr) {
    std::string varStr = "";
    std::string sectionStr = "";
    size_t pos = varNameStr.find_last_of('[');
    if (pos == varNameStr.npos) {
        varStr = varNameStr;
    } else {
        varStr = varNameStr.substr(0, pos);
        sectionStr = varNameStr.substr(pos);
    }
    if (varStr[0]=='/') varStr = varStr.substr(1);

    varInfoList.clear();
    varSpaceList.clear();
    for (VarMap::const_iterator it=varMap.begin(); it!=varMap.end(); ++it) {
        const std::string &variable = (*it).first;
        if (varPathStr.compare("") == 0 ||
            _findPathStr(variable, varPathStr) == true ) {
            if (varStr.compare("") == 0 ||
                _findNameStr(variable, varStr) == true) {
                LOGGER(ibis::gVerbose > 2)
                    << "MetadataMgr::getAllVariables("
                    << varPathStr.c_str() << ", " << varNameStr.c_str() << "):"
                    << " found variable \"" << variable << "\"";

                const VarInfo &viref = (*it).second;
                std::vector<uint64_t> offsets;
                std::vector<uint64_t> counts;
                std::vector<uint64_t> strides;
                bool ret = _parseSpaceStr
                    (sectionStr, viref.getDims(), offsets,
                     counts, strides);
                if (!ret) {
                    LOGGER(ibis::gVerbose > 2)
                        << "MetadataMgr::getAllVariables("
                        << varPathStr << ", " << varNameStr
                        << ") skip invalid variable section description \""
                        << sectionStr << "\" for variable \""
                        << viref.getPath() << "\"";
                } else {
                    varInfoList.push_back(VarInfo(viref));
                    varSpaceList.push_back(VarSpace(offsets, counts, strides));
                }
            }
        }
    }
    return varInfoList.size();
} // MetadataMgr::getAllVariables

/// Find one variable that matches the specified path and name.  If more
/// than match is found, return the first one with the least number of
/// characters in the fully qualified variable name.
int MetadataMgr::getShortMatch(VarInfo &vi, VarSpace &vs,
                               const char *path, const char *name) const {
    std::string varStr;
    std::string sectionStr;
    const char *ptr = strrchr(name, '[');
    if (ptr == 0) {
        varStr = name;
    } else {
        varStr.append(name, ptr-name);
        sectionStr = ptr;
    }
    if (varStr[0]=='/') varStr = varStr.substr(1);

    VarMap::const_iterator mit = varMap.end();
    unsigned mlen = INT_MAX;
    for (VarMap::const_iterator it = varMap.begin();
         it != varMap.end(); ++ it) {
        const std::string &var = (*it).first;
        if (_findPathStr(var, path) && _findNameStr(var, varStr)) {
            if (mit == varMap.end() || mlen > var.size()) {
                mit = it;
                mlen = var.size();
            }
        }
    }
    if (mit == varMap.end())
        return 0;

    vi.copy(mit->second);
    std::vector<uint64_t> offsets, counts, strides;
    bool ps = _parseSpaceStr(sectionStr, vi.getDims(),
                             offsets, counts, strides);
    vs.setOffsets(offsets);
    vs.setCounts(counts);
    vs.setStrides(strides);
    if (ps) return 1;
    else return -1;
} // MetadataMgr::getShortMatch

/*!
  \brief Parse the query by replacing the variable names to the ones
  used in FastBit.

  \param query IN: SQL-like query string.  Variables used in the query
  should include their full paths in the file structure.  Otherwise the
  default path "/" is given to the variable.  E.g. "a1 < 5 && b2 > 10"
  is equivalent to "/a1 < 5 && /b2 > 10".

  \param varPathStr IN: A prefix-path in the file structure for all
  variables in the query.

  \param queryFB OUT: A parsed query with the variable names used by
  FastBit.

  \param varInfoList OUT: A vector containing the basic information of
  variables in the query.

  \param varSpaceList OUT: A vector containing the subarray space
  information of variables in the query.

  \return ture if success, otherwise return false.
*/
bool MetadataMgr::parseQuery(const std::string &query,
                             const std::string &varPathStr,
                             std::string &queryFB,
                             std::vector<VarInfo> &varInfoList,
                             std::vector<VarSpace> &varSpaceList,
                             int mpi_dim) {
    varInfoList.clear();
    varSpaceList.clear();
    // parse query
    LOGGER(ibis::gVerbose > 4)
        << "MetadataMgr::partQuery starting to process \"" << query << '"';
    FQ::fqParser parser(query);
    QToken* root = parser.getTokens();
    if (root == 0) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning --- MetadataMgr::parseQuery: failed to parse \""
            << query.c_str() << "\"";
        return false;
    }
    // exam query tokens and find variables
    queryFB = "";
    std::stack<QToken*> s;
    QToken* token = root;
    while(1) {
        if (token != NULL) {
            s.push(token);
            token = token->getLeft();
        } else {
            if (s.empty()) {
                break;
            } else {
                token = s.top();
                s.pop();
                if (token->getType() != QToken::VAR) {
                    queryFB += token->getText();
                } else {
                    std::string varNameStr = token->getText();
                    // get variable information
                    std::vector<VarInfo> newVarInfoList;
                    std::vector<VarSpace> newVarSpaceList;
                    unsigned int numVar = getAllVariables
                        (newVarInfoList, newVarSpaceList, varPathStr,
                         varNameStr);
                    if (numVar == 0) {
                        LOGGER(ibis::gVerbose > 0)
                            << "Warning --- MetadataMgr::parseQuery("
                            << query.c_str() << ", " << varPathStr.c_str()
                            << ") found no variable matching \""
                            << varNameStr.c_str() << "\"";
                        return false;
                    }
                    int idx = 0;
                    if (numVar > 1) {
                        // pick the best fit variable with the shortest length
                        std::string varName = newVarInfoList[0].getPath();
                        int len = varName.size();
                        for (unsigned int i=1; i<numVar; i++) {
                            if (len > newVarInfoList[i].getPath().size()) {
                                idx = i;
                                len = newVarInfoList[i].getPath().size();
                            }
                        }
                        LOGGER(ibis::gVerbose > 1)
                            << "Warning --- MetadataMgr::parseQuery("
                            << query.c_str() << ", " << varPathStr.c_str()
                            << ") found multiple matching variables, use \""
                            << newVarInfoList[idx].getPath() << "\"";
                    }
                    // check if the variable has been appeared
                    int varIdx = -1;
                    for (unsigned int i=0; i<varInfoList.size(); i++) {
                        if (varInfoList[i].getPath().compare
                            (newVarInfoList[idx].getPath()) == 0) {
                            if (varSpaceList[i].isEqual(newVarSpaceList[idx])) {
                                varIdx = (int)i;
                            }
                        }
                    }
                    if (varIdx == -1) {
                        varIdx = varInfoList.size();

                        varInfoList.push_back(newVarInfoList[idx]);
                        varSpaceList.push_back(newVarSpaceList[idx]);
                    }
                    std::string varFBName = "Var";
                    char str[25];
                    sprintf(str, "%u",varIdx);
                    varFBName += str;
                    queryFB += varFBName;
                }
                token = token->getRight();
            }
        }
    }
    if (varSpaceList.empty()) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning --- MetadataMgr::parseQuery("
            << query.c_str() << ", " << varPathStr.c_str()
            << ") didnot find any valid variable from the query";
        return false;
    }
    // check if all variables have the same total number of elements
    uint64_t numElements = varSpaceList[0].getSize();
    for (unsigned int i=1; i<varSpaceList.size(); i++) {
        uint64_t num = varSpaceList[i].getSize();
        if (num != numElements) {
            LOGGER(ibis::gVerbose > 0)
                << "Error --- MetadataMgr::parseQuery("
                << query.c_str() << ", " << varPathStr.c_str() << "):"
                << " variable \"" << varInfoList[i].getPath()
                << "\" has different number of elements from"
                << " variable \"" << varInfoList[0].getPath() << "\"";
            return false;
        }
    }
#ifndef FQ_NOMPI
    // all variables have must have the same dimension length along the
    // splitting dimension
    std::vector<uint64_t> counts = varSpaceList[0].getCounts();
    if (counts.size() < mpi_dim) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning --- MetadataMgr::parseQuery("
            << query.c_str() << ", " << varPathStr.c_str() << "):"
            << " variable \"" << varInfoList[0].getPath()
            << "\" cannot be splitted along the " << mpi_dim << "th dimension";
        return false;
    }
    for (unsigned int i=1; i<varSpaceList.size(); i++) {
        std::vector<uint64_t> tmpCounts = varSpaceList[i].getCounts();
        if (tmpCounts.size() < mpi_dim) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning --- MetadataMgr::parseQuery("
                << query.c_str() << ", " << varPathStr.c_str() << "):"
                << " variable \"" << varInfoList[i].getPath()
                << "\" cannot be splitted along the " << mpi_dim
                << "th dimension";
            return false;
        }
        if (counts[mpi_dim] != tmpCounts[mpi_dim]) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning --- MetadataMgr::parseQuery("
                << query.c_str() << ", " << varPathStr.c_str() << "):"
                << " variable \"" << varInfoList[i].getPath()
                << "\" has different dimension length from"
                << " variable \"" << varInfoList[0].getPath() << "\"";
            return false;
        }
    }
#endif
    LOGGER(ibis::gVerbose > 2)
        << "MetadataMgr::parseQuery: successfully transformed \""
        << query << "\" into \"" << queryFB << '"';
    return true;
}

/*******************************************
 * Private Functions
 ********************************************/
bool MetadataMgr::_findPathStr(const std::string &variable,
                               const std::string &varPathStr) const {
    if (variable.empty()) return false;
    if (varPathStr.empty()) return true;

    size_t pos = 0;
    size_t len = varPathStr.length();
    if (varPathStr.compare("/")==0) return true;

    while (1) {
        pos = variable.find(varPathStr, pos);
        if (pos != variable.npos) {
            if (pos == 0 || variable[pos-1] == '/') {
                if ((pos+len) == variable.length() ||
                    variable[pos+len] == '/') {
                    return true;
                }
            }
            pos++;
        } else {
            break;
        }
    }
    return false;
} // MetadataMgr::_findPathStr

bool MetadataMgr::_findPathStr(const std::string &variable,
                               const char *vp) const {
    if (variable.empty()) return false;
    if (vp == 0 || *vp == 0) return true;

    size_t pos = 0;
    size_t len = strlen(vp);
    if (*vp == '/' && vp[1] == 0) return true;

    while (1) {
        pos = variable.find(vp, pos);
        if (pos != variable.npos) {
            if (pos == 0 || variable[pos-1] == '/') {
                if ((pos+len) == variable.length() ||
                    variable[pos+len] == '/') {
                    return true;
                }
            }
            pos++;
        } else {
            break;
        }
    }
    return false;
} // MetadataMgr::_findPathStr

bool MetadataMgr::_findNameStr(const std::string &variable,
                               const std::string &varNameStr) const {
    if (variable.empty()) return false;
    if (varNameStr.empty()) return true;

    // find the last delimiter
    size_t pos = variable.find_last_of('/');
    if (pos < variable.size())
        ++ pos;
    else
        pos = 0;
    const size_t len = varNameStr.length();
    if (pos+len != variable.size())
        return false;
    else
        return (0 == variable.compare(pos, len, varNameStr));
} // MetadataMgr::_findNameStr

bool MetadataMgr::_findNameStr(const std::string &variable,
                               const char *vn) const
{
    if (variable.empty()) return false;
    if (vn == 0 || *vn == 0) return true;

    // find the last delimiter
    size_t pos = variable.find_last_of('/');
    if (pos < variable.size())
        ++ pos;
    else
        pos = 0;
    const size_t len = strlen(vn);
    if (pos+len != variable.size())
        return false;
    else
        return (0 == variable.compare(pos, len, vn));
} // MetadataMgr::_findNameStr

bool MetadataMgr::_parseSpaceStr(const std::string &desc,
                                 const std::vector<uint64_t> &dims,
                                 std::vector<uint64_t> &offsets,
                                 std::vector<uint64_t> &counts,
                                 std::vector<uint64_t> &strides) const {
    offsets.resize(dims.size());
    counts.resize(dims.size());
    strides.resize(dims.size());

    if (desc.compare("") == 0) {
        for (unsigned int i=0; i<dims.size(); i++) {
            offsets[i] = 0;
            counts[i] = dims[i];
            strides[i] = 1;
        }
    } else {
        if (desc[0] != '[' || desc[desc.size()-1] != ']') {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "MetadataMgr::_parseSpaceStr(" << desc
                << ") subarray specification must start with '[' and end "
                "with ']'";
            return false;
        }
        // get rid of the [ ]
        std::string str = desc.substr(1, desc.length()-2);
        size_t prePos = 0;
        size_t curPos;
        // parse range for each dimension
        for (unsigned int i=0; i<dims.size(); i++) {
            if (prePos == str.npos) {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- " << "MetadataMgr::_parseSpaceStr(" << desc
                    << ") does not match to the number of data dimensions: "
                    << dims.size();
                return false;
            }
            curPos = str.find(',', prePos);
            std::string dimDesc = "";
            if (curPos != str.npos) {
                dimDesc = str.substr(prePos, curPos-prePos);
                prePos = curPos+1;
            } else {
                dimDesc = str.substr(prePos);
                prePos = curPos;
            }
            if (! _parseSpaceRange(dimDesc, dims[i], &(offsets[i]),
                                   &(counts[i]), &(strides[i]))) {
                return false;
            }
        }
        if (prePos != str.npos) {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "MetadataMgr::_parseSpaceStr("
                << desc.c_str() << "):"
                << " exceeds the number of data dimensions: " << dims.size();
            return false;
        }
    }
    return true;
}

bool MetadataMgr::_parseSpaceRange
(const std::string &desc, const uint64_t maxLen, uint64_t *offset,
 uint64_t *count, uint64_t *stride) const {
    if (desc.compare("") == 0 || desc.compare(":") == 0) {
        *offset = 0;
        *count = maxLen;
        *stride = 1;
        return true;
    }

    // get start index
    size_t prePos = 0;
    size_t curPos = desc.find(':');
    if (curPos == desc.npos) {
        curPos = desc.length();
    }
    uint64_t startIdx = 0;
    for(unsigned int i=prePos; i<curPos; i++) {
        if (desc[i] > '9' || desc[i] < '0') {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "MetadataMgr::_parseSpaceRange("
                << desc.c_str() << "):"
                << " invalid start index with non-integer number";
            return false;
        } else {
            startIdx = startIdx*10+desc[i]-'0';
        }
    }
    if (startIdx >= maxLen) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "MetadataMgr::_parseSpaceRange("
            << desc.c_str() << "):"
            << " start indx out of the range: " << maxLen;
        return false;
    }
    *offset = startIdx;
    if (curPos == desc.length()) {
        *count = 1;
        *stride = 1;
        return true;
    }
    prePos = curPos+1;

    // get end index
    uint64_t endIdx = 0;
    curPos = desc.find(':', prePos);
    if (curPos == desc.npos) {
        curPos = desc.length();
    }
    if (curPos == prePos || prePos == desc.length()) {
        endIdx = maxLen - 1; // endIdx is not specified
    } else {
        for(unsigned int i=prePos; i<curPos; i++) {
            if (desc[i] > '9' || desc[i] < '0') {
                LOGGER(ibis::gVerbose > 0)
                    << "Warning -- " << "MetadataMgr::_parseSpaceRange("
                    << desc.c_str() << "):"
                    << " invalid end index with non-integer number";
                return false;
            } else {
                endIdx = endIdx*10+desc[i]-'0';
            }
        }
        endIdx--; // to exclude upper bound value from selection
    }
    if (endIdx >= maxLen) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "MetadataMgr::_parseSpaceRange("
            << desc.c_str() << "):"
            << " end index out of the range: " << maxLen;
        return false;
    }
    if (endIdx < startIdx) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "MetadataMgr::_parseSpaceRange("
            << desc.c_str() << "):"
            << " end index is less than start index";
        return false;
    }
    if (curPos == desc.length()) {
        *count = endIdx - startIdx + 1;
        *stride = 1;
        return true;
    }
    prePos = curPos+1;

    // get stride
    if (prePos == desc.length()) {
        *stride = 1; // stride is not specified
        *count = endIdx - startIdx + 1;
        return true;
    }
    *stride = 0;
    for(unsigned int i=prePos; i<desc.length(); i++) {
        if (desc[i] > '9' || desc[i] < '0') {
            LOGGER(ibis::gVerbose > 0)
                << "Warning -- " << "MetadataMgr::_parseSpaceRange("
                << desc.c_str() << "):"
                << " invalid stride index with non-integer number";
            return false;
        } else {
            *stride = (*stride) * 10 + desc[i] - '0';
        }
    }
    if (*stride == 0 ) {
        LOGGER(ibis::gVerbose > 0)
            << "Warning -- " << "MetadataMgr::_parseSpaceRange("
            << desc.c_str() << "):"
            << " stride must be a positive number";
        return false;
    }
    *count = (uint64_t)((endIdx-startIdx)/(*stride)) + 1;
    return true;
}
