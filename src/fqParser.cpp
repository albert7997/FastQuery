// Author: Jerry Chou
//      Lawrence Berkeley National Laboratory
#include "fqParser.h"
#include "queryLexer.h"

FQ::fqParser::~fqParser() 
{
    delete qToken;
} // fqParser::~fqParser

FQ::fqParser::fqParser(const std::string &query) : qToken(0)
{
    if (query.empty()) return;

    LOGGER(ibis::gVerbose > 2)
	<< "fqParser:: to parse query \"" << query.c_str() << "\"";

    queryStr = query;
    std::istringstream iss(queryStr);
    ibis::util::logger lg;
    FQ::queryLexer lx(&iss, &(lg()));
    FQ::queryParser parser(*this);
    lexer = &lx;
#if DEBUG+0 > 0
    parser.set_debug_level(ibis::gVerbose);
    parser.set_debug_stream(lg());
#endif
    if ( parser.parse() ) {
	LOGGER(ibis::gVerbose > 0)
	    << "Warning -- fqParser::(" << query
	    << "): failed to parse the string into a list of variable names";
	if (qToken != 0) {
	    delete qToken;
	    qToken = 0;
	}
    }
} // fqParser::fqParser
