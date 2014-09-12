// Author: Jerry Chou
//      Lawrence Berkeley National Laboratory
#ifndef FQ_FQPARSER_H
#define FQ_FQPARSER_H
#include "qToken.h"
#include "ibis.h"

namespace FQ {
    class fqParser;
    class queryLexer;
    class queryParser;
}

/// A representation of the fastquery query.  It parses a string into an
/// qVar object.  
///
/// A query is a set of range conditions joined together with
/// logical operators.  The supported logical operators are
/// @code
/// &&, ||.
/// @endcode
///
/// The supported range conditions are equality conditions, 
/// and one-sided range conditions
///
/// - An equality condition is defined by the equal operator and its two
///   operands can be arithematic expressions, column names, numbers or
///   string literals.  On string valued columns, FastBit currently only
///   supports equality comparisons.  In such a case, the comparison is of
///   the form "column_name = column_value".  Internally, when FastBit
///   detect that the type of the column named "column_name" is
///   ibis::CATEGORY or ibis::TEXT, it will interpret the other side as
///   literal string value to be compared.  Note that if the left operand
///   of the equality operator is not a known column name, the evaluation
///   function will examine the right operand to see if it is a column
///   name.  If the right operand is the name of string-valued column, the
///   left operand will be used as string literal.
///
/// - A one-side range condtion can be defined with any of the following
///   operators, <, <=, >, and >=.  The two operands of the operator can be
///   any arithmetic expressions, column names or numbers.
///
/// An arithematic expression may contain operators +, -, *, / and ^
///
class FQ::fqParser {
public:
    /// Construct a query token tree from a query string.
    fqParser (const std::string &query);
    /// Destructor.
    ~fqParser();

    /// Return a pointer to the string form of the where clause.
    const char* getQueryStr(void) const {
	if (queryStr.empty())
	    return 0;
	else
	    return queryStr.c_str();
    }
    /// Return a pointer to the head of a list of qVar objects
    QToken* getTokens(void) {return qToken;}

protected:
    std::string clause_;	///< String version of the where clause.
    std::string queryStr;	///< String version of the input query
    QToken *qToken;		///< The expression tree.

private:
    FQ::queryLexer *lexer;	// hold a pointer for the parser

    friend class FQ::queryParser;
}; // class fqParser
#endif
