// Author: Jerry Chou
//      Lawrence Berkeley National Laboratory
#ifndef FQ_QUERYLEXER_H
#define FQ_QUERYLEXER_H
/** \file
    Declares the name queryLexer.  Defines the tokenizer with two
    arguments to satisfy the reentrant parser defined in queryParser.yy.
 */
#ifndef YY_DECL
// the new lex function to satisfy the reentrant parser requirement
#define YY_DECL FQ::queryParser::token_type FQ::queryLexer::lex \
    (FQ::queryParser::semantic_type* yylval, \
     FQ::queryParser::location_type* yylloc)
#endif
#include "queryParser.hh"	// class queryParser

// rename yyFlexLexer to _fqLexer
#undef yyFlexLexer
#define yyFlexLexer _fqLexer
#include <FlexLexer.h>
//#undef yyFlexLexer

namespace FQ {
    /// Defines a new class with the desired lex function for C++ output of
    /// bison.
    ///
    /// @note This version of the lexer converts hexadecimal numbers to
    /// double precision floating-point numbers, which is not suitble for
    /// handling long integers.
    /// @note This version of the lexer does not distinguish between quoted
    /// strings and unquoted strings.  In cases where a string literal is
    /// needed such as for string matches, the evaluation engine will take
    /// one string as the column name and the other as a string literal.
    /// To ensure a single string is treated as a string literal, use the
    /// expression in the form of
    /// @code
    /// column_name IN ( string_literal )
    /// @endcode
    class queryLexer : public ::_fqLexer {
    public:
	queryLexer(std::istream* in=0, std::ostream* out=0);
	virtual ~queryLexer();

	// The new lex function.  It carries the value of token and its type.
	// The value of the token is returned as the first argument and the
	// corresponding type is the return value of this function.
	virtual FQ::queryParser::token_type
	lex(FQ::queryParser::semantic_type*, FQ::queryParser::location_type*);

	void set_debug(bool);
    }; // queryLexer
} // namespace FQ
#endif
