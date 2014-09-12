/* $Id: queryLexer.ll,v 1.1 2011-04-28 22:57:48 jchou Exp $ -*- mode: c++ -*-

   Author: Jerry Chou
   Lawrence Berkeley National Laboratory
*/

%{ /* C++ declarations */
/** \file Defines the tokenlizer using Flex C++ template. */

#include "queryLexer.h"		// definition of YY_DECL
#include "queryParser.hh"	// class queryParser

typedef FQ::queryParser::token token;
typedef FQ::queryParser::token_type token_type;

#define yyterminate() return token::END
%}

/* Flex declarations and options */
/*%option noyywrap*/
%option c++
%option stack
%option nounistd
%option never-interactive
%option prefix="_queryLexer_"

VARNAME	[a-zA-Z][a-zA-Z0-9_\/]*
/*VARNAME	[\/]?[a-zA-Z][a-zA-Z0-9_\/]*  */
/*VARNAME	[a-zA-Z][a-zA-Z0-9_]*   */
DIMIDX	[0-9]*(:[0-9]*)?(:[0-9]*)?
VARIDX  {DIMIDX}(,{DIMIDX})*

%%	
	/* section defining the tokens */
({VARNAME})(\[{VARIDX}\])?				{
	yylval->stringVal = new std::string(yytext, yyleng);
	return token::NAMESTR;}
[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?		{
	yylval->stringVal = new std::string(yytext, yyleng);
	return token::NUMBER;}
"||"|"&&"					{
	yylval->stringVal = new std::string(yytext, yyleng);
	return token::LOGICOP;}
"!="|"=="|"<"|">"|"<="|">="			{
	yylval->stringVal = new std::string(yytext, yyleng);
	return token::RANGEOP;}
[-+*^]						{
	yylval->stringVal = new std::string(yytext, yyleng);
	return token::ARITHMETICOP;}
[ \t\v\n]					/*ignore whitespace*/;
.						{
	return static_cast<token_type>(*yytext);}

%%
/* additional c++ code to complete the definition of class queryLexer */
FQ::queryLexer::queryLexer(std::istream* in, std::ostream* out)
    : ::_fqLexer(in, out) {
#if defined(DEBUG) && DEBUG + 0 > 1
    yy_flex_debug = true;
#endif
}

FQ::queryLexer::~queryLexer() {
}

/* function needed by the super-class of ibis::whereLexer */
#ifdef yylex
#undef yylex
#endif

int ::_fqLexer::yylex() {
    return 0;
} // ::_fqLexer::yylex

int ::_fqLexer::yywrap() {
    return 1;
} // ::_fqLexer::yywrap
