/* $Id: queryParser.yy,v 1.1 2011-04-28 22:57:48 jchou Exp $ -*- mode: c++ -*- */
// Author: Jerry Chou
//      Lawrence Berkeley National Laboratory

%code top {
/** \file Defines the parser for the query accepted by FastQuery.
    The definitions are processed through bison.
*/

#include <iostream>
}
%code requires {
#include "fqParser.h"
}

/* bison declarations */
%require "2.3"
%debug
%error-verbose
%start START
%defines
%skeleton "lalr1.cc"
%name-prefix="FQ"
%define "parser_class_name" "queryParser"
%locations

%initial-action
{ // initialize location object
    @$.begin.filename = @$.end.filename = &(driver.clause_);
};

%parse-param {class fqParser& driver}

%union {
    std::string *stringVal;
    QToken *tokenNode;
}

%token 		   END	 0	"end of input"
%token <stringVal> LOGICOP 	"logical operator"
%token <stringVal> RANGEOP 	"range operator"
%token <stringVal> ARITHMETICOP "arithmetic operator"
%token <stringVal> NUMBER 	"number"
%token <stringVal> NAMESTR 	"variable name"

%right ARITHMETICOP

%type <tokenNode> queries queryExpr mathExpr

%destructor { delete $$; } NAMESTR NUMBER RANGEOP LOGICOP
%destructor { delete $$; } queries queryExpr mathExpr

%{
#include <queryLexer.h>

#undef yylex
#define yylex driver.lexer->lex
%}

%% /* Grammar rules */
queries:
  queryExpr LOGICOP queries {
#if defined(DEBUG) && DEBUG + 0 > 1
    LOGGER(ibis::gVerbose > 0)
        << __FILE__ << ":" << __LINE__ << " parsing -- " << *$1
        << " <LOGICOP> " << *$3;
#endif
  $$ = new QToken(QToken::LOGICOP, $2->c_str());
  $$->setRight($3);
  $$->setLeft($1);
}
| queryExpr {
  $$ = $1;
}
;

queryExpr:
  NAMESTR RANGEOP mathExpr {
#if defined(DEBUG) && DEBUG + 0 > 1
    LOGGER(ibis::gVerbose > 0)
        << __FILE__ << ":" << __LINE__ << " parsing -- " << *$1
        << " <RANGEOP> " << *$3;
#endif
  QToken* var = new QToken(QToken::VAR, $1->c_str());
  $$ = new QToken(QToken::RANGEOP, $2->c_str());
  $$->setRight($3);
  $$->setLeft(var);
}
;

mathExpr:
  mathExpr ARITHMETICOP mathExpr {
#if defined(DEBUG) && DEBUG + 0 > 1
    LOGGER(ibis::gVerbose > 0)
        << __FILE__ << ":" << __LINE__ << " parsing -- " << *$1
        << " <ARITHMETICOP> " << *$3;
#endif
  $$ = new QToken(QToken::ARTHIMETICOP, $2->c_str());
  $$->setRight($3);
  $$->setLeft($1);
}
| '(' mathExpr ')' {
  $$ = $2;
}
| NUMBER {
  $$ = new QToken(QToken::NUMBER, $1->c_str());
}
;

START : queries END { /* pass queries to the driver */
    driver.qToken = $1;
}
| queries ';' { /* pass qexpr to the driver */
    driver.qToken = $1;
}
;

%%
void FQ::queryParser::error(const queryParser::location_type& l,
                              const std::string& m) {
    LOGGER(ibis::gVerbose > 0)
        << "Warning -- queryParser encountered " << m
        << " at location " << l;
}
