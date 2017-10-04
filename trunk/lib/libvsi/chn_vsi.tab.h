/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_CHN_VSI_TAB_H_INCLUDED
# define YY_YY_CHN_VSI_TAB_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     RESIDUE = 258,
     MOLECULE = 259,
     FILENAME = 260,
     PDBFILE = 261,
     LETTER = 262,
     COLOR = 263,
     CHAIN = 264,
     ECHO_CMD = 265,
     COMMENT = 266,
     AN_ATOM = 267,
     HYDROGEN = 268,
     INTEGER = 269,
     VIEW = 270,
     CARBONYL = 271,
     JUNK = 272,
     NEWLINE = 273,
     END_OF_FILE = 274,
     CLEAR = 275
   };
#endif


#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 2058 of yacc.c  */
#line 83 "chn_vsi.y"

		int	dummy;		/* no (null) attributes */
		char 	*pcval;		/* pointer to character string */
		char 	*mval;		/* pointer to molecule character string */
		char 	cval;		/* character */
		Int4	*Vval;		/* view of protein */
		char 	aval[9];	/* atom type*/
		Int4	ival;		/* integer type attribute */
		double	rval;		/* double (real) type attribute */
		vsi_typ	Lval;	// item list 
		char	Cval[3];	// color [0]='Y'; [1]=dots? (boolean)
		char	Sval[200];	// color [0]='Y'; [1]=dots? (boolean)
	

/* Line 2058 of yacc.c  */
#line 92 "chn_vsi.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_CHN_VSI_TAB_H_INCLUDED  */
