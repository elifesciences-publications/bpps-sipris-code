/******************************************************************************************
    Copyright (C) 1997-2014 Andrew F. Neuwald, Cold Spring Harbor Laboratory
    and the University of Maryland School of Medicine.

    Permission is hereby granted, free of charge, to any person obtaining a copy of 
    this software and associated documentation files (the "Software"), to deal in the 
    Software without restriction, including without limitation the rights to use, copy, 
    modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
    and to permit persons to whom the Software is furnished to do so, subject to the 
    following conditions:

    The above copyright notice and this permission notice shall be included in all 
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
    PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
    OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
    OTHER DEALINGS IN THE SOFTWARE.

    For further information contact:
         Andrew F. Neuwald
         Institute for Genome Sciences and
         Department of Biochemistry & Molecular Biology
         University of Maryland School of Medicine
         801 West Baltimore St.
         BioPark II, Room 617
         Baltimore, MD 21201
         Tel: 410-706-6724; Fax: 410-706-1482; E-mail: aneuwald@som.umaryland.edu
 ******************************************************************************************/

/* tree2hpt program... */
%{
#include "tree2hpt.h"
/**************************** Global Variables ******************************/
int	yyparse();
int     yylex();
UInt4	PARSELINE=1;
btn_typ *Root;

extern FILE *yyin;
extern FILE *yyout;

BooLean	VERBOSE=TRUE;

extern int yyerror(char *s);

/********************* Declarations of Action Symbols ***********************/
%}
%start start
%token <pcval>	STRING			/*  name of node */
%token <ival>	NODE			/*  integers */
%token <btn>	REG_NODE		/*  a non-root node */
%token <ival>	FIRSTNODE			/*  integers */
%token <ival>	ROOT			/*  integers */
%token <btn>	ROOT_NODE		/*  a root node */
%token <cval>	JUNK			/*  syntax error in lex program */
%token <rval>   DISTANCE                    /*  float value */
%token <ival>	LEVEL_B			/*  integers */
%token <ival>	LEVEL_M			/*  integers */
%token <ival>	LEVEL_E			/*  integers */
%token <ival>	FALSE_NODE
%token <dummy>	ROOT_ZERO
%token <dummy>	IGNORE

%token <dummy>	NEWLINE
%token <dummy>	END_OF_FILE		/* require program to come to parse entire file.*/
%union	{
		int	dummy;		/* no (null) attributes */
		char 	*pcval;		/* pointer to character string */
		char 	cval;		/* character */
		double	rval;		/* float value */
		Int4	ival;		/* integer type attribute */
		btn_typ *btn;		// binary tree node
		char	Sval[200];	// string
		int	Pval[200];	// path for file nodes to root
	}
%type	<btn>	root tree child children
%type	<Sval>	junk fatal_error 
%type	<dummy>	left_paren end_file comma ignore
// e.g. input string: (((7,((15,16):11,12,(17,18):13,14):8):3,4,(9,10):5,6):1,2):0;
// each tree at least two children required
%%
start		: tree ';' end_file {
			Root=$1;
			if(VERBOSE) fprintf(stderr,"Input accepted.\n"); YYACCEPT; 
		}
		| fatal_error end_file { yyerror($1); YYERROR; }
		;

tree		: left_paren children root  { $3->AddLeftChild($2); $$ = $3; } 
		| left_paren tree root  { $3->AddLeftChild($2); $$ = $3; } 
		| left_paren child root  { $3->AddLeftChild($2); $$ = $3; } 
		;

children	: child comma children 	{ $1->AddRightChild($3); $$=$1;  }
		| tree comma children { $1->AddRightChild($3); $$=$1;  }
		| child comma child { $1->AddRightChild($3); $$=$1;  }
		| tree comma child { $1->AddRightChild($3); $$=$1;  }
		| child comma tree { $1->AddRightChild($3); $$=$1;  }
		| tree comma tree { $1->AddRightChild($3); $$=$1;  }
		;

child		: NODE { $$ = new btn_typ($1,0); if(0) fprintf(stderr,"%d",$1); }
		| NODE DISTANCE { $$ = new btn_typ($1,0); if(0) fprintf(stderr,"%d:%.2f",$1,$2); }
		| REG_NODE { $$ = $1; }
		| REG_NODE DISTANCE { $$ = $1; }
		// | NODE STRING { $$ = new btn_typ($1,$2); if(0) fprintf(stderr,"%d_%s",$1,$2); }
		// | NODE STRING DISTANCE { $$ = new btn_typ($1,$2); if(0) fprintf(stderr,"%d_%s:%.2f",$1,$2,$3); }
		;

comma		: ',' { if(0) fprintf(stderr,","); } ;
left_paren	: '(' { if(0) fprintf(stderr,"("); } ;

root		: ROOT { $$ = new btn_typ($1,0); if(0) fprintf(stderr,")%d",$1); }
		| ROOT DISTANCE { $$ = new btn_typ($1,0); if(0) fprintf(stderr,")%d:%.2f",$1,$2); }
		| ROOT_NODE { $$ = $1; }
		| ROOT_NODE DISTANCE { $$ = $1; }
		// | ROOT STRING { $$ = new btn_typ($1,$2); if(0) fprintf(stderr,")%d_%s",$1,$2); }
		// | ROOT STRING DISTANCE { $$ = new btn_typ($1,$2); if(0) fprintf(stderr,")%d_%s:%.2f",$1,$2,$3); }
		;

fatal_error	: junk { strcpy($$,$1); } 
		;

end_file	: END_OF_FILE 
		| ignore END_OF_FILE 
		;

ignore		: ignore IGNORE  
		| IGNORE
		;

junk		: junk JUNK { 
			char s[3]; s[0]=$2; s[1]=0; strcpy($$,$1); strcat($$,s); 
			if( strlen($$) > 90) yyerror($$); 
			}
		| JUNK { $$[0]=$1; $$[1]=0; $$[2]=0; }
		;

%%

int	yyerror(char *s)
{ 
	fprintf(stderr,"\nfatal error (line %d): '%s'.\n",PARSELINE,s);
	// fprintf(stderr,"%s --> %d: %s",s,yylineno,yytext);
	print_error("Input rejected.\n");
	return 1;
}

void    Close(FILE *fptr) {if(fptr != stderr && fptr != stdout) fclose(fptr); }


btn_typ	*TreeToHpt(FILE *ifp)
// e.g. input string: (((7,((15,16):11,12,(17,18):13,14):8):3,4,(9,10):5,6):1,2):0;
{ 
	yyin = ifp;
	yyout=stderr; yyparse(); // fclose(yyin);
	return Root;
}

