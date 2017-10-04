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

%e 6000
%p 6000
%n 6000
%k 6000
%a 6000
%o 6000
%START	NODE_INPUT FILE_INPUT STD_INPUT
I       (0|[1-9][0-9]*)
R       ((0|[1-9][0-9]*)[.][0-9]+)
P       ([1-9][0-9]*)
L	[A-Za-z]
C	([0-9]|{L})
STR     [^ _:,;)(]+
NSTR    ({L}({C}|[_+])+)
S       [ \t]+
SP      [ \t]*
%{
// #include "stdinc.h"
#include "btn_typ.h"
#include "y.tab.h"
#define YYLMAX	10000

extern UInt4   PARSELINE;
extern FILE *yyout;
char DEBUG=0;

%}
%%
<<EOF>>	{ return END_OF_FILE; }

<FILE_INPUT>.	|
<FILE_INPUT>[ \t\n]	{
			BEGIN FILE_INPUT;
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
			return IGNORE;
		}

[:]{R}		{ 
			BEGIN STD_INPUT;
			sscanf(yytext,":%lf",&yylval.rval); return(DISTANCE); 
		}

[)]{I}_{STR}	{
			Int4 i;
			BEGIN STD_INPUT;
	  		char *s,str[100];
			sscanf(yytext,")%d_%s",&i,str); 
			yylval.btn=new btn_typ(i,str);
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
	  		return(ROOT_NODE);
		}

[)]{NSTR}	{
			BEGIN STD_INPUT;
			yytext++;
			yylval.btn=new btn_typ(0,yytext);
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
	  		return(ROOT_NODE);
		}

[)]{I}		{
			BEGIN STD_INPUT;
	  		sscanf(yytext,")%d",&yylval.ival); 
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
	  		return(ROOT);
		}

<NODE_INPUT>{I}[_]{STR}	{
			BEGIN STD_INPUT;
			Int4 i;
	  		char *s,str[100]; sscanf(yytext,"%d_%s",&i,str); 
			yylval.btn=new btn_typ(i,str);
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
	  		return(REG_NODE);
		}


<NODE_INPUT>{NSTR}	{ 		// [_]{STR}[ ]{L}	|
			BEGIN STD_INPUT;
			yylval.btn=new btn_typ(0,yytext);
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
			return REG_NODE;
		}

<NODE_INPUT>{I}		{ 
			BEGIN STD_INPUT;
	  		sscanf(yytext,"%d",&yylval.ival); 
			if(DEBUG) fprintf(stderr,"%s\n",yytext);
	  		return(NODE);
		}


	/****************** single character tokens ***********************/
","	{ BEGIN NODE_INPUT; return(',');}
";"	{ BEGIN FILE_INPUT; return(';');}
"("	{ BEGIN NODE_INPUT; return('(');} 

[ \t]	{ BEGIN STD_INPUT; }

.	{ 
		BEGIN STD_INPUT;
		fprintf(stderr,"Invalid token (line %d): '%s'\n",PARSELINE,yytext);
		yylval.cval=yytext[0];
		return JUNK;
	}

[\n]	{
		BEGIN STD_INPUT;
		PARSELINE++; // return NEWLINE; 
	} 
%%

void	 yyerror()
{
	fprintf(stderr,"%d: %s\n",PARSELINE, yytext);
}

int yywrap() { BEGIN 0; return(1); }

