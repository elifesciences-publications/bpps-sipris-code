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
%START	CIRCUMFLEX ATSIGN PLUS_SIGN DOT SQUIGGLE UNDERLINE DASH LEFT_PARAN TURNOFF CLOSED COLON EXPLICIT LEFT_BRACKET WATER ISMOL
E	([Ee][+-]?{I})
R	([+-]?{I}"."[0-9]+{E}?|[+-]?{I}{E}?)
A       [CGASTNDEQKRHWYFVILMPX]
M       [pgms]
I       (0|[1-9][0-9]*)
N       [+-]?{I}
C       [ROYGCBMPWDLSroygcbmpwdls]
CS      [AROYGCBMPWXDLTSKaroygcbmpwdltsx]
NOFILE	[^0-9(\n]
L       [A-Z]
CH      ([A-Z]|[a-z])
MOL     ([A-Z]+|\[[A-Z0-9]+\])
ATOM1   [a-z]
ATOM2   ([a-z][0-9]|[0-9][a-z]|[a-z][a-z])
ATOM3   ([a-z][a-z][0-9]|[a-z][0-9][a-z*]|[0-9][a-z][a-z*]|[0-9][a-z][0-9]|[a-z][a-z][a-z']|[a-z][0-9][0-9']|[a-z][0-9][0-9,])
ATOM4   ([a-z0-9'*][a-z][a-z0-9'][a-z0-9'*]|[a-z0-9'*][a-z0-9][a-z,][a-z0-9'*]|[a-z][a-z0-9'][a-z0-9'*][a-z0-9'*])
ATOM4B   ([0-9'][a-z][a-z][0-9']|[0-9'*][a-z][a-z][a-z'*]|[0-9'*][a-z][0-9][a-z'*]|[0-9'][a-z][0-9][0-9']|[a-z][a-z][0-9][0-9']|[a-z][a-z][0-9][a-z']|[a-z][0-9][0-9][a-z']|[a-z][0-9]['][']|[a-z][0-9][0-9][0-9'])
ATOM	({ATOM1}|{ATOM2}|{ATOM3}|{ATOM4})
S       [ \t]+
SP      [ \t]*
%{
#include "chn_vsi.h"
#include "y.tab.h"
#define YYLMAX	10000

// MOL     ([A-Z]+|\[[A-Z]+[0-9]*\]|\[[0-9]+[A-Z]+\]|\[[A-Z]+[0-9]+[A-Z]+\]|\[[0-9]+\])

BooLean USE_TRACE=TRUE;

extern UInt4   PARSELINE;
extern FILE *yyout;
extern char *FileName;
extern Int4 TRACE_WIDTH;
extern Int4 NUMBER_OF_PDB_FILES;
extern Int4 WIRE_WIDTH;
extern Int4 THIN_WIRE_WIDTH;
extern Int4 SPACEFILL;
extern Int4 BIG_SPACEFILL;
extern char KEY_CHAIN;
extern Int4 KEY_FILE;
extern BooLean	VERBOSE;
Int4 CURRENT_FILE=0;

char DEBUG=0;
BooLean Implicit=FALSE;
BooLean IsOkay=TRUE;
BooLean PrintComment=FALSE;

Int4	PrintArchaicWarning(char *msg)
{ char str[100]; strncpy(str,msg,20); str[20]=0;
	fprintf(stderr,"line %d (archaic term): \"%s...\" (ignored)\n",PARSELINE,str); 
}

%}
%%
	/****************** PARAMETER REDEFINES ***********************/
^TraceWidth[=]{I}[.]$   |
^TraceWidth[=]{I}[.]{S} {  sscanf(yytext,"TraceWidth=%d.",&TRACE_WIDTH); }

^Wirewidth[=]{I}[.]$    |
^Wirewidth[=]{I}[.]{S}  {  sscanf(yytext,"Wirewidth=%d.",&WIRE_WIDTH); }

^wirewidth[=]{I}[.]$    |
^wirewidth[=]{I}[.]{S}  { sscanf(yytext,"wirewidth=%d.",&THIN_WIRE_WIDTH); }

^Spacefill[=]{I}[.]$    |
^Spacefill[=]{I}[.]{S}  {  sscanf(yytext,"Spacefill=%d.",&BIG_SPACEFILL); }

^spacefill[=]{I}[.]$    |
^spacefill[=]{I}[.]{S}  {  sscanf(yytext,"spacefill=%d.",&SPACEFILL); }

^text[=].*$          { if(isprint(yytext[5])) fprintf(stderr,"%s\n",yytext+5); }
^stderr[=].*$         { if(isprint(yytext[7])) fprintf(stderr,"%s\n",yytext+7); }

	/****************** File parameters ***********************/
^File{I}[=][^:]+:{L}$	|
^File{I}[=][^:]+:{L}{S} { 
			  Int4 file,i=0;
			  char str[500],c;
			  NUMBER_OF_PDB_FILES++;
			  if(VERBOSE) fprintf(stderr,"%s ",yytext);
                          if(sscanf(yytext,"File%d=%[^:]:%c",&file,str,&c) != 3)
					print_error("input error");
			  PrintComment=TRUE;
			  if(file == KEY_FILE){
				if(c=='X') KEY_CHAIN=0; else KEY_CHAIN=c;
			  	Implicit=TRUE; FileName=AllocString(str);
				BEGIN 0;
				return PDBFILE; 
			  } // else if(FileName==0) BEGIN CLOSED;
			}

<CLOSED>.+$	{ BEGIN CLOSED; if(DEBUG) fprintf(stderr,"%s",yytext); }
<CLOSED>[\n]	{ BEGIN CLOSED; PARSELINE++; if(DEBUG) fprintf(stderr,"\n"); }

^close:{I}[-]{I}[.].*$           { 
			  Int4 s,e;
			  if(sscanf(yytext,"close:%d-%d.",&s,&e) != 2) 
					print_error("close input error");
			  if(KEY_FILE >= s && KEY_FILE <= e) BEGIN CLOSED;
			  // PrintArchaicWarning(yytext); 
			}

^close:{I}[.]           { 
			  Int4 n;
			  if(sscanf(yytext,"close:%d.",&n) != 1) 
					print_error("close input error");
			  if(KEY_FILE == n) BEGIN CLOSED;
			  // if(KEY_FILE == n) Implicit=FALSE;
			  // PrintArchaicWarning(yytext); 
			}

	/****************** on/off of input file regions ***********************/
^[#]on.*$    {  BEGIN 0; if(DEBUG) fprintf(stderr,"%s\n",yytext); 
		if(Implicit) CURRENT_FILE=KEY_FILE; else CURRENT_FILE=0; }

^[#]off.*$        { BEGIN TURNOFF; if(DEBUG) fprintf(stderr,"%s\n",yytext); }

<TURNOFF>.+$	{ BEGIN TURNOFF; if(DEBUG) fprintf(stderr,"%s\n",yytext); }
<TURNOFF>[\n]	{ BEGIN TURNOFF; PARSELINE++; return NEWLINE; }

	/****************** Input file identification ***********************/
	/********** These tokens are not passed on to yacc *************/
[(]{I}[-]{I}[+][)][.]	{ Int4 start,end;
			  BEGIN 0;
                          if(sscanf(yytext,"(%d-%d+).",&start,&end) != 2)
                                print_error("Implicit input error");
			  if(KEY_FILE >= start && KEY_FILE <= end) Implicit=TRUE;
                        }

[(]{I}[+][)][.]         { Int4 n;
			  BEGIN 0;
                          if(sscanf(yytext,"(%d+).",&n) != 1) print_error("should not happen");
                          if(KEY_FILE == n) Implicit=TRUE; 
                        }

[(]{I}[-]{I}[-][)][.]   { Int4 n,start,end;
			  BEGIN 0;
                          if(sscanf(yytext,"(%d-%d-).",&start,&end) != 2) print_error("What?");
                          if(KEY_FILE >= start && KEY_FILE <= end) Implicit=FALSE;
                        }

[(]{I}[-][)][.]         { Int4 n;
			  BEGIN 0;
                          if(sscanf(yytext,"(%d-).",&n) != 1) print_error("Implicit input error");
                          if(KEY_FILE == n) Implicit=FALSE; 
                        }

^{SP}[(]{I}[-]{I}[)]/{I}	{
		   	 Int4 s,e;
			 BEGIN 0;
		   	 while(isspace(yytext[0])) yytext++;
		   	 sscanf(yytext,"(%d-%d)",&s,&e);
		   	 if(KEY_FILE >= s && KEY_FILE <= e){
			   CURRENT_FILE=KEY_FILE;
		   	 } else CURRENT_FILE=0;
			}

^{SP}[(]{I}[)]/{I}	{  Int4 i;
		   BEGIN 0;
		   while(isspace(yytext[0])) yytext++;
		   sscanf(yytext,"(%d)",&i);
		   CURRENT_FILE=i;
		}

<WATER>HOH/{I}	{   BEGIN 0; 
		    // fprintf(stderr,"%s (%d; %d)\n",yytext,CURRENT_FILE,KEY_FILE);
		    if(CURRENT_FILE==KEY_FILE){
		        yylval.pcval=AllocString("HOH");
			return MOLECULE;
		    }
		}

^{SP}HOH/{I}	{   BEGIN 0; 
		    if(CURRENT_FILE==KEY_FILE){
		        yylval.pcval=AllocString("HOH"); return MOLECULE;
		    }
		}

","/HOH{I}	{ BEGIN WATER; if(CURRENT_FILE==KEY_FILE) return(','); }

^{SP}{I}/HOH{I}	{ Int4 file;  
		   BEGIN WATER;
		   while(isspace(yytext[0])) yytext++;
		   sscanf(yytext,"%d",&file); CURRENT_FILE=file;
		   //  fprintf(stderr,"found: %s (%d)\n",yytext,file);
		}

^{SP}{I}[-]{I}/[!]	|
^{SP}{I}[-]{I}/{A}{I}	{
		   Int4 s,e;
		   BEGIN 0;
		   while(isspace(yytext[0])) yytext++;
		   sscanf(yytext,"%d-%d",&s,&e);
		   if(KEY_FILE >= s && KEY_FILE <= e){
			   CURRENT_FILE=KEY_FILE;
		   } else CURRENT_FILE=0;
		}

^{SP}{I}/[!]	|
^{SP}{I}/{A}{I}	{ Int4 n,i,file;  
		   BEGIN 0;
		  // fprintf(stderr,"%s\n",yytext);
		   while(isspace(yytext[0])) yytext++;
		   sscanf(yytext,"%d",&file);
		   CURRENT_FILE=file;
		}

	/****************** Ignored strings ***********************/
^[#].*$         { BEGIN 0; if(DEBUG) fprintf(stderr,"%s\n",yytext); }

[/][/].*$               { BEGIN 0;
				if(VERBOSE && PrintComment) fprintf(stderr,"%s",yytext);
				yylval.pcval=AllocString(yytext); return COMMENT; 
			}

^define_.*$		{ BEGIN 0; PrintArchaicWarning(yytext); }
^[!]File{I}[-]{I}[=].*$ { BEGIN 0; PrintArchaicWarning(yytext); }
^[!]File{I}[=].*$       { BEGIN 0; PrintArchaicWarning(yytext); }
^Dots[=].*$             { BEGIN 0; PrintArchaicWarning(yytext); }
^Dots{I}[=].*$          { BEGIN 0; PrintArchaicWarning(yytext); }

	/****************** Command strings ***********************/
^clear[.]               { BEGIN 0;  if(CURRENT_FILE==KEY_FILE) return CLEAR; }

^clear{I}[.]            { BEGIN 0;
			  Int4 file; sscanf(yytext,"clear%d.",&file);
			  // fprintf(stderr,"%s\n",yytext);
			  if(file==KEY_FILE) return CLEAR; 
			}

^echo[=].*$             { BEGIN 0; 
			   if(CURRENT_FILE==KEY_FILE){
		    	      while(yytext[0] != '=') yytext++;
			      yytext++; yylval.pcval=AllocString(yytext);
			      return ECHO_CMD;
			   }
			}

^echo{I}[=].*$          { BEGIN 0;
			   Int4 file; sscanf(yytext,"echo%d=",&file);
			   if(file==KEY_FILE){
		    	      while(yytext[0] != '=') yytext++;
			      yytext++; yylval.pcval=AllocString(yytext);
			      return ECHO_CMD;
			   }
			}

^echo{I}[-]{I}[=].*$	{ Int4 start,end;
			  BEGIN 0;
                          if(sscanf(yytext,"echo%d-%d=).",&start,&end) != 2)
                                print_error("Implicit input error");
			  if(KEY_FILE >= start && KEY_FILE <= end){
		    	      while(yytext[0] != '=') yytext++;
			      yytext++; yylval.pcval=AllocString(yytext);
			      return ECHO_CMD;
			  }
                        }

^view{I}[:]R[=]{N}[,]{N}[,]{N}[;]T[=]{R}[,]{R}[,]{N}[.]            |
^file{I}[:]R[=]{N}[,]{N}[,]{N}[;]T[=]{R}[,]{R}[,]{N}[.]            { BEGIN 0; 
                          Int4  i,i1,i2,rx,ry,rz,tz;
			  double tx,ty;
			  Int4	*view;
			  if(yytext[0] == 'v'){
                            if(sscanf(yytext,"view%d:R=%d,%d,%d;T=%lf,%lf,%d."
                                ,&i,&rx,&ry,&rz,&tx,&ty,&tz)!= 7){
					print_error("input error");
			    }
			  } else { 
                            if(sscanf(yytext,"file%d:R=%d,%d,%d;T=%lf,%lf,%d."
                                ,&i,&rx,&ry,&rz,&tx,&ty,&tz)!= 7){
					print_error("input error");
			    }
			  }
			  if(i==KEY_FILE){
				NEW(view,10,Int4);
			        view[0]=rx; view[1]=ry; view[2]=rz; 
				view[3]=(Int4) (100.0* tx);
				view[4]=(Int4) (100.0 *ty);
				view[5]=tz;
				yylval.Vval=view;
				return VIEW;
			  }
			}

	/****************** Molecule information ***********************/
^water{I}[=].*$ { BEGIN 0; 
		}
<ISMOL>{CH}/[_@:.~^+]	{ BEGIN 0; yylval.cval=yytext[0]; if(CURRENT_FILE==KEY_FILE) return CHAIN; }

{L}/[_@:.~^+]	{ BEGIN 0; yylval.cval=yytext[0]; if(CURRENT_FILE==KEY_FILE) return CHAIN; }

<SQUIGGLE>{CS}	|  
<PLUS_SIGN>{CS}	|  
<COLON>{CS}	|  
<DOT>{CS}	|  
<ATSIGN>{CS}	|  
<CIRCUMFLEX>{CS}	|  
<LEFT_BRACKET>{CS} { 
		BEGIN 0; yylval.cval=yytext[0]; if(CURRENT_FILE==KEY_FILE) return COLOR; 
		}
<LEFT_PARAN>{CS} { BEGIN 0; yylval.cval=yytext[0]; if(CURRENT_FILE==KEY_FILE) return COLOR; }

{M}{A}		{  BEGIN 0; }

<INITIAL>{A}/{I} {  BEGIN 0; yylval.cval=yytext[0]; 
		if(CURRENT_FILE==KEY_FILE) return RESIDUE; }

[!]{MOL}	{ BEGIN ISMOL; 
		  NEW(yylval.pcval,yyleng+4,char); // fprintf(stderr,"%s\n",yytext);
		  yytext++; sscanf(yytext,"%s",yylval.pcval);
		  if(CURRENT_FILE==KEY_FILE) return MOLECULE;
		}

<UNDERLINE>[c][=][o]/[.] { BEGIN 0; if(CURRENT_FILE==KEY_FILE) return CARBONYL; }

<UNDERLINE>{ATOM}[\*]	|
<UNDERLINE>{ATOM}	{ BEGIN 0; sscanf(yytext,"%s",yylval.aval); 
				if(CURRENT_FILE==KEY_FILE) return AN_ATOM; }

<DASH>{ATOM}[\*]	|
<DASH>{ATOM}		{ BEGIN 0; sscanf(yytext,"%s",yylval.aval); 
				if(CURRENT_FILE==KEY_FILE) return HYDROGEN; }

<ISMOL>{I}/{CH}[_@:.~^+]	{ 
	  sscanf(yytext,"%d",&yylval.ival); // fprintf(stderr,"%s\n",yytext);
	  if(CURRENT_FILE==KEY_FILE) return(INTEGER);
	}

{I}	{ BEGIN 0; 
	  sscanf(yytext,"%d",&yylval.ival); // fprintf(stderr,"%s\n",yytext);
	  if(CURRENT_FILE==KEY_FILE) return(INTEGER);
	}

"~$="	{ 
		BEGIN CLOSED; 	// skip the rest of the input file.
		// return END_OF_FILE; 
	}
	/****************** single character tokens ***********************/
"."	{ BEGIN DOT; if(CURRENT_FILE==KEY_FILE) return('.');}
"!"	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return('!');}
"*"	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return('*');}
"~"	{ BEGIN SQUIGGLE; if(CURRENT_FILE==KEY_FILE) return('~');}
"+"	{ BEGIN PLUS_SIGN; if(CURRENT_FILE==KEY_FILE) return('+');}
"{"	{ BEGIN LEFT_BRACKET; if(CURRENT_FILE==KEY_FILE) return('{');}
"}"	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return('}');}
","	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return(',');}
"@"	{ BEGIN ATSIGN; if(CURRENT_FILE==KEY_FILE) return('@');}
"^"	{ BEGIN CIRCUMFLEX; if(CURRENT_FILE==KEY_FILE) return('^');}
"-"	{ BEGIN DASH; 
		// fprintf(stderr,"%s\n",yytext);
		if(CURRENT_FILE==KEY_FILE) return('-'); }
"_"	{ BEGIN UNDERLINE; if(CURRENT_FILE==KEY_FILE) return('_');}
":"	{ BEGIN COLON; if(CURRENT_FILE==KEY_FILE) return(':');}
"="	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return('=');}
"/"	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return('/');}
["]	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return('"');}
"("	{ BEGIN LEFT_PARAN; if(CURRENT_FILE==KEY_FILE) return('(');}
")"	{ BEGIN 0; if(CURRENT_FILE==KEY_FILE) return(')');}

[ \t]	{ BEGIN 0; }
.	{ BEGIN 0;
		fprintf(stderr,"Invalid token (line %d): '%s'\n",PARSELINE,yytext);
		yylval.cval=yytext[0];
		return JUNK;}

<<EOF>>	{ return END_OF_FILE; }

[\n]	{ BEGIN 0; PARSELINE++; 
		if(VERBOSE && PrintComment) { PrintComment=FALSE; fprintf(stderr,"\n"); }
		if(Implicit) CURRENT_FILE=KEY_FILE; else CURRENT_FILE=0;
		return NEWLINE; } 
%%

int yywrap() { BEGIN 0; return(1); }

