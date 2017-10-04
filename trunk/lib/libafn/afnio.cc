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

#include "afnio.h"

FILE	*open_file(const char *fstring,const char *subfile,const char *cmnd)
{
	FILE	*fptr;
	char	s[502];

	s[500]=0;
	while(fstring[0] == ' ') fstring++;
	strcpy(s,fstring); strcat(s,subfile);
	if(s[500] != 0) print_error("File name length is > 500 characters...aborting");
	if((fptr = fopen(s,cmnd)) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",s);
		print_error("File does not exist!\n");
	}
	return(fptr);
}

Int4     ParseReals(char *str, double *values, const char *msg)
/** WARNING: index starts at 0 not 1 as for ParseIntegers() **/
{
        Int4	n;
	double	k;

	if(!isdigit(str[0]) && str[0] != '-') print_error(msg);
        for(n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0]) || str[0] == '-'){
                if(sscanf(str,"%lf", &k) != 1){ print_error(msg); }
                else { 
		   values[n-1] = k; 
		   if(str[0] == '-') str++;
		   while(isdigit(str[0]) || str[0] == '.') str++; 
		}
           } else print_error(msg);
        }
        return n;
}

Int4     ParseRegions(char *str, Int4 **start, Int4 **end, const char *msg)
// input string: "5..7,9..17,30..55" 
// returns:      3 (regions) and sets start=[0,5,9,30,0] & end=[0,7,17,55,0]
//                                    0 1 2 3 ...
{
        Int4     N,n,w,*Start,*End,S,E;

	if(!isdigit(str[0])) print_error(msg);
	N=strlen(str); N = (N/2) + 1;
	NEW(Start,N+3,Int4); NEW(End,N+3,Int4);
        for(n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%d..%d", &S,&E) != 2) print_error(msg);
                else if(S > E) print_error(msg);
                else if(S < End[n-1]) print_error(msg);
                else if(E < Start[n-1]) print_error(msg);
                else { 
		   Start[n] = S; End[n] = E;
		   while(str[0] && str[0] != ',') str++; 
		}
           } else print_error(msg);
        } *start=Start; *end=End;
	return n;
}

Int4     ParseIntegers(char *str, Int4 *values, const char *msg)
// input string: "3,5-7,9,11-17" 
// returns:      12 and sets values=[12,3,5,6,7,9,11,12,13,14,15,16,17]
//                                    0 1 2 3 ...
{
        Int4     n,v,w;

	if(!isdigit(str[0])) print_error(msg);
        for(n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%d", &v) != 1){ print_error(msg); }
                else { 
		   values[n] = v; while(isdigit(str[0])) str++; 
		   if(str[0] == '-'){
			str++; if(!isdigit(str[0])) print_error(msg);
			if(sscanf(str,"%d", &w) != 1) print_error(msg); 
			if(w <= v) print_error(msg);
			for(Int4 i=v+1; i <= w; i++){ n++; values[n] = i; }
			while(isdigit(str[0])) str++;
		   }
		}
           } else print_error(msg);
        } return n;
}

double	RealOption(char *argv, char option, double min, double max, const char *usage)
{
        double	d;
	char	str[50];

	if(argv[0] != '-' || !isalpha(argv[1])) print_error(usage); 
	if(argv[2] == '=') sprintf(str,"-%c=%%lf",option);
	else sprintf(str,"-%c%%lf",option);
	if(sscanf(argv,str, &d) != 1) print_error(usage);
	if(d < min || d > max) print_error(usage);
        return d;
}

Int4	IntOption(char *argv, char option, const Int4 min, const Int4 max, const char *usage)
{
        Int4     n;
	char	str[50];

	if(argv[0] != '-' || !isalpha(argv[1])) print_error(usage); 
	if(argv[2] == '=') sprintf(str,"-%c=%%d",option);
	else sprintf(str,"-%c%%d",option);
	if(sscanf(argv,str, &n) != 1) print_error(usage);
	if(n < min || n > max) print_error(usage);
        return n;
}

void    print_bits(FILE *fptr, char type, UInt4 n)
{
        Int4    i,end;

	if(type == 'l') end = 31;
	else if(type == 's') end = 15;
	else if(type == 'b') end = 7;
	else print_error("print_bits( ) type error");
        for(i=end; i>=0; i--){
           if(n & (1 << i)){ fprintf(fptr,"1"); }
           else fprintf(fptr,"0");
        } fprintf(fptr,"\n"); 
}

Int4	string2argv(char *argv[], char *s)
/*** convert the string s into argv format. ****/
/*** WARNING: assumes argv array is large enough!! ****/
{
	Int4	n,i,argc;

	for(argc=0; *s != 0; ){
	   while(isspace(*s)) s++;
	   if(*s != 0){
	     for(n=i=0; s[i] != 0 && !isspace(s[i]); i++) n++;
	     NEW(argv[argc],n+2,char); 
	     for(i=0 ; i < n; s++,i++) argv[argc][i] = *s;
	     argc++; 
	   }
	}
	argv[argc] = NULL;
	return argc;
}

char	*AllocString(const char *s)
{
	char *S;
	if((S=(char*) calloc((strlen(s)+2),sizeof(char)))==NULL)
		print_error("Out of Memory");
	strcpy(S,s); return S; 
}

char	*NewString(const char *s)
{
	Int4	n=strlen(s);
	char *S = new char[n+2];
	if(!S) print_error("Out of Memory");
	strcpy(S,s); return S; 
}
