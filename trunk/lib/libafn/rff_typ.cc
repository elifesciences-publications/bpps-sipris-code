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

#include "rff_typ.h"

rff_typ::rff_typ( ) { init( ); }

rff_typ::rff_typ(FILE *fptr) { init( ); Read(fptr); }

void    rff_typ::init( ) { len_ss=0; ss=0; mode=0; }

void    rff_typ::Free( ) { if(ss) free(ss); }

rff_typ::rff_typ(const rff_typ *rff,UInt4 start,UInt4 end)
{ init(); copy2(rff,start,end); }

rff_typ& rff_typ::operator=(const rff_typ& rff) // called for rff_typ rff2; rff2=rff;
{ if (this != &rff) { Free(); init(); copy(rff); } return *this; }

void    rff_typ::copy2(const rff_typ *rff) { copy2(rff,1,rff->len_ss); }

void    rff_typ::copy2(const rff_typ *rff,UInt4 start,UInt4 end) 
// private function that assumes 'this' has been initialized and
// copies rff to 'this'.
{
	if(rff->ss){
	  Int4	i,j;
	  assert(start < end && rff->len_ss >= end);
	  mode=rff->mode;
          len_ss=end-start+1;
          NEW(ss,len_ss+3,char);
	  for(j=1,i=start; i <= end; i++,j++) ss[j]=rff->ss[i];
	}
}

void    rff_typ::copy(const rff_typ& rff) 
// private function that assumes 'this' has been initialized and
// copies rff to 'this'.
{
	if(rff.ss){
          len_ss=rff.len_ss;
	  mode=rff.mode;
          NEW(ss,len_ss+3,char);
	  for(Int4 i=0; i <= len_ss; i++) ss[i]=rff.ss[i];
	}
}

rff_typ::rff_typ(char *str) { init( ); Read(str); }

void    rff_typ::Read(char *str)
{
        char    c,tmpstr[4]; 
        UInt4    i,j=0,len;

	c=str[j]; j++;
	if(c != '+') print_error("rff_typ() input error");
	mode=str[j]; j++;
        switch(mode){
          case 'S':     // secondary structure information.
                if(sscanf(str,"+S(%u)=",&len) != 1)
                   print_error("rff_typ() file format error");
                else {
                  assert(ss==0 && len > 0);
		  tmpstr[1]=0;
		  while(str[j] != '=') j++; j++;
                  len_ss=len; 
                  NEW(ss,len+3,char);
                  for(i=1; i <= len; i++){
                        c=str[j]; j++; 
			tmpstr[0]=c;
                        if(strstr("xHhSsCcu",tmpstr)==0) print_error("rff_typ() format error");
                        ss[i]=c;
                  }
                  c=str[j]; j++;; if(c != ';') print_error("rff_typ() format error");
                }
           break;
          default: print_error("rff_typ() format error"); break;
        } return;
}

void    rff_typ::Read(FILE *fptr)
// Just saw a '^+' in input stream; scan in information;
{
        char    c,str[3];
        UInt4    i,len;

	c=fgetc(fptr);
	if(c != '+') print_error("rff_typ() input error");
	mode=fgetc(fptr);
        switch(mode){
          case 'S':     // secondary structure information.
                if(fscanf(fptr,"(%u)=",&len) != 1)
                   print_error("rff_typ() file format error");
                else {
                  assert(ss==0 && len > 0);
                  len_ss=len; str[1]=0;
                  NEW(ss,len+3,char);
                  for(i=1; i <= len; i++){
                        c=fgetc(fptr); 
			if(isspace(c)){ i--; continue; }
			str[0]=c;
                        if(strstr("xHhSsCcu",str)==0) print_error("rff_typ() format error");
                        ss[i]=c;
                  }
                  c=fgetc(fptr); if(c != ';'){
			for( ; c!=EOF && c != '\n'; c=fgetc(fptr)) fprintf(stderr,"%c",c); 
			print_error("rff_typ() format error");
		  }
                  c=fgetc(fptr); if(c != '\n') print_error("rff_typ() format error");
// fprintf(stderr,"len = %d: ss = %s\n",len_ss,ss+1);
                }
           break;
          default: print_error("rff_typ() format error"); break;
        } return;
}

void    rff_typ::Put(FILE *fptr,UInt4 start,UInt4 end)
{
	assert(ss);
	assert(start < end && len_ss >= end);
	fprintf(fptr,"+S(%d)=",end-start+1);
	for(UInt4 i=start; i <= end; i++){
		if(i%50==0 && i < len_ss) fprintf(fptr,"\n");
		fprintf(fptr,"%c",ss[i]);
	} fprintf(fptr,";\n");
}

char	ScanOverRFF(FILE *fptr)
{
	char	c;
	while((c=fgetc(fptr))!=EOF){ if(c==';') break; }
	if(c == EOF) print_error("rff_typ() format error");
	c=fgetc(fptr); if(c != '\n') print_error("rff_typ() format error");
	return c;
}

