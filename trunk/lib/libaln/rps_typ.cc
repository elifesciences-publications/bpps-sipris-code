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

#include "rps_typ.h"

rps_typ::rps_typ(FILE *fp)
{
	init( );
	char	str[200];
	assert(fscanf(fp,"<Domains(%u)=\n",&NumSeqs) == 1);
	domains = new pdl_typ*[NumSeqs+2];
	sq_ids = new char*[NumSeqs+2];
	NumDoms = new unsigned short[NumSeqs+2];
	for(Int4 s=0; s <= NumSeqs; s++){ domains[s]=0; sq_ids[s]=0; }

	UInt4   strt_sq,end_sq;
        UInt4	strt_dom,end_dom,len_dom,extend;
	float	evalue;
	Int4	i=1,err,SqNum;
	char	c;
	do {
	    err = fscanf(fp,"->(%d)%[^=]={",&SqNum, str); 
	    if(err != 2){
		if(fscanf(fp,">.\n") == 0) break;
		else print_error("rps_typ input error");
	    } 
	    assert(SqNum > 0 && SqNum <= NumSeqs);
	    assert(i == SqNum);
	    sq_ids[i] = AllocString(str);
	    Int4 J=0;
	    do {
	      if((c=getc(fp)) == '}'){
		ungetc(c,fp);
		if(fscanf(fp,"};\n") == 0) break; 
		else print_error("rps_typ input error");
	      } 
	      ungetc(c,fp);
              err=fscanf(fp,"[%u-%u->%[^:]:%u-%u(%u);%f]",
                   &strt_sq,&end_sq,str,&strt_dom,&end_dom,&extend,&evalue);
                // evalue == inf???
	      if(err != 7){ print_error("rps_typ input error"); }
	      J++;
	      len_dom=end_dom+extend;
	      if(domains[i]==0){
		domains[SqNum] = 
		  new pdl_typ(strt_sq,end_sq,strt_dom,end_dom,len_dom,str,evalue);
	      } else {
		domains[SqNum]->Append(strt_sq,end_sq,strt_dom,end_dom,
			len_dom,str,evalue);
	      }
	    } while(TRUE);
	    NumDoms[SqNum] = J;
	    i++;
	} while(TRUE);
}

Int4    rps_typ::SqStart(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->SqStart(d);
	else return 0; // this should not happen.
}

Int4    rps_typ::SqEnd(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->SqEnd(d);
	else return 0; // this should not happen.
}

char    *rps_typ::DomainName(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->DomainName(d);
	else return 0; // this should not happen.
}

float	rps_typ::Evalue(UInt4 s, unsigned short d)
{
	if(s < 1 || s > NumSeqs || d < 1 || d > NumDoms[s]) return 0;
	else if(domains[s]) return domains[s]->Eval(d);
	else return 0; // this should not happen.
}

Int4	rps_typ::NumDomains(UInt4 s)
{ if(s < 1 || s > NumSeqs) return 0; else return NumDoms[s]; }

void    rps_typ::Put(FILE *fp, BooLean *skip, char *skipdom)
{
	Int4	s,S,n;
	if(skip){
	   for(n=0,s=1; s <= NumSeqs; s++) if(!skip[s]) n++;
	} else n=NumSeqs;
	fprintf(fp,"<Domains(%d)=\n",n);
	for(s=0,S=1; S <= NumSeqs; S++){
	   if(skip && skip[S]) continue;
	   else s++;
	   fprintf(fp,"->(%d)%s={",s,sq_ids[S]);
	   if(domains[S]) domains[S]->Put(fp,skipdom);
	   fprintf(fp,"};\n");
	} fprintf(fp,">.\n");
}

void    rps_typ::Free()
{
	for(Int4 s=1; s <= NumSeqs; s++){ 
	   if(domains[s]) delete domains[s];
	   if(sq_ids[s]) free(sq_ids[s]);
	}
	if(sq_ids) delete sq_ids;
	if(domains) delete domains;
	if(NumDoms) delete NumDoms;
	init();
}

void    rps_typ::init( )
{ domains = 0; sq_ids = 0; NumSeqs=0; NumDoms=0; }

rps_typ& rps_typ::operator=(const rps_typ& dom)
// called for rps_typ dom2; dom2=dom;
{ if (this != &dom) { Free(); copy(dom); } return *this; }

void    rps_typ::copy(const rps_typ& dom) 
// private function that assumes 'this' has been initialized and
// copies dom to 'this'.
{
	NumSeqs=dom.NumSeqs;
	
	domains = new pdl_typ*[NumSeqs+2];
	sq_ids = new char*[NumSeqs+2];
	NumDoms = new unsigned short[NumSeqs+2];
	for(Int4 s=1; s <= NumSeqs; s++){
	   if(dom.domains[s]){ 
		domains[s] = new pdl_typ;
		domains[s]->Copy(*(dom.domains[s]));
		// domains[s] = dom.domains[s]; 
	   }
	   else domains[s] = 0;
	   if(dom.sq_ids[s]){ sq_ids[s] = AllocString(dom.sq_ids[s]); }
	   else sq_ids[s]=0; 
	   NumDoms[s] = dom.NumDoms[s];
	}
}

