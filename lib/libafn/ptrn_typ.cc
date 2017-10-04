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

#include "ptrn_typ.h"

ptrn_typ MakePtrn(Int4 N, Int4 Len, sst_typ **sst, double **lpr)
{
	ptrn_typ ptrn;
	MEW(ptrn,1,pattern_type);
	ptrn->N=N;
	ptrn->Len=Len;
	ptrn->SST=sst;
	ptrn->LPR=lpr;
	return ptrn;
}

void	PutPtrn(FILE *fp,ptrn_typ ptrn, Int4 i, a_type AB)
{
	sst_typ	sst;
	Int4	j,r;
	assert(i > 0 && i <= ptrn->N);
	if(ptrn->SST[i] != 0){ 
	// fprintf(fp,"%d: \n",i);
	  for(j=1; j<=ptrn->Len; j++){
	    sst=ptrn->SST[i][j];
	    if(sst){ // then print out common residue(s):
		for(r=1; r <= nAlpha(AB); r++){
		   if(MemSset(r,sst)){ fprintf(fp,"%c",AlphaChar(r,AB)); }
		} fprintf(fp,"%d ",j);
	        // fprintf(fp,"  %.1f:\n",ptrn->LPR[i][j]);
	    }
	  }
	} fprintf(fp,"\n");
}

void	PutPtrn(FILE *fp,ptrn_typ ptrn, a_type AB)
{ for(Int4 i=1; i<=ptrn->N; i++){ fprintf(fp,"column %d: \n",i); PutPtrn(fp,ptrn,i,AB); } }

void	WritePtrn(FILE *fp, ptrn_typ ptrn)
{
        Int4    i,rtn;

        rtn=fwrite(ptrn,sizeof(pattern_type),1,fp);
        if(rtn != 1) print_error("WritePtrn() input error 1");
#if 0	// this stuff is written above...
        rtn=fwrite(&ptrn->N,sizeof(Int4),1,fp);
        if(rtn != 1) print_error("WritePtrn() input error 2");
        rtn=fwrite(&ptrn->Len,sizeof(Int4),1,fp);
        if(rtn != 1) print_error("WritePtrn() input error 3");
#endif
	for(i=1; i<=ptrn->N; i++){
	  // if(ptrn->LPR[i] == 0){ continue; }
          rtn=fwrite(ptrn->LPR[i],sizeof(double),ptrn->Len+2,fp);	// made in cmc_typ
          if(rtn != ptrn->Len+2) print_error("WritePtrn() input error 4");
          rtn=fwrite(ptrn->SST[i],sizeof(sst_typ),ptrn->Len+2,fp);	// made in cmc_typ
          if(rtn != ptrn->Len+2) print_error("WritePtrn() input error 4");
	}
}

void    NilPtrn(ptrn_typ ptrn)
{
        Int4    i,rtn;
	for(i=1; i<=ptrn->N; i++){ free(ptrn->LPR[i]); free(ptrn->SST[i]); }
	free(ptrn->LPR); free(ptrn->SST);
	free(ptrn);
}

ptrn_typ ReadPtrn(FILE *fp)
{
        Int4    i,rtn;
	ptrn_typ ptrn;

	MEW(ptrn,1,pattern_type);
        rtn=fread(ptrn,sizeof(pattern_type),1,fp);
        if(rtn != 1) print_error("ReadPtrn() input error 1");
#if 0	// this stuff is written above...
        rtn=fread(&ptrn->N,sizeof(Int4),1,fp);
        if(rtn != 1) print_error("ReadPtrn() input error 2");
        rtn=fread(&ptrn->Len,sizeof(Int4),1,fp);
        if(rtn != 1) print_error("ReadPtrn() input error 3");
#endif
	NEWP(ptrn->LPR,ptrn->N+3,double);
	NEWP(ptrn->SST,ptrn->N+3,sst_typ);
	for(i=1; i<=ptrn->N; i++){
	  NEW(ptrn->LPR[i],ptrn->Len+3,double);
	  NEW(ptrn->SST[i],ptrn->Len+3,sst_typ);
          rtn=fread(ptrn->LPR[i],sizeof(double),ptrn->Len+2,fp);	// made in cmc_typ
          if(rtn != ptrn->Len+2) print_error("ReadPtrn() input error 4");
          rtn=fread(ptrn->SST[i],sizeof(sst_typ),ptrn->Len+2,fp);	// made in cmc_typ
          if(rtn != ptrn->Len+2) print_error("ReadPtrn() input error 4");
	}
	return ptrn;
}

