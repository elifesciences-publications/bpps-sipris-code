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

#include "cmsa.h"

Int4    *FirstRptFullCMSA(cma_typ cma)
{ assert(cma && cma->FullSeq); return cma->FirstRpt; }

Int4    *first_rpt_full_cmsa(cma_typ cma)
// return an array with the number of the subseq first 
{
	Int4	*first,sq,lsq;
	assert(cma && cma->FullSeq);
	NEW(first, NSeqsSeqSet(cma->FullSeq) + 3, Int4);
	first[1] = 1;
	for(lsq=1,sq=2; sq <= NSeqsSeqSet(cma->FullSeq); lsq++,sq++){
	   if(cma->FullRpts[sq]) first[sq]=cma->FullRpts[lsq]+first[lsq];
	} return first;
}

unsigned short *FullRptsCMSA(cma_typ cma)
{ assert(cma && cma->FullSeq); return cma->FullRpts; }

Int4    RptNumCMSA(Int4 sq, cma_typ cma)
// return the number of the repeat in full seq corresponding to sq.
{
	if(sq==0) return 0;
	assert(cma && cma->FullSeq && sq > 0 && sq <= NumSeqsCMSA(cma));
	Int4 fsq,sq0;
	fsq = SubToFullCMSA(sq,cma);
	sq0 = cma->FirstRpt[fsq];
	return (sq-sq0+1);
}

Int4    FullRptsCMSA(Int4 sq, cma_typ cma)
// return number of repeats in full sequence sq.
{
	assert(cma && cma->FullSeq && sq > 0 && sq <= NSeqsSeqSet(cma->FullSeq));
	return cma->FullRpts[sq];
}

void    CopyFullCountsCMSA(cma_typ cma2, cma_typ cma)
// Copy full counts info from cma to cma2;
{
   RmFullCountsCMSA(cma2);
   if(cma->FullSeq){
	Int4 N=NSeqsSeqSet(cma->FullSeq);
	NEW(cma2->FullRpts, N + 3, unsigned short);
	NEW(cma2->FirstRpt, N + 3, Int4);
	for(Int4 n=1; n <= N; n++){
	   cma2->FullRpts[n]=cma->FullRpts[n];
	   cma2->FirstRpt[n]=cma->FirstRpt[n];
	}
	cma2->FullSeq = CopySeqSet(cma->FullSeq);
	ss_type	data = TrueDataCMSA(cma);
	Int4 S = NSeqsSeqSet(data); 
	NEW(cma2->SubToFull, S+3, Int4);
	for(Int4 s=1; s<=S; s++){ cma2->SubToFull[s] = cma->SubToFull[s]; }
   }
}

void    AddFullCountsCMSA(ss_type FullSeq, unsigned short  *FullRpts, cma_typ cma)
// Add information about full sequence counts when subseqs are used.
// WARNING: eventually need to replace with full sequence representations.
{
	assert(cma && FullSeq && FullRpts);
	if(cma->FullSeq) RmFullCountsCMSA(cma);
	cma->FullSeq = FullSeq; cma->FullRpts = FullRpts;
	Int4	r,f,s,S,F;
	ss_type	data = TrueDataCMSA(cma);
	S = NSeqsSeqSet(data); F = NSeqsSeqSet(FullSeq);
	NEW(cma->SubToFull, S +3, Int4);
	for(s=0,f=1; f<=F; f++){
	    for(r=FullRpts[f]; r > 0; r--){
		s++; cma->SubToFull[s] = f; 
#if 0	// debug: 
	   	fprintf(stderr,"%d (%d): ",f,r);
		e_type E = SeqSetE(s,data);
		PutShortSeqID(stderr,E);
	   	fprintf(stderr,"\n");
#endif
	    }
	   	// fprintf(stderr,"\n");
	} 
	if(s != S){
	   fprintf(stderr,"s = %d != S = %d\n",s,S);
	   fprintf(stderr,"May need to remove concensus sequence\n");
	   print_error("Fatal: fullseq information inconsistent with cma file.\n");
	   assert(s == S);
	}
	cma->FirstRpt = first_rpt_full_cmsa(cma);
}

void    RmFullCountsCMSA(cma_typ cma)
{
        if(cma->FullSeq) NilSeqSet(cma->FullSeq); cma->FullSeq=0;
	if(cma->FullRpts) free(cma->FullRpts); cma->FullRpts=0;
	if(cma->SubToFull) free(cma->SubToFull); cma->SubToFull=0;
	if(cma->FirstRpt) free(cma->FirstRpt); cma->FirstRpt=0;
	if(cma->Domains) RmDomainsCMSA(cma);
}

Int4    SubToFullCMSA(Int4 sq, cma_typ cma)
{
	assert(cma->SubToFull != 0);
	ss_type data=DataCMSA(cma);
	Int4 N=NSeqsSeqSet(data);
	assert(sq > 0 && sq <= N);
	return cma->SubToFull[sq];
}

