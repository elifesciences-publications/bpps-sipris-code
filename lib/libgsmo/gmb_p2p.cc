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

#include "gmb_typ.h"

//************************ Single Sequence Sampling Routines ***********************
char	*gmb_typ::AlignP2P_Operation(cma_typ child,Int4 &start)
{
	Int4	i,trace_length,score,*newpos,*sites;
	cma_typ CMA=SSX->RtnCMA( ); 
	char	dm=SSX->RtnDmsMode();
	gssCMSA(child)->FakeSqSet(); // creates fake seq set...
	ssx_typ *ssx = new ssx_typ(500,5,500,2,(double)pernats,child,dm);
	// cma_typ csq_cma=AddConsensusCMSA(child);
	gss_typ *gss=gssCMSA(CMA);
	char *operation= SSX->AlignP2P(ssx,start,score);
	delete ssx;
	return operation;
}

BooLean	gmb_typ::AlignP2P(FILE *fp,cma_typ child)
// Optimally align two profiles.
{
	Int4	i,start,trace_length,score,*newpos,*sites;
	cma_typ CMA=SSX->RtnCMA( ); 
	char	dm=SSX->RtnDmsMode();
	gssCMSA(child)->FakeSqSet(); // creates fake seq set...
	ssx_typ *ssx = new ssx_typ(500,2,500, 1,(double)pernats,child,dm);
	// don't use gap penalties, so can ignore...
	// cma_typ csq_cma=AddConsensusCMSA(child);

	gss_typ *gss=gssCMSA(CMA);
	cma_typ rcma=EmptyCMSA(nBlksCMSA(CMA),LengthsCMSA(CMA),TrueDataCMSA(CMA),
			gss->GapOpen(),gss->GapExtend(),PerNatsCMSA(CMA),0,0);

	char *operation= SSX->AlignP2P(ssx,start,score);

	if(score <= 0){
		fprintf(stderr,"score = %d; skipping...\n",score);
		return FALSE;
	}
	e_type  csqE=MkConsensusCMSA(child);
fprintf(stderr,"operation = %s\n",operation);
 this->PutGappedAln(stderr,csqE,operation,start);
fprintf(stderr,"score = %d; operation = %s\n",score,operation);
#if 0
	// PutSeq(stderr,csqE,AB);
	gsq_typ *gsq0 = new gsq_typ[1];
	NEW(sites,nBlksCMSA(CMA)+3,Int4);
	gsq0->initialize(operation,csqE,sites);
	free(operation);
	
	fprintf(stderr,"sites[1]= %d\n",sites[1]);
	// gsq0->Put_cma_format(fp,1,nBlksCMSA(CMA),sites,LengthsCMSA(CMA),AB);
	// Put_cma_format(FILE *fp,Int4 J0,Int4 nblk,Int4 *start, Int4 *blk_len,a_type A);
	delete []gsq0;
#endif

operation=SSX->GapAlnTrace(csqE,0,start,score);
 this->PutGappedAln(stderr,csqE,operation,start);
fprintf(stderr,"score = %d; operation = %s\n",score,operation);
	free(operation);
#if 0
	ss_type data=SeqSet1(NameCMSA(child),csqE,AB);
	cma_typ xcma=this->CreateAlign(data);
	// PutCMSA(stderr,xcma); 
	TotalNilCMSA(xcma);
#endif
	return TRUE;
}

