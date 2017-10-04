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

#include "rsq_typ.h"

Int4	PutBestRepsCMSA(FILE *fp, Int4 num_phyla, BooLean IgnoreGaps, cma_typ cma)
{ PutBestRepsCMSA(fp, num_phyla, IgnoreGaps, cma,FALSE,0); }

Int4	PutBestRepsCMSA(FILE *fp, Int4 num_phyla, BooLean IgnoreGaps, cma_typ cma,BooLean UseAllPDB,
	cma_typ csq_cma)
{
	dh_type dH=NULL;
	Int4    hits,KeySeq=0,Score,I,J,N=NumSeqsCMSA(cma); // set KeySeq == to the best sequence...
	assert(num_phyla > 0);
	if(nBlksCMSA(cma) != 1) print_error("this option requires only one blk");
	if(N < 1){ print_error("no sequences in cma file"); }
	dH=dheap(N+2,4);
	if(csq_cma == 0){
	    double prob,best=-9999999999.0;
	    for(J=1; J <= N; J++){
       		if(IgnoreGaps) prob = GetProbCMSA(1,J,cma);
		else prob = GetGappedProbCMSA(1,J,cma);
		if(prob > best) { best=prob; KeySeq=J; }
	    }
	    // put out best in each phylum only...
	    // fprintf(stderr,"KeySeq=%d\n",KeySeq);
	    Int4 SelfScore=PseudoAlnScoreCMSA(KeySeq, KeySeq,cma);
	    insrtHeap(KeySeq,-((keytyp)SelfScore + 0.01),dH); 
	    // fprintf(stderr,"KeySeq=%d; Self Score = %d\n",KeySeq,Score);
	} BooLean *skip; NEW(skip,N+3,BooLean); 
	for(J=1; J <= N; J++){
		skip[J]=TRUE;
		if(csq_cma){ Score = PseudoAlnScoreTwoCMSA(1,csq_cma,J,cma); }
		else if(J==KeySeq) continue; else Score=PseudoAlnScoreCMSA(KeySeq,J,cma);
		e_type sE=SeqSetE(J,TrueDataCMSA(cma));
		if(Score <= 0) continue;
		if(UseAllPDB && PdbSeq(sE)){ Score = Score*10; }
		if(UseAllPDB && EST_Seq(sE)){ Score = Score/2; }
		insrtHeap(J,-(keytyp)Score,dH); 
	}
	Int4	NumPhyla=0,*phyla=GetPhylaSeqSet(0, &NumPhyla, TrueDataCMSA(cma));
	// phyla=GetPhylaSeqSet(stderr, &NumPhyla, TrueDataCMSA(cma));
	BooLean	*Printed=0; NEW(Printed,NumPhyla+3,BooLean);
	Int4 *list; NEW(list, N+3,Int4); 
	for(hits=0,I=0; I < num_phyla && !emptyHeap(dH); ){
		assert((J=delminHeap(dH)) != 0);
		if(!Printed[phyla[J]]) {
			skip[J]=FALSE; Printed[phyla[J]]=TRUE;
			I++; list[I]=J; hits++;
		}
	} Nildheap(dH); 
	for(J=1; J <= N; J++){ if(skip[J]){  I++; assert(I <= N); list[I]=J; } }
	assert(I == N);
	PutSelectOneCMSA(fp,skip,list,cma);
	free(list); free(skip); free(phyla); free(Printed);
	return hits;
}

void	PutBestCMSA(FILE *fp, Int4 put_best, cma_typ cma)
{
	a_type	AB=AlphabetCMSA(cma);
	Int4	N,I,J,NumPhyla=0,NumHits,*phyla;
	BooLean	*Printed=0,*skip;
	double	score;

	N=NumSeqsCMSA(cma);
	dh_type dH=dheap(N+2,4);
	if(nBlksCMSA(cma) != 1) print_error("this option requires only one blk");
	if(N < 1){ print_error("cdh_typ.cc error: no sequences in cma file"); }
	NEW(skip,N+3,BooLean); 
	NumPhyla=0,NumHits;Printed=0;
	// phyla=GetPhylaSeqSet(stderr,&NumPhyla,TrueDataCMSA(cma));
	phyla=GetPhylaSeqSet(0,&NumPhyla,TrueDataCMSA(cma));

	e_type csqE = MkConsensusCMSA(cma);	// for alternative method below.
	h_type HG=Histogram("pseudoscores", 0,5000,25.0);
	for(J=1; J <= N; J++){
#if 0	// insert prob score..
               	score= GetProbCMSA(1,J,cma);
		// prob = GetGappedProbCMSA(1,J,cma);
#else	// insert scort to consensus sequence.
		// e_type sE=TrueSeqCMSA(J,cma);
		// ResidueCMSA(register Int4 t, register Int4 n, register Int4 s, cma_typ cma);
		// PseudoAlnScoreSqToCMSA(e_type E, Int4 sq2, cma_typ cma);
		score = PseudoAlnScoreSqToCMSA(csqE,J,cma);
#endif
		insrtHeap(J,-(keytyp)score,dH); 
		IncdHist((double)score,HG);
		skip[J]=TRUE;
	}
	NilSeq(csqE);
        PutHist(stderr,50,HG); NilHist(HG); 
	Int4 *list; NEW(list, N+3,Int4); NEW(Printed,NumPhyla+3,BooLean);
	for(NumHits=I=0; I < put_best && !emptyHeap(dH); ){
                assert((J=delminHeap(dH)) != 0);
		if(!Printed[phyla[J]]) {
		    skip[J]=FALSE; NumHits++;
		    // fprintf(stderr,"Seq Num: %d\n",J);
		    Printed[phyla[J]]=TRUE;
		    I++; list[I]=J;
		}
	} Nildheap(dH); 
	//************ put out best in each phylum only... ***********
#if 0
		for(J=1; J <= N; J++){
			e_type sE=TrueSeqCMSA(J,cma);
			if(DistinctPhyla){
			  if(KingdomSeq(sE) == 'U') continue;	// kingdom is unknown...
			  if(KingdomSeq(sE) == 'X') continue;	// kingdom is unknown...
			  if(KingdomSeq(sE) == 0) continue;	// kingdom is unknown...
			}
			Score=PseudoAlnScoreCMSA(KeySeq, J,cma);
			// IncdHist((double)Score,HG);
			insrtHeap(J,-(keytyp)Score,dH); 
		}
#endif
        // PutHist(stdout,50,HG); NilHist(HG); 
	for(I=0, J=1; J <= N; J++){ if(skip[J]){  I++; assert(I <= N); list[I]=J; } }
	// assert(I == N);
	// if(NumHits > 1) PutSelectOneCMSA(fp,skip,list,cma);
	// if(NumHits > 1) PutSelectOneCMSA(stdout,skip,list,cma);
	PutSelectCMSA(fp,skip,cma); 
	free(list); free(skip); free(phyla); free(Printed);
	return;
}

