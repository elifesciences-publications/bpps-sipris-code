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

#if !defined (_AFNGBLAST_)
#define	_AFNGBLAST_

#include "stdinc.h"
#include "my_ncbi.h"
#include "my_posit.h"	// contains all the functions I need...
#include "seqset.h"	// contains all the functions I need...
#include "histogram.h"
#include "random.h"
#include "gblast.h"
#include "gpsiblst.h"
#include "smatrix.h"
#include "karlin.h"
#include "dsets.h"
#include "marg_prob.h"
#include "set_typ.h"
#include "cmsa.h"	// add cma_typ...

class gpsi_type {	// AFN BLASTPGP routines
public:
			gpsi_type(e_type,ss_type,double,double,double,Int4,Int4);
			gpsi_type(e_type,ss_type,double,double,double,Int4,Int4,char*);
			gpsi_type(e_type,ss_type,double,double,double,Int4,Int4,FILE*);
			gpsi_type(e_type,ss_type,double,double,double,Int4,Int4,FILE*,Int4);
			~gpsi_type( );
	sap_typ		Search(Int4,BooLean,char);
	sap_typ		Search(Int4 T,char mode){ return Search(T,FALSE,mode); }
	sap_typ		TrimmedSearch(Int4,double,char);
	sap_typ		SearchWithPrior(Int4 Threshold,double misalncut,char mode,
			                cma_typ mtfcma);
	double		*FastSearch(BooLean,Int4,Int4,char);
	double		*FastSearch(BooLean,Int4,Int4,char,Int4 *,Int4 IgnoreThis=0);
	void		ShowAlign(FILE *,Int4,Int4);
	void		ShowAlign(FILE *,Int4,Int4,sap_typ);
	void		PutHitList(FILE *,Int4);
	void		PutDistHits(FILE *,Int4);
	Int4		DatabaseSize(){ return dbsLength; }
	BooLean		ComputeMatrix(char *);
	BooLean		ComputeMatrix(FILE *,char *);
//********************* NEW: Marginal Probability PSI-BLAST ************************
	BooLean		ComputeMatrix(char *, double ***);
//********************* NEW: Marginal Probability PSI-BLAST ************************
//********************* NEW: Prior motif alignment PSI-BLAST ************************
	void    	AppendFakeResults(gpsi_type *fake_gpsi);
//********************* NEW: Prior motif alignment PSI-BLAST ************************
	BooLean		NotConverged( );
	void		PutSeqs(FILE *fp){ if(head) PutSeqsListGSAP(fp,head,AB); }
	void		PutSameSeqs(FILE *fp,Int4 MinLen,double IDs){
				if(head) PutSameSeqsListGSAP(fp,head,MinLen,IDs,AB); 
			}
	sap_typ		GetSAP(){ return (head); }
	BooLean		*GetHitList( ); // Return Boolean array of hits.
	BooLean		*GetSeqHitArray(Int4 N){ 
	   BooLean *skip;
	   NEW(skip,N+3,BooLean);
	   if(head){
             for(sap_typ sap=head; sap != 0; sap=sap->next){
	      e_type sE=SubjectSeqGSAP(sap);
	      Int4 i=SeqI(sE);
	      assert(i > 0 && i <= N);
	      skip[i]=TRUE;
             }
	   } return (skip); 
	}
	void		PutSeqsAsCMA(FILE *fp,Int4 left,Int4 right){
			  if(head) SeqAlignToCMA(fp,head,left,right,AB);}
			// NOTE: the above routine is not working correctly!
			// THIS NEEDS TO BE FIXED!!!
	void		PutMissedSeq(FILE *fp); // output rest of seqs.
	void		PutMergedSubSeqs(FILE *fp,Int4 lt, Int4 rt){
				if(head) PutMergedGSAP(fp,head,lt,rt,AB); 
			}
	void		PutSubSeqs(FILE *fp,Int4 lt, Int4 rt){
				if(head) PutSubSeqsListGSAP(fp,head,lt,rt,AB); 
			}
	void		Test( ) { Test('A'); }
	void		Test(char mode);
	void		SeeMatrix(FILE *fp);
	void		PutRasMol(FILE *, sap_typ);
	void		PutRasMol(FILE *fp){ PutRasMol(fp,head); }
	sap_typ		SampleAlign(Int4 **,sap_typ);
	sap_typ		SampleMatrix(char *);
	void		KeepQuite( ){ verbose=FALSE; }
	void		SpeakUp( ){ verbose=TRUE; }
	Int4		ThisPass( ) { return thisPassNum; }
	Int4		DelimitSeqAlign(Int4 **SS,Int4 **SP,Int4 **ES,Int4 **EP);
	Int4		**posMatrix;		// for psi-blast
#if 1	// TEST ROUTINES
	void		TweakAlignment(sap_typ sap){
			// Need to reset (but don't free) head (i.e., list tweaked).
			// WARNING: THIS ROUTINE IS NOT TESTED AT ALL.
			   	assert(search); head=sap;
				search->hitlist_count=NumSeqsListGSAP(head);
			}
	void		SwapAlignment(sap_typ sap){
			   	assert(search);
				FreeHead( ); head=sap;
				search->hitlist_count=NumSeqsListGSAP(head);
				// WARNING: THIS NEEDS TO BE CHANGED TO 
				// Realign( ), which calls Fred's HMM_typ...
			}
	Int4		**ExchangeMatrix(Int4 **newMatrix){
				Int4 **old=posMatrix;
				posMatrix=newMatrix+1;		// for psi-blast
				return old;
			}
	double		**posMatchWeights( ){
			   if(posSearch && posSearch->posMatchWeights){
				return posSearch->posMatchWeights -1;  // start at 1
			   } else return 0;
			}
	double		**CopyPosPseudoFreq( );
	double		**CopyObsWtCounts( );
	double		**CopyPosMatchWeights( );
	double		*SegWeight( ){
			   if(posSearch && posSearch->posA){
				return posSearch->posA -1;	// start at 1
			   } else return 0;
			}
	Int4		RightExtent(Int4 q){
			   if(posSearch && posSearch->posExtents){
				return posSearch->posExtents[q-1].rightExtent;
			   } return 0;	// WARNING: not sure what this means!!!!!!!!!!!
			}
#endif
	Int4		*CopyOfIntervalSizes(){
			  if(posMatrix){
				Int4 *is;
				NEW(is,LenSeq(queryE)+2,Int4);
				for(Int4 i=0; i < LenSeq(queryE); i++){
					is[i+1]=posSearch->posIntervalSizes[i];
				} return is;
			  } else return 0;
			}
	Int4		**CopyOFposMatrix(){
			  if(posMatrix){
				Int4 **pm;
				NEWP(pm,LenSeq(queryE)+2,Int4);
				for(Int4 i=0; i< LenSeq(queryE); i++){
				  NEW(pm[i],nAlpha(AB)+2,Int4);
				  for(Int4 j=0; j <= nAlpha(AB); j++){
					pm[i][j] = posMatrix[i][j];
				  }
				} return pm;
			  } else return 0;
			}
	double		Lambda( ){
			   // if(compactSearch && compactSearch->search->sbp->kbp[0]->Lambda
			   if(search && search->sbp && search->sbp->kbp){
				return search->sbp->kbp[0]->Lambda;
			    } else return 0;
			}
	sap_typ 	FakeSearchCMSA(Int4 Threshold,cma_typ cma);
	sap_typ		FakeSearchGSAP(Int4 Threshold,sap_typ sap);
	sap_typ 	StealSAP( ){ sap_typ sap = head; head = 0; return sap; }
	double		*AveSeqWeights( );
	BooLean		FractionSeqsAligned(double *fract_aligned);
	void		CreateHitSet(){
				HitSet=MakeSet(NSeqsSeqSet(seqset)+2); ClearSet(HitSet);
			}
	set_typ		RtnHitSet(){ set_typ rtn = HitSet; HitSet=0; return rtn; }
	void		CreateExtendSet(){
				ExtendSet=MakeSet(NSeqsSeqSet(seqset)+2); ClearSet(ExtendSet);
			}
	set_typ		RtnExtendSet(){ set_typ rtn = ExtendSet; ExtendSet=0; return rtn; }
	void		ProvideSkipSet(set_typ skip_set){
			   SkipSet = skip_set;
			   if(SetN(SkipSet) <= NSeqsSeqSet(seqset)){
				fprintf(stderr,"SkipSet_N = %d <= NSeqsSeqSet = %d\n",
					SetN(SkipSet),NSeqsSeqSet(seqset));
				print_error("SkipSet too small for gpsi_typ");
			   }
			}
	void		CreateMarginalSet(){
				MarginalSet=MakeSet(NSeqsSeqSet(seqset)+2); ClearSet(MarginalSet);
			}
	void		CreateMarginalSet(double fraction){
				CreateMarginalSet();
				if(fraction > 1.0 || fraction < 0.01){
					print_error("MarginalGapTrigger out of range");
				}
				MarginalGapTrigger = fraction;
			}
	set_typ		RtnMarginalSet(){ set_typ rtn = MarginalSet; MarginalSet=0; return rtn; }
	void		SetGapTriggerInBits(double value){
			   gap_trigger_bits=value;
			}
private:
	Boolean		FakeGappedAlignmentWithTracebackCMSA2(gab_typ gap_align,
        			Int4 sq, e_type qE, e_type sE, Int4 motif, 
				Int4 score, cma_typ cma);	// Old routine...
	Boolean		FakeGappedAlignmentWithTracebackCMSA(gab_typ gap_align,
        			Int4 sq, e_type qE, e_type sE, Int4 motif, 
				Int4 score, cma_typ cma);
	Boolean		FakeIdenticalTracebackCMSA(gab_typ gap_align, Int4 score, e_type qE);
	sap_typ 	FakeGBlastEngineCoreCMSA(Int4 Threshold,BooLean srch_subset, cma_typ cma);
	sap_typ		FakeGBlastEngineCoreGSAP(Int4 Threshold,BooLean srch_subset,sap_typ sap);

	Int4  		GapAlignStart(unsigned char *,unsigned char *,Int4,
			Int4,Int4,Int4,Int4 **);
	sap_typ		GBlastEngineCore(Int4,BooLean,char,sap_typ);
	void		init(e_type,ss_type,double,double,double,Int4,Int4);
	void		SetUpSearch(FILE *);
	void		FreeHead( ){ FreeGSeqAlignList(head); head=0; }
	a_type		AB;	// alphabet
	set_typ		SubSet;
	set_typ		HitSet;		// returns sequence hits by seq id.
	set_typ		ExtendSet;	// returns sequences extended by seq id.
	set_typ		MarginalSet;	// returns seqs with ungapped score >= MarginalGapTrigger.
	double		MarginalGapTrigger;
	double		gap_trigger_bits;
	set_typ		SkipSet;	// sequences in this set are skipped.
	ss_type		seqset;
	Int4		dbsLength;
	e_type		queryE;
        bsb_typ		search;
        sap_typ		head;
        psim_typ	posSearch;
        csi_typ		*compactSearch;
        Boolean		recoverCheckpoint,alreadyRecovered;
	Boolean		verbose;
        Int4		thisPassNum,hitlist_size;
	Int4		pseudoCountConst,maxNumPasses;
	// Int4		number_of_cpus;
	// Int4		window; // 2 hit window size (default 40)
	// Int4		hsp_max_window; // (default 11)
	Int4		gap_open,gap_extend;
	Int4		gap_x_dropoff,gap_x_dropoff_final,gap_trigger;
	double		Ecut,ethresh,searchsp_eff;
	double		word_extend_dropoff_in_bits;  // 10.0 or 7.0 == default
	double		x_parameter,x_parameter_final;
};

#endif

