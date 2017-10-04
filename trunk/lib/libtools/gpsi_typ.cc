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

#include "gpsi_typ.h"

gpsi_type::gpsi_type(e_type qE,ss_type data,double Ethresh,double Ecutoff,
	double xdrop_final,Int4 hpsz,Int4 maxNumPass)
{
	FILE *fp=0;
	init(qE,data,Ethresh,Ecutoff,xdrop_final,hpsz,maxNumPass); SetUpSearch(fp); 
}

gpsi_type::gpsi_type(e_type qE,ss_type data,double Ethresh,double Ecutoff,
	double xdrop_final,Int4 hpsz,Int4 maxNumPass, char *chckptfile)
{
	FILE    *fp=0;
	if(chckptfile != NULL){ 
		// std::cerr << "checkin file = ";
		// std::cerr << chckptfile; std::cerr << std::endl;
		fp= open_file(chckptfile,"","r");
	}
	init(qE,data,Ethresh,Ecutoff,xdrop_final,hpsz,maxNumPass); 
	SetUpSearch(fp); if(fp) fclose(fp);
}

gpsi_type::gpsi_type(e_type qE,ss_type data,double Ethresh,double Ecutoff,
	double xdrop_final,Int4 hpsz,Int4 maxNumPass, FILE *fp,Int4 dbl)
{
	init(qE,data,Ethresh,Ecutoff,xdrop_final,hpsz,maxNumPass);
	dbsLength=dbl; SetUpSearch(fp); 
}

gpsi_type::gpsi_type(e_type qE,ss_type data,double Ethresh,double Ecutoff,
	double xdrop_final,Int4 hpsz,Int4 maxNumPass, FILE *fp)
{
	init(qE,data,Ethresh,Ecutoff,xdrop_final,hpsz,maxNumPass);
	SetUpSearch(fp); 
}

void	gpsi_type::init(e_type qE,ss_type data,double Ethresh,double Ecutoff,
	double xdrop_final,Int4 hpsz,Int4 maxNumPass)
{
	// recoverCheckpoint=FALSE;
	verbose=TRUE;
	ethresh=0.001;		// expect cutoff for inclusion in the profile.
	Ecut = 0.001;		// expect cutoff for reporting a hit.
	dbsLength=0;		// 
	seqset = data; queryE = qE;
	maxNumPasses=maxNumPass; hitlist_size=hpsz;
	HitSet=0; ExtendSet=0; SkipSet=0; MarginalSet=0;
	MarginalGapTrigger=0.90;
	gap_trigger_bits = 22.0;
	SubSet=MakeSet(NSeqsSeqSet(data)+2); ClearSet(SubSet); 
	ethresh=Ethresh; Ecut=Ecutoff;
	gap_open=11; gap_extend=1;	// default gap params for "BLOSUM62"
	word_extend_dropoff_in_bits=7.0;  // or 10.0?? for ungapped....
	x_parameter=15.0; x_parameter_final=xdrop_final;
}

void	gpsi_type::SetUpSearch(FILE *chkpt_fp)
{
        thisPassNum = 0; head=NULL; posMatrix=NULL;
	AB = SeqSetA(seqset);
        if(maxNumPasses==1) ethresh=0.0; // zero out e-value threshold if not used
	search=GBlastSearchBlkNew(hitlist_size);
	search->positionBased=FALSE;		// FALSE by default
        search->query_length=LenSeq(queryE);
        // search->query_sequence=SeqPtr(queryE)+1;
        search->query_sequence=XSeqPtr(queryE)+1;
	search->ethresh=ethresh;	
	search->pseudoCountConst=10;	
	search->maxNumPasses=maxNumPasses;	
        search->posConverged=FALSE;
	// OBTAIN GBlastKarlin PARAMETERS.
        // search->sbp = GBLAST_ScoreBlkNew(AB,SeqPtr(queryE)+1,LenSeq(queryE),gap_open,
	//	gap_extend, TotalSeqSet(seqset), NSeqsSeqSet(seqset));
#if 0
        search->sbp = GBLAST_ScoreBlkNew(AB,XSeqPtr(queryE)+1,LenSeq(queryE),gap_open,
		gap_extend, TotalSeqSet(seqset), NSeqsSeqSet(seqset));
#else
	if(dbsLength <= 0) dbsLength=TotalSeqSet(seqset);
        search->sbp = GBLAST_ScoreBlkNew(AB,XSeqPtr(queryE)+1,LenSeq(queryE),gap_open,
		gap_extend, dbsLength, NSeqsSeqSet(seqset));
#endif
        alreadyRecovered = FALSE; posSearch = NULL; compactSearch = NULL;
// fprintf(stderr,"DEBUG Z\n");
	if(chkpt_fp != NULL){ 
		recoverCheckpoint=TRUE; search->positionBased = TRUE; 
		posSearch=MakePSIM(NSeqsSeqSet(seqset), AB);
		compactSearch = compactSearchNew(search);
		posInitializeInformation(posSearch,search);
		
		// chkpt_fp= open_file(chckptfile,"","r");
		// if(!posReadCheckpoint(posSearch,compactSearch,chckptfile))
// PutSeq(stderr,queryE,AB);
		if(!posReadCheckpoint(posSearch,compactSearch,chkpt_fp))
		  print_error("gpsi_typ: Error recovering from checkpoint");
		else posMatrix = posSearch->posMatrix;
		// fclose(chkpt_fp);
		// std::cerr << chckptfile; std::cerr << std::endl;
		// outputPosMatrix(stderr,posSearch,compactSearch,AB);
// fprintf(stderr,"DEBUG ZY\n");
	} else {
// fprintf(stderr,"DEBUG ZZ\n");
		recoverCheckpoint=FALSE; search->positionBased = FALSE; 
	}
}

sap_typ	gpsi_type::TrimmedSearch(Int4 Threshold,double misalncut,char mode)
{
	Search(Threshold,FALSE,mode);
        if(posMatrix){
           sap_typ      new_sap,this_sap,nxt_sap;
           Int4         TrimWindow=30;
           double       PerNats=2/log(2);
           for(this_sap=head,nxt_sap=head->next; nxt_sap != 0;
                                this_sap=nxt_sap,nxt_sap=nxt_sap->next){
	      if(new_sap=MargProbTrimSAP(nxt_sap,posMatrix-1,PerNats,misalncut,
        				gap_open,gap_extend,TrimWindow,AB)){
                 sap_typ nxt_nxt_sap=nxt_sap->next;
                 this_sap->next = new_sap;
                 new_sap->next = nxt_sap->next;
                 nxt_sap->next=0;
                 FreeGSeqAlignList(nxt_sap);
                 nxt_sap=new_sap;
              }
           } TweakAlignment(head);
        } return head;
}

sap_typ	gpsi_type::Search(Int4 Threshold,BooLean second_pass, char mode)
{
        thisPassNum++;
	assert(!(second_pass && thisPassNum==1));
	assert(!search->posConverged);
// fprintf(stderr,"thisPassNum=%d; recoverCheckpoint=%d\n",thisPassNum,recoverCheckpoint);
        if(1==thisPassNum 
		&& (!recoverCheckpoint))
		posSearch=MakePSIM(NSeqsSeqSet(seqset), AB);
// fprintf(stderr,"DEBUG ZX\n");	/// purify complaining on some sequences
// assert(thisPassNum==24);
        if(compactSearch != NULL) compactSearchDestruct(compactSearch);
	if(recoverCheckpoint) {
                    recoverCheckpoint=FALSE; alreadyRecovered=TRUE;
        }
	// if(verbose) fprintf(stderr,"Psi-BLAST iteration %d\n",thisPassNum);
	head=GBlastEngineCore(Threshold,second_pass,mode,0);
        compactSearch = compactSearchNew(search);
	// The next two calls (after the "if") are diagnostics for Stephen. */
	// Don't perform this if only one pass will be made (i.e., standard BLAST) */
        if (1 != search->maxNumPasses) {
              if ((1 == thisPassNum)  && (!alreadyRecovered))
                   posInitializeInformation(posSearch, search);
              posPrintInformation(posSearch, search, thisPassNum);
	}
	if(1 != search->maxNumPasses) { // Have things converged? 
		  posConvergencedTest(posSearch, search, head);
        } 
	if(second_pass) { thisPassNum--; search->posConverged = FALSE; }
// fprintf(stderr,"DEBUG ZW\n");
	return head;
}

static const double gpsi_ln2 = 0.693147180559945309;

Int4  gpsi_type::GapAlignStart(unsigned char *query,unsigned char *subject,
	Int4 q_start,Int4 s_start,Int4 q_end, Int4 s_end, Int4 **matrix)
//	Function to look for the highest scoring window (of size window)
//	in an HSP and return the middle of this.  Used by the gapped-alignment
//	functions to start the gapped alignments.
{
	Int4		i,max_offset,score,max_score,hsp_end;
	unsigned char 	*q_var,*s_var;
	const short 	Window=11;
	Int4		hsp_len = q_end - q_start + 1;

	if(hsp_len <= Window) { max_offset = hsp_len/2; }
	else {
	    hsp_end = q_start + Window;
	    q_var = query + q_start; s_var = subject + s_start;
	    for(score=0,i=q_start; i<hsp_end; i++) {
		  if(!(search->positionBased)) score += matrix[*q_var][*s_var];
		  else score += posMatrix[i][*s_var];
		  q_var++; s_var++;
	    }
	    max_score = score; max_offset = hsp_end - 1;
	    hsp_end = q_end;
	    for (i=q_start + Window; i<hsp_end; i++) {
		if(!(search->positionBased)){
		    score -= matrix[*(q_var-Window)][*(s_var-Window)];
		    score += matrix[*q_var][*s_var];
		} else {
		    score -= posMatrix[i-Window][*(s_var-Window)];
		    score += posMatrix[i][*s_var];
                }
		if(score > max_score) { max_score = score; max_offset = i; }
		q_var++; s_var++;
	    }
	    max_offset -= Window/2;
	    max_offset -=q_start; 
	    // max_offset++; // AFN: Tom appears to be too low by one??
	} return max_offset;
}

BooLean	*gpsi_type::GetHitList( )
// This is currently very inefficient... Should be fixed...
{
	Int4	N=NSeqsSeqSet(seqset);
	UInt4 I;
	BooLean	*hits;
	e_type	sE;

   if(head){
	NEW(hits,N+3,BooLean);
        for(Int4 sq=1; sq<=NSeqsSeqSet(seqset); sq++){
	   I = SeqI(SeqSetE(sq,seqset));
	   for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	      dsp_typ dsp = gsap->segs;
	      sE = dsp->sE;
	      if(I == SeqI(sE)){ hits[I] = TRUE; break; }
	   }
	}
	return hits;
   } else return 0;
}

void	gpsi_type::PutMissedSeq(FILE *fp)
{
  e_type sE;
  assert(NSeqsSeqSet(seqset) >= 1);
  BooLean *hits=GetHitList( );
  for(Int4 sq=1; sq<=NSeqsSeqSet(seqset); sq++){
    if(!hits || !hits[sq]){ sE = SeqSetE(sq,seqset); PutSeq(fp,sE,AB); }
  }
  free(hits);
}

sap_typ	gpsi_type::SearchWithPrior(Int4 Threshold,double misalncut,char mode,
		cma_typ mtfcma)
{
	print_error("SearchWithPrior( ) not yet working....");
	BooLean second_pass=FALSE;
        thisPassNum++;
	assert(!(second_pass && thisPassNum==1));
	assert(!search->posConverged);
        if(1==thisPassNum && (!recoverCheckpoint))
		posSearch=MakePSIM(NSeqsSeqSet(seqset), AB);
        if(compactSearch != NULL) compactSearchDestruct(compactSearch);
	if(recoverCheckpoint) { recoverCheckpoint=FALSE; alreadyRecovered=TRUE; }
	if(mtfcma){ 
		head=FakeGBlastEngineCoreCMSA(Threshold,FALSE,mtfcma);
        	compactSearch = compactSearchNew(search);
#if 0
        if (1 != search->maxNumPasses) {
              if ((1 == thisPassNum)  && (!alreadyRecovered))
                   posInitializeInformation(posSearch, search);
              posPrintInformation(posSearch, search, thisPassNum);
	}
	if(1 != search->maxNumPasses) { // Have things converged? 
		  posConvergencedTest(posSearch, search, head);
	}
#endif
  		if(verbose) PutMultiGSeqAlign(stderr, head, 60, AB);
		// exit(1);
	} else head = 0;
	// go right into full search without reinitializing...
	head=GBlastEngineCore(Threshold,second_pass,mode,head);	// pass in head or 0.
        compactSearch = compactSearchNew(search);
	// The next two calls (after the "if") are diagnostics for Stephen. */
	// Don't perform this if only one pass will be made (i.e., standard BLAST) */
        if (1 != search->maxNumPasses) {
              if ((1 == thisPassNum)  && (!alreadyRecovered))
                   posInitializeInformation(posSearch, search);
              posPrintInformation(posSearch, search, thisPassNum);
	}
	if(1 != search->maxNumPasses) { // Have things converged? 
		  posConvergencedTest(posSearch, search, head);
        } 
	if(second_pass) { thisPassNum--; search->posConverged = FALSE; }
	// end of Search procedure...
        if(misalncut < 1.0 && posMatrix){
	   print_error("This option may work, but it is not yet tested! I'll have to kill it...");
           sap_typ      new_sap,this_sap,nxt_sap;
           Int4         TrimWindow=30;
           double       PerNats=2/log(2);
           for(this_sap=head,nxt_sap=head->next; nxt_sap != 0;
                                this_sap=nxt_sap,nxt_sap=nxt_sap->next){
	      if(new_sap=MargProbTrimSAP(nxt_sap,posMatrix-1,PerNats,misalncut,
        				gap_open,gap_extend,TrimWindow,AB)){
                 sap_typ nxt_nxt_sap=nxt_sap->next;
                 this_sap->next = new_sap;
                 new_sap->next = nxt_sap->next;
                 nxt_sap->next=0;
                 FreeGSeqAlignList(nxt_sap);
                 nxt_sap=new_sap;
              }
           } TweakAlignment(head);
        } return head;
}



sap_typ	gpsi_type::GBlastEngineCore(Int4 Threshold,BooLean srch_subset,char mode,
		sap_typ old_head)
{
  Int4		r,s,e,d,sq2,max,n_ug;
  Int4		nhits,GapCut,total_hits,word_extend_dropoff,opt_offset;
  BooLean	debug=FALSE; // debug=TRUE;
  e_type	sE;
  sap_typ	sap,sap2;
  hhp_typ	hhp;		// Temporary heap for HSPs.
  hsp_typ	hsp;
  gpb_typ	gpb;
  gb_typ	GB; 
  double	key;
  // double	gap_trigger_bits = 22.0;  // ungapped cutoff for triggering gapxdrop
  Int4		marginal_gap_tigger;
  // AFN: need to modify gap_trigger_bits when evalue is high & database is small!!!!!
  h_type	HG=0,xHG=0;
  static UInt4	num_calls=0;

  num_calls++;
  assert(NSeqsSeqSet(seqset) >= 1);
  // time_t  time1=time(NULL); 
  // change to: GABNew(queryE,gap_open,gap_extend,SMatrixSBP(search->sbp)); ???
  gab_typ gap_align= GABNew(gap_open, gap_extend);
  gap_align->include_query = 0;	// query length from q_start that MUST be aligned.
  // gap_align->query = SeqPtr(queryE) + 1;
  gap_align->query = XSeqPtr(queryE) + 1;
  gap_align->query_length = LenSeq(queryE);
  gap_align->decline_align = (-(GBLAST_SCORE_MIN));
  gap_align->matrix=SMatrixSBP(search->sbp);  // NOTE: may need to change sbp->matrix!!!

  if(posMatrix != NULL){		// set to posMatrix for PSI-BLAST
	SetPsiStatsSBP(search->sbp);	// use PSI-BLAST statistics
	word_extend_dropoff=WordXDropoffSBP(word_extend_dropoff_in_bits,search->sbp);
	gap_trigger=GapTriggerSBP(gap_trigger_bits,search->sbp);
	marginal_gap_tigger = (Int4) ceil(MarginalGapTrigger * gap_trigger);
     	gap_align->posMatrix=posMatrix; gap_align->positionBased = TRUE;	
	gpb=MakeGPsiBlst(500,gap_trigger,Threshold,queryE,AB,posMatrix-1,word_extend_dropoff);
	// subtract 1 because GPsiBlst starts at 1 while posMatrix starts at 0.
  } else { 
	SetStdStatsSBP(search->sbp);	// use standard statistics
	word_extend_dropoff=WordXDropoffSBP(word_extend_dropoff_in_bits,search->sbp);
	gap_trigger = GapTriggerSBP(gap_trigger_bits,search->sbp);
	marginal_gap_tigger = (Int4) ceil(MarginalGapTrigger * gap_trigger);
#if 1	// AFN: modification (not sure this is correct)
	//gap_trigger /= 2.0; // Fudge?? half bits???
#endif
#if 0	// AFN: modification (not sure this is correct)
	Int4 s2;
	double Ecut2;
	GBlastCutoffs_simple(&s2, &Ecut2,compactSearch->kbp_gap_std[0],
		search->sbp->searchsp_eff,TRUE);
	s2 = MAXIMUM(Int4,s2, 1);
	gap_trigger = MINIMUM(double,s2,gap_trigger);
#endif
	gap_align->posMatrix=NULL; gap_align->positionBased = FALSE; 
	GB=MkGBlast2(500,gap_trigger,Threshold,queryE,AB,word_extend_dropoff);
  }
  gap_x_dropoff_final = GapXDropoffSBP(x_parameter_final,search->sbp);
  gap_x_dropoff = GapXDropoffSBP(x_parameter,search->sbp);
  gap_align->x_parameter=gap_x_dropoff_final;
  GapCut=GappedEvalueToScoreSBP(Ecut, search->sbp);
  if(verbose){
   if(num_calls == 1){
    fprintf(stderr,"word_extend_dropoff = %d; gap_trigger = %d; GapCut = %d (E <= %g)\n",
	word_extend_dropoff,gap_trigger,GapCut,Ecut);
    fprintf(stderr,"============ gap_x_drop = %d (%d)\n",gap_x_dropoff,gap_x_dropoff_final);
   }
  }
  if (thisPassNum > 1 && old_head==0) { ResetGBlastSearchBlk(search); FreeHead( ); }
  rhp_typ rhp = MakeBRHHeap(hitlist_size);
  brh_typ results; total_hits=0;
  if(debug){
	 HG=Histogram("Blast MSP scores",0,100,2.0);
	 xHG=Histogram("Blast gapped extension scores",0,100,2.0);
  }
  for(n_ug=0,sq2=1; sq2<=NSeqsSeqSet(seqset); sq2++){
    // if(srch_subset && sq2%100 == 0) fprintf(stderr,".");
    if(srch_subset && !MemberSet(sq2,SubSet)) continue;  // subset search mode 
    if(SkipSet && !MemberSet(sq2,SkipSet)) continue;  // skip these sequences
    sE = SeqSetE(sq2,seqset);
#if 0
    if(mode=='x' ProcessSeqPSeg(17,2.2,2.5,100,sE,AB);  // AFN: TEST...
#else
    if(mode=='x' || mode=='b' || mode=='B') ProcessSeqPSeg(17,2.2,2.5,5,sE,AB);  // AFN: TEST...
    if(mode=='B') ProcessSeqPSeg(45,3.4,3.75,5,sE,AB);  // Mask out coiled coils (too much?) as well...
    if(mode=='b') ProcessSeqPSeg(45,3.3,3.65,5,sE,AB);  // Mask out coiled coils as well...
// PutXSeq(stderr,sE,AB); fprintf(stderr,"mode = %c\n",mode);
#endif
    if(posMatrix != NULL){
	max=MatcherGPsiBlstStr(LenSeq(sE),XSeqPtr(sE),gpb); nhits=nHitsGPsiBlst(gpb);
    } else {
     	max=MatcherGBlastStr(LenSeq(sE),XSeqPtr(sE),GB,AlphaR(AB));
	nhits=nHitsGBlast(GB);
    }
    if(mode=='x' || mode=='b') UnXSeq(sE);	// Seg only for initial hits...
// PutXSeq(stderr,sE,AB);
    if(debug && max > 0) IncdHist(max,HG);
#if 1	// 12_8_08: keep track of ungapped hits over 
    if(MarginalSet && max >= marginal_gap_tigger){
      AddSet(sq2,MarginalSet);	// Keep track of marginal hits.
    }
    if(ExtendSet && nhits > 0){
      AddSet(sq2,ExtendSet);	// Keep track of extensions.
    }
#endif
    // if(max==0) PutSeqSetE(stderr,sq2,seqset);  // are typically short fragments
    if(nhits > 0){
      total_hits+=nhits;
// std::cerr << "nhits = "; std::cerr << nhits; std::cerr << std::endl;
      gap_align->subject = SeqPtr(sE) + 1;
      gap_align->subject_length = LenSeq(sE);
      hhp = MkHSPHeap(nhits+2); 
      for(max=0,r=1; r<=nhits; r++){
  	if(posMatrix != NULL){
		s=StartGPsiBlst(r,gpb); e=EndGPsiBlst(r,gpb); d=DiagGPsiBlst(r,gpb);
	} else { s=StartGBlast(r,GB); e=EndGBlast(r,GB); d=DiagGBlast(r,GB); }
	// if(debug) PutDiagonalSeq(stderr,d,queryE, sE, AB);
	s--; e--;	// adjustment for ncbi starts at 0 instead of 1.
	// opt_offset=GapAlignStart(SeqPtr(queryE)+1,SeqPtr(sE)+1,
	// 			s+d,s,e+d,e,gap_align->matrix);
	opt_offset=GapAlignStart(XSeqPtr(queryE)+1,SeqPtr(sE)+1,
				s+d,s,e+d,e,gap_align->matrix);
	gap_align->q_start=s+d+opt_offset;
	gap_align->s_start=s+opt_offset;
   	// == Starting points (on original HSP) for a gapped extension with X dropoff.
#if 0	// Potentially faster route for a large database search???
  	gap_align->x_parameter=gap_x_dropoff;
	PerformGappedAlignment(gap_align);
	if(gap_align->score >= GapCut){
  	   gap_align->x_parameter=gap_x_dropoff_final;
	   PerformGappedAlignmentWithTraceback(gap_align);
	   hsp=MakeHSPTyp(gap_align,search->sbp);
	   assert(InsertHSPHeap(hsp, hhp) != NULL);
	}
#endif
#if 1
        if(!WithinHSPHeap(s+d,s,e+d,e,hhp)){  // if not already within a previous hsp.
	  PerformGappedAlignmentWithTraceback(gap_align);
	  max = MAXIMUM(Int4,max,gap_align->score);
	  if(gap_align->score >= GapCut){
	        hsp=MakeHSPTyp(gap_align,search->sbp);
		assert(InsertHSPHeap(hsp, hhp) != 0);
	  } else { GXEBDelete(gap_align->edit_block); gap_align->edit_block=0; }
        }
#endif
      }
      AddSet(sq2,SubSet);	// Keep track of potential hits.
      // if(max > gap_trigger+10) AddSet(sq2,SubSet);	// Keep track of potential hits.
      if(debug && max > 0) IncdHist(max,xHG);
// std::cerr << sq2; std::cerr << std::endl;
      Int4 n=PurgeHSPHeap(hhp); 
      if(n > 0){
	results = MakeGBLASTResultHitlist(n+1,sE);
	// double max_score=-99999;
	while((hsp=DelMinHSPHeap(&key, hhp)) != 0){
		assert(results->subject_id == sq2);
		// max_score = MAXIMUM(double, max_score,key);
 		AddHspBRH(hsp,results);
	}
        if(InsertBRHHeap(results,rhp)==NULL) results=NilBRH(results);
#if 1	// 12_8_08: keep track of ungapped hits 
	if(HitSet) { AddSet(sq2,HitSet); }
#else	// 
	if(HitSet) {
		// bin HitSet by scores.
		// or by E-values:
		// HitSet[Eval]
		AddSet(sq2,HitSet); 
	}
#endif
      } NilHSPHeap(hhp);
    }
  }
  if(debug){ 
	PutHist(stderr,60,HG); NilHist(HG);
	fprintf(stderr,"sequence in database: %d\n",NSeqsSeqSet(seqset));
	PutHist(stderr,60,xHG); NilHist(xHG);
  }
  head=old_head; 
  while((results=DelMinBRHHeap(&key,rhp)) != 0){
	sap = ExtractAlnBRH(results,queryE);
        if(head == 0) { head = sap; }
        else{ for(sap2=head;sap2->next;) sap2=sap2->next; sap2->next=sap; }
  	Int4 seqIndex=search->hitlist_count;
	assert(search->hitlist_max > search->hitlist_count);
	search->hitlist_count++;
	search->results[seqIndex]=results;
  } NilBRHHeap(rhp); GABDelete(gap_align);
  // if(verbose) fprintf(stderr,"total_hits = %d\n",total_hits);
  if(verbose){
     if(search->hitlist_count > 0){
	PutHitList(stdout,60);
	fprintf(stdout,"\n"); 
     } else { 
	fprintf(stderr,"Query: "); 
	PutSeqID(stderr,queryE); 
	fprintf(stderr," (no hits)\n"); 
     }
  }
  if(posMatrix != NULL) NilGPsiBlst(gpb); else NilGBlast(GB);
  if(!srch_subset){
	// if(verbose && num_calls == 1) 
	if(verbose) 
	 fprintf(stderr,"%d sequences extended; %d detected\n",
		CardSet(SubSet),search->hitlist_count);
  } else fprintf(stderr,"second pass: %d sequences detected\n",
		search->hitlist_count);
#if 0
  if(verbose){
	fprintf(stderr,"time = %.1f seconds\n",difftime(time(NULL),time1));
  }
#endif
  return head;
}

void	gpsi_type::AppendFakeResults(gpsi_type *fake_gpsi)
// WARNING: This assumes that a search has been completed...
{
  // fake_gpsi->FakeSearchCMSA(T,mtfcma);

  sap_typ sap=fake_gpsi->StealSAP();
  sap_typ sap2;

  if(head == 0) { head = sap; }
  else{ for(sap2=head;sap2->next;) sap2=sap2->next; sap2->next=sap; }

  Int4 seqIndex=this->search->hitlist_count;
  Int4 fakeseqIndex=fake_gpsi->search->hitlist_count;
  this->search->hitlist_count += fake_gpsi->search->hitlist_count;
  assert(search->hitlist_max > this->search->hitlist_count);
  for(Int4 r=0; r < fakeseqIndex;r++){
	this->search->results[seqIndex]=fake_gpsi->search->results[r];
	seqIndex++;
  } assert(seqIndex == this->search->hitlist_count);
}

#if 1	// NEW...
void	gpsi_type::PutDistHits(FILE *fp,Int4 len_line)
{
  Int4    i,maxscore=INT4_MIN,minscore=INT4_MAX,score;
  if(search==NULL || search->results == NULL) return;
  for(i=0; i< search->hitlist_count; i++){
    brh_typ results = search->results[i];
    if(results){
	// Int4 bit_score = (Int4) GappedBitScoreSBP(results->high_score, search->sbp);
	Int4 bit_score = (Int4) -log10(results->best_evalue);
	if(bit_score > maxscore) maxscore=bit_score;
	if(bit_score < minscore) minscore=bit_score;
    }
  }
  h_type HG=Histogram("-log10(E-value)",minscore,maxscore,1.0);
  for(i=0; i< search->hitlist_count; i++){
    brh_typ results = search->results[i];
    if(results){
	score = (Int4) -log10(results->best_evalue);
	for(Int4 s=score; s >= minscore; s--){
		IncdHist(s,HG);
	}
	// IncdHist(-log10(results->best_evalue),HG);
   }
  } PutHist(fp,len_line,HG); NilHist(HG);
}
#endif	// END NEW...
#if 0
void	gpsi_type::PutDistHits(FILE *fp,Int4 len_line)
{
  if(search==NULL || search->results == NULL) return;
  // h_type HG=Histogram("HSP bit scores",0,2000,1.0);
  h_type HG2=Histogram("-log10(E-value)",-100,1000,0.25);
  for(Int4 x,i=0; i< search->hitlist_count; i++){
    brh_typ results = search->results[i];
    if(results){
	double bit_score = GappedBitScoreSBP(results->high_score, search->sbp);
	// IncdHist(bit_score,HG);
	IncdHist(-log10(results->best_evalue),HG2);
   }
  }
  // PutHist(fp,len_line,HG); NilHist(HG);
  PutHist(fp,len_line,HG2); NilHist(HG2);
}
#endif

void	gpsi_type::PutHitList(FILE *fp,Int4 len_line)
{
  static UInt4 num_calls=0;

  num_calls++;
  if(search==NULL || search->results == NULL) return;
  char str[200];
  // fprintf(fp, "\nResults from round %d\n", thisPassNum);
  fprintf(fp,"\nQuery: ");
  PutSeqID(fp,queryE); fprintf(fp,"\n");
  // PutSeqInfo(fp,queryE); fprintf(fp,"\n");
  sprintf(str,"%c%ds  Score    E%c",'%',len_line,'\n');
  fprintf(fp,str,"   ");
  sprintf(str,"%c-%ds  (bits) Value%c%c",'%',len_line,'\n','\n');
  fprintf(fp,str,"Sequences producing significant alignments:");
  for(Int4 x,i=0; i< search->hitlist_count; i++){
    brh_typ results = search->results[i];
    if(results){
	StrSeqInfo(str,results->sE);
	for(x=0; x < len_line; x++){
		if(str[x] == 0) break;
		else fprintf(fp,"%c",str[x]);
	}
	while(x < len_line) { fprintf(fp,"."); x++; }
	double bit_score = GappedBitScoreSBP(results->high_score, search->sbp);
	fprintf(fp,"  %5.0f %1.0lg\n",bit_score,results->best_evalue);
   }
  } if(num_calls == 1) PutKarlinAltschulSBP(stderr,search->sbp);
}

gpsi_type::~gpsi_type( )
{
	if(posSearch && posSearch->posMatrix) posCleanup(posSearch,compactSearch);
	NilPSIM(posSearch); posSearch = 0; posMatrix=0;
	if(compactSearch != NULL) compactSearchDestruct(compactSearch);
	if(search!=NULL){
		GBLAST_ScoreBlkDestruct(search->sbp); GBlastSearchBlkDelete(search);
	} FreeHead( ); NilSet(SubSet); 
	if(HitSet) NilSet(HitSet);
	if(ExtendSet) NilSet(ExtendSet);
}

//********************* NEW: Marginal Probability PSI-BLAST ************************
BooLean gpsi_type::ComputeMatrix(char *checkpoint, double ***MatrixMP)
{
        search->positionBased = TRUE;
        if(alreadyRecovered) {
                  posCheckpointFreeMemory(posSearch, compactSearch->qlength);
		  posMatrix=0; alreadyRecovered = FALSE;
        }
        if(posMatrix) { posCleanup(posSearch, compactSearch); posMatrix=0; }
	if(verbose && search->maxNumPasses > 2) {
	    fprintf(stderr,"thisPassNum = %d; posConverged=%d; maxNumPasses = %d\n",
 		thisPassNum,search->posConverged,search->maxNumPasses);
	}
        if (!search->posConverged && 
	   (search->maxNumPasses == 0 || (thisPassNum < search->maxNumPasses))){
		// fprintf(stderr,"taking checkpoint...(%s)\n",checkpoint);
            posMatrix=CposComputation(posSearch,search,compactSearch,head,
			checkpoint,MatrixMP);
        } return TRUE;
}
//********************* END: Marginal Probability PSI-BLAST ************************

BooLean gpsi_type::ComputeMatrix(char *checkpoint)
{
	BooLean rtn=0;
	FILE	*chkfp=0;
	if(checkpoint && !search->posConverged && 
	   (search->maxNumPasses == 0 || (thisPassNum < search->maxNumPasses))){
	    fprintf(stderr,"creating profile...(%s)\n",checkpoint);
	    // fprintf(stderr,"taking checkpoint...(%s)\n",checkpoint);
	    chkfp=open_file(checkpoint,"","w");
	}
	rtn=ComputeMatrix(chkfp,checkpoint);
	if(chkfp) fclose(chkfp);
	return rtn;
}

BooLean gpsi_type::ComputeMatrix(FILE *chkpt_fp,char *checkpoint)
{
        search->positionBased = TRUE;
        if(alreadyRecovered) {
                  posCheckpointFreeMemory(posSearch, compactSearch->qlength);
		  posMatrix=0; alreadyRecovered = FALSE;
        }
        if(posMatrix) { posCleanup(posSearch, compactSearch); posMatrix=0; }
	if(verbose && search->maxNumPasses > 2){
	    fprintf(stderr,"thisPassNum = %d; posConverged=%d; maxNumPasses = %d\n",
 		thisPassNum,search->posConverged,search->maxNumPasses);
	}
        if ((!search->posConverged && 
	   (search->maxNumPasses == 0 || (thisPassNum < search->maxNumPasses)))
#if 0	// AFN: added but don't know if this is needed... if run j==2 should do the same?
	   || (search->maxNumPasses == 1 && checkpoint)
#endif
	   			){
	    // if(checkpoint && verbose) fprintf(stderr,"taking checkpoint...(%s)\n",checkpoint);
	    if(chkpt_fp && verbose) fprintf(stderr,"taking checkpoint...(%s)\n",checkpoint);
            posMatrix=CposComputation(posSearch,search,compactSearch,head,chkpt_fp);
#if 0
  Int4 querySize,c,n;
  querySize = search->query_length;

  h_type HG=Histogram("Residue counts by position",0,querySize,1.0);
  for(c = 0; c < querySize; c++) {
        n = posSearch->posCount[c];
        if(n > 0) IncdMHist((double)c+1,n ,HG);
  } PutHist(stderr,60,HG); NilHist(HG);
#endif
	    return TRUE;
        } else {
#if 0
	    fprintf(stderr,"search->posConverged = %d; search->maxNumPasses = %d; thisPassNum = %d\n",
			search->posConverged,search->maxNumPasses, thisPassNum);
#endif
	    return FALSE;
	}
}

BooLean	gpsi_type::NotConverged( )
{
#if 0
	if (!search->posConverged && (0 == search->maxNumPasses 
			|| thisPassNum < search->maxNumPasses)){
	    // outputPosMatrix(stderr,posSearch, compactSearch,AB); 
	} // a diagnostic, possibly make an option. 
#endif
	if((0 == search->maxNumPasses || thisPassNum < (search->maxNumPasses))
                && (! (search->posConverged))) return TRUE; else return FALSE; 
}

void	gpsi_type::SeeMatrix(FILE *fp)
{
     if(posMatrix){
#if 0
        smx_typ smx=GPSI2SMatrix(LenSeq(queryE),posMatrix,AB);
        PutSMatrix(fp,smx); NilSMatrix(smx);
#endif
	outputPosMatrix(fp,posSearch,compactSearch,AB);
     } else fprintf(fp,"posMatrix does not exist!\n");
}

void	gpsi_type::Test(char mode)
{
   if (head) {
	Int4	querySize=LenSeq(queryE);
	double	inc=2.0;
	h_type	iHG=0,dHG=0,mHG=0;
	switch(mode){
	  case 'A':
  	    iHG=Histogram("Inserts by position",0,querySize,inc);
  	    dHG=Histogram("Deletions by position",0,querySize,inc);
  	    mHG=Histogram("Matches by position",0,querySize,inc);
	    break;
	  case 'M':
  	    mHG=Histogram("Matches by position",0,querySize,inc);
	    break;
	  default: print_error("gpsi_type::Test( ) input error");
	    break;
	}
	Int4	*numins,*numdel,*numaln;
	NEW(numins,querySize+1, Int4); NEW(numdel,querySize+1, Int4);
	NEW(numaln,querySize+1, Int4);
	Int4	qs,ss,len,q,s;
	for(sap_typ gsap=head; gsap!=NULL; gsap=gsap->next){
	   dsp_typ dsp = gsap->segs;
	   q = dsp->starts[0]+1; s = dsp->starts[1]+1;
	   for(Int4 i=0; i < dsp->numseg; i++){
		qs = dsp->starts[2*i]; ss = dsp->starts[2*i+1];
		len = dsp->lens[i];
		if(qs == -1) numins[q]+=len;  // insertion relative to query
		// if(qs == -1) numins[q]++;  // insertion relative to query
		else if(ss == -1){
			for(Int4 j=0; j < len; j++){ numdel[q]++; q++; }
		} else { for(Int4 j=0; j < len; j++){ numaln[q]++; q++; } }
	   }
	}
 	for(Int4 c=1; c <= querySize; c++){
          if(iHG && numins[c] > 0) IncdMHist((double)c,numins[c],iHG);
          if(dHG && numdel[c] > 0) IncdMHist((double)c,numdel[c],dHG);
          if(mHG && numaln[c] > 0) IncdMHist((double)c,numaln[c],mHG);
  	} 
	if(iHG){ PutHist(stderr,60,iHG); NilHist(iHG); }
	if(dHG){ PutHist(stderr,60,dHG); NilHist(dHG); }
	if(mHG){ PutHist(stderr,60,mHG); NilHist(mHG); }
	free(numins); free(numdel); free(numaln);
   }
}

Int4	gpsi_type::DelimitSeqAlign(Int4 **SS,Int4 **SP,Int4 **ES,Int4 **EP)
{ return DelimitGSeqAlignList(head,SS,SP,ES,EP); }

void	gpsi_type::ShowAlign(FILE *fp,Int4 width, Int4 option)
{ ShowAlign(fp,width, option, head); }

void	gpsi_type::ShowAlign(FILE *fp,Int4 width, Int4 option, sap_typ sap)
{
      if(sap) {
	if(option == 4) PutMultiGSeqAlign(fp, sap, width, AB);
	else if(option == 1) PutGSeqAlignList(fp, sap, width, AB);
	else if(option == 0) ; // don't do anything...
      } else { fprintf(fp, "\n\n ***** No hits found ******\n\n"); }
}

double	*gpsi_type::FastSearch(BooLean top_only,Int4 end_sq,Int4 Threshold,
	char mode)
{ return FastSearch(top_only,end_sq,Threshold,mode,0); }

double	*gpsi_type::FastSearch(BooLean top_only,Int4 end_sq,Int4 Threshold,
	char mode,Int4 *Alnlen, Int4 IgnoreThis)
// faster version for clustering and purging sequences...
{
  Int4		r,s,e,d,time1,sq2,n,max,n_ug;
  Int4		nhits,total_hits;
  Int4		word_extend_dropoff,opt_offset;
  e_type	sE;
  gb_typ	GB; 
  // double	gap_trigger_bits = 22.0;
  double	*scores;

  thisPassNum++; assert(thisPassNum==1);
  assert(NSeqsSeqSet(seqset) >= 1);
  // time1=time(NULL); 
  gab_typ gap_align= GABNew(gap_open, gap_extend);
  gap_align->include_query = 0;	// query length from q_start that MUST be aligned.
  // gap_align->query = SeqPtr(queryE) + 1;
  gap_align->query = XSeqPtr(queryE) + 1;
  gap_align->query_length = LenSeq(queryE);
  gap_align->decline_align = (-(GBLAST_SCORE_MIN));
  gap_align->matrix=SMatrixSBP(search->sbp);  // NOTE: may need to change sbp->matrix!!!

  word_extend_dropoff=WordXDropoffSBP(word_extend_dropoff_in_bits,search->sbp);
  gap_trigger=GapTriggerSBP(gap_trigger_bits,search->sbp);
  gap_align->posMatrix=NULL; gap_align->positionBased = FALSE; 
  GB=MkGBlast2(500,gap_trigger,Threshold,queryE,AB,word_extend_dropoff);
  gap_x_dropoff_final = GapXDropoffSBP(x_parameter_final,search->sbp);
  gap_align->x_parameter=gap_x_dropoff_final;
  total_hits=0;
  NEW(scores,NSeqsSeqSet(seqset)+2,double);
  end_sq = MINIMUM(Int4,NSeqsSeqSet(seqset),end_sq);
  for(n=n_ug=0,sq2=1; sq2<=end_sq; sq2++){
    if(sq2 == IgnoreThis){ scores[sq2]=-1; continue; }
    sE = SeqSetE(sq2,seqset);
    max=MatcherGBlastStr(LenSeq(sE),XSeqPtr(sE),GB,AlphaR(AB));
    nhits=nHitsGBlast(GB);
    if(nhits > 0){
      total_hits+=nhits;
      gap_align->subject = SeqPtr(sE) + 1;
      gap_align->subject_length = LenSeq(sE);
      if(mode=='I'){	// count percent identities...
	Int4 alnlen;
	s=StartGBlast(1,GB); e=EndGBlast(1,GB); d=DiagGBlast(1,GB); 
	s--; e--;	// adjustment for ncbi starts at 0 instead of 1.
	// opt_offset=GapAlignStart(SeqPtr(queryE)+1,SeqPtr(sE)+1,
	//			s+d,s,e+d,e,gap_align->matrix);
	opt_offset=GapAlignStart(XSeqPtr(queryE)+1,XSeqPtr(sE)+1,
				s+d,s,e+d,e,gap_align->matrix);
	gap_align->q_start=s+d+opt_offset;
	gap_align->s_start=s+opt_offset;
   	// == Starting points (on original HSP) for a gapped extension with X dropoff.
	PerformGappedAlignmentWithTraceback(gap_align);
	sap_typ	sap = GXEBToGSeqAlign(gap_align->edit_block,sE,queryE);
	Int4 ident,min_leng,inserts;
	min_leng = IdentitiesGSAP(&alnlen,&ident,&inserts,sap); 
	scores[sq2] = 100.0*(double) ident/(double) min_leng;
	if(Alnlen) Alnlen[sq2] = alnlen;
	// scores[sq2] = 100.0*PercentIdentGSAP(&alnlen,sap);
	GSeqAlignFree(sap); GXEBDelete(gap_align->edit_block); 
	gap_align->edit_block=0;
      } else {
        for(max=0,r=1; r<=nhits; r++){
	  s=StartGBlast(r,GB); e=EndGBlast(r,GB); d=DiagGBlast(r,GB); 
	  s--; e--;	// adjustment for ncbi starts at 0 instead of 1.
	  // opt_offset=GapAlignStart(SeqPtr(queryE)+1,SeqPtr(sE)+1,
	// 			s+d,s,e+d,e,gap_align->matrix);
	  opt_offset=GapAlignStart(XSeqPtr(queryE)+1,XSeqPtr(sE)+1,
				s+d,s,e+d,e,gap_align->matrix);
	  gap_align->q_start=s+d+opt_offset;
	  gap_align->s_start=s+opt_offset;
   	  // == Starting points (on original HSP) for a gapped extension with X dropoff.
	  PerformGappedAlignment(gap_align);
#if 1
	  if(gap_align->score > max) { max = gap_align->score; }
#else	// in bit scores.  uses half bits.
	  // double       PerNats=2/log(2);
	  // Int4 bit_score=(Int4) floor((1.442695*(gap_align->score/PerNat))+0.5);
	  Int4 bit_score=(Int4) floor(((double)gap_align->score/(double)2)+0.5);
	  if(bit_score > max) { max = bit_score; }
#endif
	  if(top_only) break;	// just take highest hit...
        } 
	if(mode == 'S') scores[sq2] = (double) max; 
	else if(mode=='E') scores[sq2]=-log10(GappedScoreToEvalueSBP(max,search->sbp));
	else print_error("FastSearch( ) input error");
      }
    }
  } GABDelete(gap_align); NilGBlast(GB);
  // fprintf(stderr,"total_hits = %d; time = %d\n",total_hits,(clock( )-time1)/10000);
  return scores;
}

//=========================== Sampling routines ===========================

sap_typ	gpsi_type::SampleMatrix(char *checkout)
{
  Int4 **posMtrx = posMatrix; 
  sap_typ sap=SampleAlign(posMtrx,head);
  FreeGSeqAlignList(head); head=sap;  // set head to new alignment...
  ComputeMatrix(checkout);	// compute matrix with new alignment (Calls posCleanup( )).
  return head;
}

sap_typ	gpsi_type::SampleAlign(Int4 **posMtrx, sap_typ sap)
// Sample a realignment of the seq align 'sap'
// converts posMatrix into smx for now...
// WARNING: this is untested and may have lots of bugs!!
{
   Int4		numopers,score;
   char         *operation;
   sap_typ	sap2=NULL,rsap,saplist=NULL,sap_tmp;
   e_type	qE,sE;
   dsp_typ	segs;
   smx_typ	M[2];

  assert(posMtrx);
  M[1] = GPSI2SMatrix(LenSeq(queryE),posMtrx,AB);
  // PutSMatrix(stdout, M[1]);
  for(sap2=sap; sap2 != NULL; sap2=sap2->next){
	segs = sap2->segs;
	qE=segs->qE; sE=segs->sE;
	operation=SampleOperationSeqAlnSMatrixSW(10, 1, LenSeq(sE),XSeqPtr(sE),
        		0, 1, M, NULL, &numopers, &score);
   // fprintf(stderr,"score = %d\n",score);
	rsap=ToGSeqAlign(numopers, operation, qE, sE, 0, 0);
	if(saplist == NULL) { saplist=rsap; saplist->next=NULL; }
	else{ 
		for(sap_tmp=saplist;sap_tmp->next;) sap_tmp=sap_tmp->next; 
		sap_tmp->next=rsap; 
	} // printf("%s\n",operation);
	free(operation);
  } NilSMatrix(M[1]); 
  return saplist;
}

void	gpsi_type::PutRasMol(FILE *fp, sap_typ sap)
{
   Int4 ncolors=12;
   const char	*Colors[12]={"magenta","magenta","red","red",
                "yellow","yellow","green","green","cyan","cyan","blue","blue"};
   PutRasMolGSAP(fp, sap, ncolors, Colors);
}

//********************* NEW FAKE CMA SEARCH ROUTINES **********************

sap_typ	gpsi_type::FakeSearchCMSA(Int4 Threshold,cma_typ cma)
// Trick BLAST into thinking it has done a search by returning an alignment 
// derived from a cma file.
// 
// Need to initialize gpsi object using the query and a database that includes 
// the aligned sequences (put in front of database???)
// create database using MergeSeqSets( ) function.
{
        thisPassNum++;
	assert(thisPassNum==1);
	assert(!search->posConverged);
        if(1==thisPassNum && (!recoverCheckpoint))
		posSearch=MakePSIM(NSeqsSeqSet(seqset), AB);
        if(compactSearch != NULL) compactSearchDestruct(compactSearch);
	if(recoverCheckpoint) {
                    recoverCheckpoint=FALSE; alreadyRecovered=TRUE;
        }
	head=FakeGBlastEngineCoreCMSA(Threshold,FALSE,cma);
        compactSearch = compactSearchNew(search);
	// The next two calls (after the "if") are diagnostics for Stephen. */
	// Don't perform this if only one pass will be made (i.e., standard BLAST) */
        if (1 != search->maxNumPasses) {
              if ((1 == thisPassNum)  && (!alreadyRecovered))
                   posInitializeInformation(posSearch, search);
              posPrintInformation(posSearch, search, thisPassNum);
	}
	if(1 != search->maxNumPasses) { // Have things converged? 
		  posConvergencedTest(posSearch, search, head);
        } 
	return head;
}

Boolean	gpsi_type::FakeGappedAlignmentWithTracebackCMSA(gab_typ gap_align,
	Int4 sq, e_type qE, e_type sE, Int4 motif, Int4 score, cma_typ cma)
{
	assert(motif > 0 && motif <= nBlksCMSA(cma));
	Int4	i,j,s,pos[4],len=LengthCMSA(motif,cma);
	gxeb_typ edit_block;
	gss_typ *gss=gssCMSA(cma);

        Int4	buffer_size = LenSeq(qE) + LenSeq(sE);
	buffer_size+=gss->NumIns(sq)+gss->NumDel(sq);
	Int4	*tback;
	NEW(tback,buffer_size+8,Int4);

	assert(PosSiteCMSA(motif,1,pos,cma));  // position of motif in first seq...
	gap_align->query_start=pos[1]-1;	// need to adjust to PSI-BLAST zero start...
	// Look for insertions and deletions in cma file for 
	// query sequence to determine the query stop position!
	gap_align->query_stop=pos[1]+len-2;	// Can't assume no gaps!!

	assert(PosSiteCMSA(motif, sq, pos, cma)); 
	gap_align->subject_stop=gss->TrueSite(sq,pos[1]+len-2);	// from zero; not 1.
	// gap_align->subject_stop=gss->TrueSite(sq,pos[1]+len-1);

	Int4	del,ins;

// NEED TO TRIM DELETIONS AND INSERTIONS FROM ENDS...
// NOTE: gap_align starts sequences at 0 to LenSeq-1
// 1. Find start sites...
#if 0	// This was in error.
	for(i=0,s=1; IsDeletedCMSA(sq,s,cma); ) { s++; }
	gap_align->query_start += s-1;
#else	// this fixes the problem (afn: 12/17/09).
	for(i=0,s=1,j=pos[1]; IsDeletedCMSA(sq,j,cma); ) { j++; s++; }
	gap_align->query_start += s-1;
#endif
// fprintf(stderr,"%d: start query =%d; ndel=%d; pos=%d\n",sq,gap_align->query_start,s,pos[1]);
	gap_align->subject_start=(gss->TrueSite(sq,pos[1]+s-1)-1);
	// gap_align->subject_start=(gss->TrueSite(sq,pos[1]+s-2));
// fprintf(stderr,"  start subject =%d\n",gap_align->subject_start);

// 2. Find matching regions...
	for( ; s <= LengthCMSA(motif,cma); s++){
#if 0	//*********************** check for IIIddd (afn 2-20-08) ****************************
		if(s < LengthCMSA(motif,cma) && InsertionCMSA(sq,pos[1]+s-1,cma) && 
				IsDeletedCMSA(sq,pos[1]+s,cma)){
				PutSeqInfo(stderr,sE);
				fprintf(stderr,"WARNING: Insertion-to-Deletion transition.\n");
				fprintf(stderr,"sq=%d; site=%d.\nSkipped this seq.\n",
					sq,pos[1]+s-1);
				return FALSE;
		}
#endif	//***********************************************************************************
		for(del=0; IsDeletedCMSA(sq,pos[1]+s-1,cma); ){
			del++; 
			if(s >= LengthCMSA(motif,cma)){ s++; break; }
#if 1	// Is this a necessary requirement??? 
	// Yes: core dump if ignored in my_posit.cc at free(posSearch->pde_typMatrix[i]);
			// assert(!InsertionCMSA(sq,pos[1]+s-1,cma));
			if(InsertionCMSA(sq,pos[1]+s-1,cma)){
				PutSeqInfo(stderr,sE);
				fprintf(stderr,"WARNING: Deletion-to-insertion");
				fprintf(stderr,"sq=%d; site=%d.\nSkipped this seq.\n", sq,pos[1]+s-1);
				return FALSE;
			}
#endif
			s++;
		}
		if(s > LengthCMSA(motif,cma)) break;
		// at this point !IsDeletedCMSA(sq,s,cma).
		if(del){ tback[i] = -del; i++; }
		// currently at a match position...
		tback[i] = 0; i++;
		gap_align->query_stop=s-1;
		gap_align->subject_stop=gss->TrueSite(sq,pos[1]+s-1)-1;
		if(s >= LengthCMSA(motif,cma)) break;
		ins=InsertionCMSA(sq,pos[1]+s-1,cma);
		if(ins){ tback[i] = ins; i++; }
		// else { tback[i] = 0; i++; }
	} tback[i]=AlnEndFlagGXEB();
	// 0 = match; -k = k deletions; +k = k insertions...
#if 1
	gap_align->subject_start -= OffSetSeq(sE);
	gap_align->subject_stop -= OffSetSeq(sE);
#endif
	edit_block=TracebackToGXEB(0,0,
                gap_align->query_stop - gap_align->query_start+1,
                gap_align->subject_stop-gap_align->subject_start+1,
                tback, gap_align->query_start, gap_align->subject_start);
        gap_align->score = score; free(tback);
	gap_align->edit_block=edit_block;
        gap_align->edit_block->length1 = gap_align->query_length;
        gap_align->edit_block->length2 = gap_align->subject_length;
#if 0	// DEBUG...
   if(sq == 5 || sq == 6){
	gsq_typ *gsq=gss->GetGSQ(sq);
	gsq->Put(stderr, AB);
	PutGapAlign(stderr,gap_align);	// DEBUG...
	// gss->Put(stderr);
   }
#endif
	return TRUE;
}

Boolean	gpsi_type::FakeGappedAlignmentWithTracebackCMSA2(gab_typ gap_align,
	Int4 sq, e_type qE, e_type sE, Int4 motif, Int4 score, cma_typ cma)
{
	assert(motif > 0 && motif <= nBlksCMSA(cma));
	Int4	i,s,pos[4],len=LengthCMSA(motif,cma);
	gxeb_typ edit_block;
	gss_typ *gss=gssCMSA(cma);

// if(sq==15 || 1) {gss->Put(stderr,sq);  }
        Int4	buffer_size = LenSeq(qE) + LenSeq(sE);
	buffer_size+=gss->NumIns(sq)+gss->NumDel(sq);
	Int4	*tback;
	NEW(tback,buffer_size+1,Int4);

#if 0
	gap_align->query_start=1;
	gap_align->query_stop=len;
#else
	assert(PosSiteCMSA(motif,1,pos,cma));  // position of motif in first seq...
	gap_align->query_start=pos[1]-1;	// need to adjust to PSI-BLAST zero start...
	// Look for insertions and deletions in cma file for 
	// query sequence to determine the query stop position!
	gap_align->query_stop=pos[1]+len-2;	// Can't assume no gaps!!
#endif
	assert(PosSiteCMSA(motif, sq, pos, cma)); 
	gap_align->subject_stop=gss->TrueSite(sq,pos[1]+len-2);	// from zero; not 1.
	// gap_align->subject_stop=gss->TrueSite(sq,pos[1]+len-1);

	Int4	del,ins;
// NEED TO TRIM DELETIONS AND INSERTIONS FROM ENDS...
// NOTE: gap_align starts sequences at 0 to LenSeq-1
// 1. Find start sites...
	for(i=0,s=1; IsDeletedCMSA(sq,s,cma); ) { s++; }
	gap_align->query_start += s-1;
#if 0
	gap_align->subject_start=(gss->TrueSite(sq,pos[1]+s-2));
#else
	gap_align->subject_start=(gss->TrueSite(sq,pos[1]+s-1)-1);
#endif
// fprintf(stderr,"start subject =%d\n",gap_align->subject_start);
// 2. Find matching regions...
	for( ; s <= LengthCMSA(motif,cma); s++){
#if 1	//*********************** check for IIIddd (afn 2-20-08) ****************************
		if(s < LengthCMSA(motif,cma) && InsertionCMSA(sq,pos[1]+s-1,cma) && 
				IsDeletedCMSA(sq,pos[1]+s,cma)){
				PutSeqInfo(stderr,sE);
				fprintf(stderr,"WARNING: Insertion-to-Deletion transition.\n");
				fprintf(stderr,"sq=%d; site=%d.\nSkipped this seq.\n",
					sq,pos[1]+s-1);
				// return FALSE;
		}
#endif	//***********************************************************************************
		for(del=0; IsDeletedCMSA(sq,pos[1]+s-1,cma); ){
			del++; 
			if(s >= LengthCMSA(motif,cma)){ s++; break; }
#if 1	// Is this a necessary requirement??? 
	// Yes: core dump if ignored in my_posit.cc at free(posSearch->pde_typMatrix[i]);
			// assert(!InsertionCMSA(sq,pos[1]+s-1,cma));
			if(InsertionCMSA(sq,pos[1]+s-1,cma)){
				PutSeqInfo(stderr,sE);
				fprintf(stderr,"WARNING: Deletion-to-insertion");
				fprintf(stderr,"sq=%d; site=%d.\nSkipped this seq.\n", sq,pos[1]+s-1);
				return FALSE;
			}
#endif
			s++;
		}
		if(s > LengthCMSA(motif,cma)) break;
		// at this point !IsDeletedCMSA(sq,s,cma).
#if 1
if(del) fprintf(stderr,"del(%d): s=%d; i=%d; sE=-; qE=%c; tb=%d\n",sq,s,i,
		AlphaChar(ResSeq(s,qE),AB),-del); 
#endif
		if(del){ tback[i] = -del; i++; }
		// currently at a match position...
#if 0
fprintf(stderr,"s=%d; i=%d; sE=%c; qE=%c; tb=%d\n",s,i,
		AlphaChar(ResSeq(gss->TrueSite(sq,pos[1]+s-1)-OffSetSeq(sE),sE),AB),
		AlphaChar(ResSeq(s,qE),AB),0); 
#endif
		tback[i] = 0; i++;
		gap_align->query_stop=s-1;
		gap_align->subject_stop=gss->TrueSite(sq,pos[1]+s-1)-1;
		if(s >= LengthCMSA(motif,cma)) break;
		ins=InsertionCMSA(sq,pos[1]+s-1,cma);
#if 1
if(ins) fprintf(stderr,"ins(%d): s=%d; i=%d; sE=%c; qE=%c; tb=%d\n",sq,s,i,
		AlphaChar(ResSeq(gss->TrueSite(sq,pos[1]+s-1)-OffSetSeq(sE),sE),AB),
		AlphaChar(ResSeq(s,qE),AB),ins); 
#endif
		if(ins){ tback[i] = ins; i++; }
		// else { tback[i] = 0; i++; }
	} tback[i]=AlnEndFlagGXEB();
	// 0 = match; -k = k deletions; +k = k insertions...
#if 1
	gap_align->subject_start -= OffSetSeq(sE);
	gap_align->subject_stop -= OffSetSeq(sE);
#endif
	edit_block=TracebackToGXEB(0,0,
                gap_align->query_stop - gap_align->query_start+1,
                gap_align->subject_stop-gap_align->subject_start+1,
                tback, gap_align->query_start, gap_align->subject_start);
        // gap_align->score = (2*NumSeqsCMSA(cma))-sq; free(tback);
        gap_align->score = score; free(tback);
	gap_align->edit_block=edit_block;
        gap_align->edit_block->length1 = gap_align->query_length;
        gap_align->edit_block->length2 = gap_align->subject_length;
#if 1	// DEBUG...
   if(sq == 5 || sq == 6){
	gsq_typ *gsq=gss->GetGSQ(sq);
	gsq->Put(stderr, AB);
	// gsq=gss->GetGSQ(6);
	// gsq->Put(stderr, AB);
	PutGapAlign(stderr,gap_align);	// DEBUG...
	// gss->Put(stderr);
   }
#endif
	return TRUE;
}

Boolean	gpsi_type::FakeIdenticalTracebackCMSA(gab_typ gap_align, Int4 score, e_type qE)
// align the query against itself...
{
	Int4	i,s,*tback;
	gxeb_typ edit_block;
	Int4    buffer_size = LenSeq(qE) + LenSeq(qE);

	NEW(tback,buffer_size+1,Int4);
	//gap_align->subject_start=gap_align->query_start=1;
	// gap_align->query_stop=gap_align->subject_stop=LenSeq(qE);
	gap_align->subject_start=gap_align->query_start=0;
	gap_align->query_stop=gap_align->subject_stop=LenSeq(qE)-1;
	for(i=0,s=1; s <= LenSeq(qE); s++){
		tback[i] = 0; i++;
	} tback[i]=AlnEndFlagGXEB();
	// 0 = match; -k = k deletions; +k = k insertions...
	// gap_align->subject_start -= OffSetSeq(qE);
	// gap_align->subject_stop -= OffSetSeq(qE);
	edit_block=TracebackToGXEB(0,0,
                gap_align->query_stop - gap_align->query_start+1,
                gap_align->subject_stop-gap_align->subject_start+1,
                tback, gap_align->query_start, gap_align->subject_start);
        gap_align->score = score; // Not sure that this will work...
	free(tback);
	gap_align->edit_block=edit_block;
        gap_align->edit_block->length1 = gap_align->query_length;
        gap_align->edit_block->length2 = gap_align->subject_length;
	return TRUE;
}


sap_typ	gpsi_type::FakeGBlastEngineCoreCMSA(Int4 Threshold,BooLean srch_subset, cma_typ cma)
// use cma_typ to initialize the sap_typ
{
  Int4		r,sq2,n_ug;
  time_t	time1;
  Int4		nhits,total_hits,word_extend_dropoff;
  e_type	sE;
  sap_typ	sap,sap2;
  hhp_typ	hhp;		// Temporary heap for HSPs.
  hsp_typ	hsp;
  double	key;
  // double	gap_trigger_bits = 22.0;  // ungapped cutoff for triggering gapxdrop
  // AFN: need to modify gap_trigger_bits when evalue is high & database is small!!!

  assert(NSeqsSeqSet(seqset) >= 1);
  time1=time(NULL); 
  // PutSeq(stdout,queryE,AB);
  // change to: GABNew(queryE,gap_open,gap_extend,SMatrixSBP(search->sbp)); ???
  gab_typ gap_align= GABNew(gap_open, gap_extend);
  gap_align->include_query = 0;	// query length from q_start that MUST be aligned.
  gap_align->query = SeqPtr(queryE) + 1;
  gap_align->query_length = LenSeq(queryE);
  gap_align->decline_align = (-(GBLAST_SCORE_MIN));
  gap_align->matrix=SMatrixSBP(search->sbp);  // NOTE: may need to change sbp->matrix!!!

  SetStdStatsSBP(search->sbp);	// use standard statistics
  word_extend_dropoff=WordXDropoffSBP(word_extend_dropoff_in_bits,search->sbp);
  gap_trigger = GapTriggerSBP(gap_trigger_bits,search->sbp);
  gap_align->posMatrix=NULL; gap_align->positionBased = FALSE; 
  gap_x_dropoff_final = GapXDropoffSBP(x_parameter_final,search->sbp);
  gap_x_dropoff = GapXDropoffSBP(x_parameter,search->sbp);
  gap_align->x_parameter=gap_x_dropoff_final;
  if (thisPassNum > 1) { ResetGBlastSearchBlk(search); FreeHead( ); }
  rhp_typ rhp = MakeBRHHeap(hitlist_size);
  brh_typ results; total_hits=0;
  for(n_ug=0,sq2=1; sq2<=NSeqsSeqSet(seqset); sq2++){
      if(sq2 > NumSeqsCMSA(cma)) break;	// search only the sequences in cma.
      sE = SeqSetE(sq2,seqset);
      // fprintf(stderr,"*********** Fake search of sequence %d ************\n",sq2);
      if(sq2==1){
	nhits=1; // HITS always == 1 for cma with one block (for now...)
      	total_hits+=nhits;
      	assert(total_hits <= hitlist_size); // AFN: my addition.
      	gap_align->subject = SeqPtr(sE) + 1;
        gap_align->subject_length = LenSeq(sE);
        hhp = MkHSPHeap(nhits+2); 
        BooLean success=TRUE;
	assert(FakeIdenticalTracebackCMSA(gap_align,(2*NumSeqsCMSA(cma)),sE));
	// hsp=MakeHSPTyp(gap_align,search->sbp);
	hsp=MakeHSPTyp(gap_align,1e-20,(2*NumSeqsCMSA(cma)));
	// PutHSPTyp(stderr,hsp);
	assert(InsertHSPHeap(hsp, hhp) != NULL);
        AddSet(sq2,SubSet);	// Keep track of potential hits.
        Int4 n=PurgeHSPHeap(hhp); 
        if(n > 0){
	  results = MakeGBLASTResultHitlist(n+1,sE);
	  while((hsp=DelMinHSPHeap(&key, hhp)) != 0){
		assert(results->subject_id == sq2);
 		AddHspBRH(hsp,results);
	  }
          if(InsertBRHHeap(results,rhp)==NULL) results=NilBRH(results);
        }
      } else {
	nhits=nBlksCMSA(cma); // HITS == number of motifs in cma...
        total_hits+=nhits;
      	assert(total_hits <= hitlist_size); // AFN: my addition.
      	gap_align->subject = SeqPtr(sE) + 1;
      	gap_align->subject_length = LenSeq(sE);
      	hhp = MkHSPHeap(nhits+2); 
      	BooLean success=TRUE;
      	for(r=1; r<=nhits; r++){
#if 0	// Doesn't look like gap_align->q_start & gap_align->s_start are needed...
	  gap_align->q_start=s+d+opt_offset;
	  gap_align->s_start=s+opt_offset;
#endif 
	// KEY FUNCTION IS NEXT; IF I CAN GET THIS RIGHT THEN ALL ELSE SHOULD WORK.
#if 1
	 if(FakeGappedAlignmentWithTracebackCMSA(gap_align,sq2,queryE,sE,r,
				 (2*NumSeqsCMSA(cma))-sq2,cma)){
	  // Added this to accommodate Viterbi insert followed by delete and vvs.
	  //hsp=MakeHSPTyp(gap_align,search->sbp);
	  hsp=MakeHSPTyp(gap_align,1e-20,(2*NumSeqsCMSA(cma)-sq2));
	  // PutHSPTyp(stderr,hsp);
	  assert(InsertHSPHeap(hsp, hhp) != NULL);
	 } else {
		fprintf(stderr,"sequence %d failed to gap align\n",sq2);
		success=FALSE;
	 }
#else
	 assert(FakeGappedAlignmentWithTracebackCMSA(gap_align,sq2,queryE,sE,r,
				 (2*NumSeqsCMSA(cma))-sq2,cma));
	 // Added this to accommodate Viterbi insert followed by delete and vvs.
	 hsp=MakeHSPTyp(gap_align,search->sbp);
	 // PutHSPTyp(stderr,hsp);
	 assert(InsertHSPHeap(hsp, hhp) != NULL);
#endif
        }
        if(success){
          AddSet(sq2,SubSet);	// Keep track of potential hits.
	  // fprintf(stderr,"sequence %d within PurgeHSPHeap\n",sq2);
          Int4 n=PurgeHSPHeap(hhp); 
          if(n > 0){
	    // fprintf(stderr,"sequence %d put into Results\n",sq2);
	    results = MakeGBLASTResultHitlist(n+1,sE);
	    while((hsp=DelMinHSPHeap(&key, hhp)) != 0){
		// assert(results->subject_id == sq2); // don't worry about this...
 		AddHspBRH(hsp,results);
	    }
            if(InsertBRHHeap(results,rhp)==NULL) results=NilBRH(results);
          } else fprintf(stderr,"sequence %d NOT put into Results\n",sq2);
        }
       } NilHSPHeap(hhp);  // delete hhp regardless (sq2==1 or sq2 == subject...)
     } head=0; 
  while((results=DelMinBRHHeap(&key,rhp)) != 0){
        sap = ExtractAlnBRH(results,queryE);
	// PutGSeqAlign(stderr, sap, 60, AB);
        if(head == 0) { head = sap; }
        else{ for(sap2=head;sap2->next;) sap2=sap2->next; sap2->next=sap; }
        Int4 seqIndex=search->hitlist_count;
        assert(search->hitlist_max > search->hitlist_count);
        search->hitlist_count++;
        search->results[seqIndex]=results;
  } NilBRHHeap(rhp); GABDelete(gap_align);
  if(verbose) PutHitList(stderr,60);
  // PutGSeqAlignList(stderr, head, 60, AB);	// DEBUG...
  // PutMultiGSeqAlign(stderr, head, 60, AB);	// DEBUG...
  if(verbose){
    if(!srch_subset) fprintf(stderr,"%d sequences extended; %d detected\n",
		CardSet(SubSet),search->hitlist_count);
     else fprintf(stderr,"second pass: %d sequences detected\n",
		search->hitlist_count);
     fprintf(stderr,"time = %0.1f seconds\n",difftime(time(NULL),time1));
  }
  return head;
}

//********************************* Information routines ********************
double          **gpsi_type::CopyPosPseudoFreq( )
{
       if(posSearch && posSearch->posFreqs){
           double **pm,**PM;
           NEWP(PM,LenSeq(queryE)+3,double);
           pm=PM+1;                // start at 1
           for(Int4 j,i=0; i< LenSeq(queryE); i++){
               NEW(PM[i+1],nAlpha(AB)+2,double);
               double total=0.0;
               for(j=0; j <= nAlpha(AB); j++){
                       total+=pm[i][j] = posSearch->posFreqs[i][j];
               }
               if(total==0.0) { free(pm[i]); pm[i]=0; }
               else { // renormalize...
                // WARNING: PSI-BLAST is corrently not normalizing this properly!!
                    for(j=0; j <= nAlpha(AB); j++){
                       pm[i][j] = pm[i][j]/total;
                    }
                    //assert(total > 0.999 && total < 1.001);
               }
           } return PM;
       } else return 0;
}

double   **gpsi_type::CopyObsWtCounts( )
{
         if(posSearch && posSearch->posMatchWeights){
             double **pm,**PM;
             NEWP(PM,LenSeq(queryE)+3,double);
             pm=PM+1;                // start at 1
             for(Int4 i=0; i< LenSeq(queryE); i++){
                 NEW(PM[i+1],nAlpha(AB)+2,double);
                 for(Int4 j=0; j <= nAlpha(AB); j++){
                     pm[i][j] = posSearch->posWtCount[i] *
                         posSearch->posMatchWeights[i][j];
                 }
#if 0   // DEBUG
        if(i==96){
                for(Int4 x=0; x <= nAlpha(AB); x++){
                        fprintf(stderr,"96 %c: WtCounts = %.3f ; posMatchWts =%.3f \n",
                                AlphaChar(x,AB),
                                posSearch->posWtCount[i],
                                posSearch->posMatchWeights[i][x]);
                }                  
                
        }
#endif
             } return PM;
        } else return 0;
}

double  **gpsi_type::CopyPosMatchWeights( )
{
       if(posSearch && posSearch->posMatchWeights){
          double **pm,**PM;
          NEWP(PM,LenSeq(queryE)+3,double);
          pm=PM+1;                // start at 1
          for(Int4 i=0; i< LenSeq(queryE); i++){
               NEW(PM[i+1],nAlpha(AB)+2,double);
               for(Int4 j=0; j <= nAlpha(AB); j++){
                   pm[i][j] = posSearch->posMatchWeights[i][j];
               }
          } return PM;
       } else return 0;
}

BooLean  gpsi_type::FractionSeqsAligned(double *fract_aligned)
{
       if(posMatrix == 0) return FALSE;
       else {
            // compute number of sequences actually used...
            Int4 N=1; // includes query...
            for(Int4 i=1; i<=posSearch->posNumSequences; i++){
                if(posSearch->posUseSequences[i]) N++;
            }
            for(Int4 s=1; s <= LenSeq(queryE); s++){
                double tmp=((double) posSearch->posCount[s-1]/(double) N);
                // fprintf(stderr,"%d: %.3f\n",s,tmp);
                fract_aligned[s]=tmp;
            }
       } return TRUE;
}

double  *gpsi_type::AveSeqWeights( )
{
	double	*avesqwt=0;
	if(posMatrix){
            // compute number of sequences actually used...
	    NEW(avesqwt,posSearch->posNumSequences+3,double);
            Int4 N=1; // includes query...
            for(Int4 i=1; i<=posSearch->posNumSequences; i++){
                if(posSearch->posUseSequences[i]) N++;
            }
            for(Int4 s=1; s <= LenSeq(queryE); s++){
                double tmp=((double) posSearch->posCount[s-1]/(double) N);
                // fprintf(stderr,"%d: %.3f\n",s,tmp);
                avesqwt[s]=tmp;
            }
	} return avesqwt;
}

