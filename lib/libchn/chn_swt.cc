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

#include "chn_typ.h"

void	chn_typ::GetMergedIntegerSWT(cma_typ qcma, cma_typ fg_cma, cma_typ bg_cma)
{
	// 0. Make sure that this is being called from chn_see with a FG & BG aln...
	assert(fg_cma && bg_cma);
	// 1. Merge FG and BG alignments into one Main alignment.
	cma_typ	tmpCMAs[5];
#if 0	// if switch off, then need 'N=NumSeqsCMSA(fg_cma);' below.
	tmpCMAs[1]=qcma; tmpCMAs[2]=fg_cma; tmpCMAs[3]=bg_cma; 
	FILE *fp = tmpfile(); PutMergedCMSA(fp,3,tmpCMAs); rewind(fp);
#elif 1
	tmpCMAs[1]=fg_cma; tmpCMAs[2]=bg_cma; 
	FILE *fp = tmpfile(); PutMergedCMSA(fp,2,tmpCMAs); rewind(fp);
#else
	tmpCMAs[1]=fg_cma; tmpCMAs[2]=bg_cma; 
	FILE *fp = tmpfile(); PutMergedCMSA(fp,2,tmpCMAs); rewind(fp);
#endif
	MainCMA=ReadCMSA(fp,AB); fclose(fp);

	// 2. create new swt for the full alignment
	swt[0] = new swt_typ(qcma,MainCMA,use_ocma_pseudo,GlobalSqWts);
	// swt[0]->Put(stderr);
	
	// 3. compute sequence weights as for pps_typ for Main align...
	SqIWtMain=swt[0]->GetIntegerWts(&MainSqWt);

	// 4. Split Array into two for Observed counts for FG & BG using these weights.
	Int4	N,M;
	Int4	s,m,f,b,r;
	// N=NumSeqsCMSA(qcma) + NumSeqsCMSA(fg_cma); 
	N=NumSeqsCMSA(fg_cma); 
	M=NumSeqsCMSA(bg_cma);
	if(N+M != NumSeqsCMSA(MainCMA)){
		fprintf(stderr,"N=%d; M=%d; N+M=%d\n",N,M,NumSeqsCMSA(MainCMA));
		print_error("this should not happen");
	}
	TotalFG_N=TotalBG_M=0;
	double wt_factor=1.0/(double)swt[0]->WtFactor();
	// Initialize Obs[1][s][r] and Obs[2][s][r]
	UInt4	**ResWtFG,**ResWtBG;
	// NEWP(ResWtFG,LengthCMSA(1,MainCMA)+3,UInt4); NEWP(ResWtBG,LengthCMSA(1,MainCMA)+3,UInt4);
	for(s=1; s <= LengthCMSA(1,MainCMA); s++){
		for(r=0; r <= nAlpha(AB); r++) Obs[1][s][r] = Obs[2][s][r]=0.0;
		// NEW(ResWtFG[s],nAlpha(AB)+5,UInt4); NEW(ResWtBG[s],nAlpha(AB)+5,UInt4);
	}
	for(m=f=1; m <= N; f++,m++){
	   TotalFG_N += SqIWtMain[m]; 
	   // Recalculate FG Observed Weighted counts
	   unsigned char *seq=GetAlnResInSiteCMSA(1,m,MainCMA);
	   for(s=1; s <= LengthCMSA(1,MainCMA); s++){
		Obs[1][s][seq[s]] += wt_factor*MainSqWt[s][m]; // in che_typ.cc : wt = SqWt[s][sq];
	   }
	}
	for(b=1; b <= M; b++,m++){
	   TotalBG_M+=SqIWtMain[m]; 
	   // Recalculate BG Observed Weighted counts
	   unsigned char *seq=GetAlnResInSiteCMSA(1,m,MainCMA);
	   for(s=1; s <= LengthCMSA(1,MainCMA); s++){
		Obs[2][s][seq[s]] += wt_factor*MainSqWt[s][m]; // in che_typ.cc : wt = SqWt[s][sq];
	   }
	}
	// Compute new WtFrq[s][r];
	double verytiny=0.0000000000001;
	double tiny = verytiny * (double) nAlpha(AB);
	for(s=1; s <= LengthCMSA(1,MainCMA); s++){
	   for(r=StartAlpha; r<=nAlpha(AB); r++){
#if 0
               WtFrq[1][s][r] = (Obs[1][s][r]+1.0)/(wt_factor*TotalFG_N);
               WtFrq[2][s][r] = (Obs[2][s][r]+1.0)/(wt_factor*TotalBG_M);
#else		// afn 6/24/09: add a tiny number of pseudocounts to avoid problems with sum freq >> 1.0
               WtFrq[1][s][r] = (Obs[1][s][r] + verytiny)/((wt_factor*TotalFG_N)+tiny);
               WtFrq[2][s][r] = (Obs[2][s][r] + verytiny)/((wt_factor*TotalBG_M)+tiny);
#endif
	   }
        }
}

#if 0	// moved to cmsa.h
static Int4 pseudo_aln_score_chn_see(register Int4 Length, register unsigned char *seq1,
        register unsigned char *seq2, register char **R)
{
        register Int4 score;
        for(score=0; Length > 0; Length--,seq1++,seq2++){
                // Skip over gap residues...
                if(*seq1 && *seq2) score += R[*seq1][*seq2];
        }
        return score;
}

Int4    PseudoAlnScoreCMSA_chn_see(e_type E, Int4 sq2, cma_typ cma)
// use cma alignment to obtain pseudo pairwise scores for two aligned sequences.
{
        assert(nBlksCMSA(cma) ==1);
        assert(LenSeq(E) == LengthCMSA(1,cma));
        Int4 N=NumSeqsCMSA(cma);
        assert(sq2 > 0 && sq2 <= N);
        a_type  A=AlphabetCMSA(cma);
        unsigned char *seq1 = SeqPtr(E);
        unsigned char *seq2 = SeqSeqSet(sq2,DataCMSA(cma));
        Int4 s2 = SitePos(1,sq2,1,SitesCMSA(cma));
        return pseudo_aln_score_chn_see(LengthCMSA(1,cma),seq1+1,seq2+s2,AlphaR(A));
}
#endif

// void	chn_typ::GetSWT(char *NameFReadSWT)
void	chn_typ::GetSWT(hsw_typ *HSW)
// afn: 11_9_07.
// get information and compute sequence weights...
{

  double	**Observed=0,**nullfreq=0,**WtFreq=0,**NewFreq=0,**OFreq=0;
  Int4		i,j,s,m,n,r;
  FILE		*fp;
 	cma_typ cma=IN_CMA[1];

//******************* Initialize primary alignments *************************
      fp=tmpfile();
      Int4 length = LengthCMSA(1,CMA[1]);
      PutAlnCMSA(fp,CMA[1]); rewind(fp); MA=ReadSMA(fp); fclose(fp);
      if(MA == NULL){
	print_error("SMA read error: check input cma file for syntax errors (e.g., defline)");
      }
      if(Number > 2){
            for(j=1,s=3; j <= NumSeqAln; j++,s++){
                if(extra_files && verbose) fprintf(stderr,"len(%d) = %d (%d)\n",
                        s,LengthCMSA(1,IN_CMA[s]),LenSeq(TrueSeqCMSA(j,cma)));
            }
      }
      NEW(Kingdom,NumSeqsCMSA(CMA[1])+3,char);
      NEWP(Phylum,NumSeqsCMSA(CMA[1])+3,char);
      assert(NumSeqsCMSA(CMA[1]) == nseqSMA(MA));
      for(j=1; j<=NumSeqsCMSA(CMA[1]); j++){
                Kingdom[j]=KingdomSeq(TrueSeqCMSA(j,CMA[1]));
                Phylum[j]=PhylumSeq(TrueSeqCMSA(j,CMA[1]));
      }
      Int4 gapo=11,gapx=1,minscore=INT4_MAX;
      if(extra_files && verbose) fprintf(stderr,"seq  score\n");
      for(s=1; s <= NumSeqsCMSA(cma); s++){
                e_type tmpE=TrueSeqCMSA(s,cma);
                Int4 score = FastAlnSeqSW(gapo,gapx,qE,tmpE,AB);
                if(score < minscore) minscore=score;
                if(extra_files && verbose) fprintf(stderr,"%3d: %d\n",s,score);
      }

//***************** Main loop for gathering info for each analysis *****************
	// fprintf(stderr,"Computing relevant information.\n");
//***************** Main loop for gathering info for each analysis *****************
//***************** Main loop for gathering info for each analysis *****************

  NEWP(swt,NumberMCMA+3,swt_typ);
  for(Int4 b=1,anal=0; anal <= NumberMCMA; b++,anal++){
     if(anal > 0){	//****** Deal with hierarchical alignments first.
	Int4 num_rm=0;
	BooLean	*skipseq=0;
        //************* remove sequences close to consensus (query) sequence ************
        //************* remove sequences close to consensus (query) sequence ************
        //************* remove sequences close to consensus (query) sequence ************
	if(percent_cut > 0.0 && percent_cut <= 100.0){
	  NEW(skipseq,NumSeqsCMSA(MCMA[anal])+2,BooLean);
	  for(Int4 sq=1; sq <=NumSeqsCMSA(MCMA[anal]); sq++){
	   Int4 ident=0;
	   double per_id;
#if 1	// Speed up ResidueCMSA( )...
	   Int4 s0=SitePos(1,sq,1,SitesCMSA(MCMA[anal])); // start of fake sequence.
           unsigned char *fsq = SeqSeqSet(sq,DataCMSA(MCMA[anal]));
	   unsigned char *qsq = SeqPtr(qE);
	   for(s=1; s <= LenSeq(qE); s++,s0++){
		if(qsq[s] != 0 && qsq[s] == fsq[s0]) ident++;
	   } 
#else 
	   for(s=1; s <= LenSeq(qE); s++){
		Int4 r1 = ResSeq(s,qE);
		Int4 r2 = ResidueCMSA(1,sq,s,MCMA[anal]);
		if(r1 != 0 && r1 == r2) ident++;
	   } 
#endif
	   per_id=100.0*(double)ident/(double)LenSeq(qE);
	   if(per_id >= percent_cut){ num_rm++; skipseq[sq]=TRUE; }
	  }
	  if(extra_files) fprintf(stderr,"%.1f%% identical sequences removed or %d/%d\n",
			percent_cut,num_rm,NumSeqsCMSA(MCMA[anal]));

	} else if(percent_cut == 0.0){
	  h_type HG=Histogram("scores of aligned sequences",-100,1000,5.0);
          NEW(skipseq,NumSeqsCMSA(MCMA[anal])+2,BooLean);
	  for(Int4 sq=1; sq <=NumSeqsCMSA(MCMA[anal]); sq++){
		Int4 score = ConsensusScoreCMSA(MCMA[anal],qE,gapo,gapx,sq);
		if(score >= minscore){ num_rm++; skipseq[sq]=TRUE; }
		IncdHist(score,HG);
	  }
	  if(extra_files) fprintf(stderr,"%d/%d sequences removed at or above a score of %d\n",
                        num_rm,NumSeqsCMSA(MCMA[anal]),minscore);
	  // PutHist(stdout,60,HG); 
	  NilHist(HG);

	}
        if((NumSeqsCMSA(MCMA[anal])-num_rm) < 1) print_error("too many sequences removed!\n");
	if(skipseq){ // new file read in.
          fp=tmpfile(); PutSelectCMSA(fp,skipseq,MCMA[anal]); free(skipseq);
          TotalNilCMSA(MCMA[anal]);
          rewind(fp); MCMA[anal]=ReadCMSA(fp,AB); fclose(fp);
	}
        //********** END: remove sequences close to consensus (query) sequence *********
        //********** END: remove sequences close to consensus (query) sequence *********

        //************* BEGIN: Get weighted sequence information here. ************
        //************* BEGIN: Get weighted sequence information here. ************
        //************* BEGIN: Get weighted sequence information here. ************
	// fprintf(stderr,"Computing sequence weights for alignment #%d.\n",anal);
	// NOTE: swt_typ won't work for background block-based cma files!
	if(anal==1){
	   CMA[b]=MCMA[anal];
	   if(IN_CMA[2]!=MCMA[anal]) IN_CMA[2]=MCMA[anal];
#if 0
	   swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,GlobalSqWts);
#else	   //
	   // swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,GlobalSqWts,NameFReadSWT);
	   // swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,GlobalSqWts,HSW);
	   if(HSW && HSW[anal]){
		swt[anal] = new swt_typ(HSW[anal]); 
		assert(swt[anal]->RtnHSW_CMA() == MCMA[anal]);	// make sure this is the same Main aln.
	   } else { swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,GlobalSqWts); }
#endif
	   if(!extra_files) swt[anal]->Silent( );
	} else {  // i.e., anal >= 2.
#if 0	   // Test out using pseudo alignment scores...
	   e_type	tmpE=MkConsensusCMSA(MCMA[anal]); // 
	   Int4		max_i=0,max_score=0;
	   for(i=1; i <= NumSeqsCMSA(MCMA[anal-1]); i++){
		Int4 scr = PseudoAlnScoreCMSA_chn_see(tmpE,i,MCMA[anal-1]);
		if(scr > max_score){ max_score=scr; max_i=i; }
	   }
	   swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,
				FakeSeqCMSA(max_i,MCMA[anal-1]),GlobalSqWts);
#else
	   if(HSW && HSW[anal]){
		swt[anal] = new swt_typ(HSW[anal]); 
		assert(swt[anal]->RtnHSW_CMA() == MCMA[anal]);	// make sure this is the same Main aln.
	   } else {
	        e_type tmpE=MkConsensusCMSA(MCMA[anal-1]); // 
	   	swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,tmpE,GlobalSqWts);
	   	NilSeq(tmpE);
	   }
#endif
	   if(!extra_files) swt[anal]->Silent( );
	   CMA[b]=IN_CMA[NumSeqAln+1+anal]=MCMA[anal];
	}
        //************* END: Get weighted sequence information here. ************
        //************* END: Get weighted sequence information here. ************
	// if(extra_files) { swt[anal]->Put(stdout); swt[anal]->Put(stderr,180); }
	// swt[anal]->Test(stdout);
	//************** Compute Observed counds and background frequencies *******************
	//************** Compute Observed counds and background frequencies *******************
	// fprintf(stderr,"Computing residue frequencies.\n");
	if(NoWeights){
	  BooLean *skip=0;
	  WtFreq=ColResFreqsCMSA(1,skip,&Observed,MCMA[anal]);
	} else {
	    if(compute_mpwf){	// Get MargProb WtFreq 
		WtFreq=swt[anal]->MargProbWtFreq( );
	        if(WtFreq==0){
			fprintf(stderr,"WARNING: marginal probability numerical error\n");
			fprintf(stderr,"  ...using column specific weights instead.\n");
			WtFreq=swt[anal]->WeightedFreq( );
		}
	    } else {
// ******************** KEY SECTION RIGHT HERE... ************************
		WtFreq=swt[anal]->WeightedFreq( ); 
	    }
	    Observed=swt[anal]->ObsWtCnts( ); 
// ******************** END KEY SECTION. ********************
	    // fprintf(stderr,"getting ObsWtCnts( ) for analysis %d...\n",anal);
	}
#if 0
	MargProb[b]=swt[anal]->MargProbAln( ); // b == anal+1;
#else 
	// NOTE that swt[anal+1] returns MargProb for consensus of anal versus anal+1.
        // e_type tmpE=MkConsensusCMSA(MCMA[anal-1]); //
        // swt[anal] = new swt_typ(cma,MCMA[anal],use_ocma_pseudo,tmpE);
	// AFN: 3/9/06: this MargProbAln() routine is broke (core dumps)
	// problems in both chn_cnt and chn_see.
	// Fix this part of code when get a chance.
	if(compute_mpwf){	// Get MargProb WtFreq 
	  print_error("Fatal: MargProbAln() broken");
	  Int4 anal_m1=anal-1;
	  switch(AlignCode[anal_m1])
	  {	// Note: MargProb[anal] actually corresponds to anal-1;
	   case 'Q': MargProb[anal]=0; break;
	   case 'F': MargProb[anal]=swt[anal_m1]->MargProbAln( );
	     // this uses consensus of anal-1 (i.e., MargProb[anal]) versus alignment of anal?
	     break;
	   case 'B': MargProb[anal]=0; break;
	   case 'C': MargProb[anal]=swt[anal]->MargProbAln( );
	     // this uses consensus of anal-1 (i.e., MargProb[anal]) versus alignment of anal.
	     break;
	   case 'S': MargProb[anal]=0; break; // swt[anal+1] doesn't exist for 'S'.
	   default: MargProb[anal]=0; break;
	  }
	}
#endif
	//***************** Base Freqs on Marginal Probabilities **********************
	//***************** Base Freqs on Marginal Probabilities **********************
	//***************** Base Freqs on Marginal Probabilities **********************
	// WARNING: The following 'if ' needs lots of work!
	if(!NoWeights && mp_cutoff > 0.0){
	   print_error("This MargProb section still under development!\n");
	   double *FractAln=swt[anal]->FractSeqAln();
	   NEWP(NewFreq,LenSeq(qE)+3,double);
	   NEWP(OFreq,LenSeq(qE)+3,double);
	   double **tmpFreq=0;
	   tmpFreq=WtFreq;
	   double **obs;
	   OFreq=ColResFreqsCMSA(1,0,&obs,cma);
	   for(s=1; s<=LenSeq(qE); s++){
		if(MargProb[b] && MargProb[b][s] >= mp_cutoff 
		   && FractAln[s] >= fsa_cutoff) NewFreq[s] = tmpFreq[s];
		else { NewFreq[s] = OFreq[s]; }
	   } WtFreq=NewFreq;
	   free(obs);
	}
	if(verbose) CheckWeightedMatrix(qE,WtFreq, AB);
	//*************** Set Main Set Foreground & background frequencies ******************
	//*************** Set Main Set Foreground & background frequencies ******************
// ******************* KEY SECTION RIGHT HERE... *********************************
    	NEWP(nullfreq,length+3,double);
    	for(i=1; i <= length; i++) { 
      	  nullfreq[i]=tFreqSeqSet(TrueDataCMSA(MCMA[anal]));
	} // nullfreq[0] == NULL indicates that arrays should not be freed.
	Obs[anal]=Observed; NullFreq[anal]=nullfreq;
	WtFrq[anal]=WtFreq; // WtFrq is passed to rtf_typ for background frequencies...
// ******************** END KEY SECTION. *******************************
        // fprintf(stderr,"**** output_name(5)='%s' ****\n",output_name);
        //*********** Modify background frequencies if block-based file provided ************
        //*********** Modify background frequencies if block-based file provided ************
        if(BG_CMA && anal == NumberMCMA) GetSwtBG(anal);
   } else { //****  anal == 0: Deal with regular family alignment here.
	// Appears to be for first (query family) alignment.
       double	**observed; 
      if(NumAnalysis > 1){
        NEWP(observed,length+2,double);
        NEWP(nullfreq,length+3,double);
    	for(i=1; i <= length; i++) { 
      	  nullfreq[i]=tFreqSeqSet(TrueDataCMSA(CMA[1]));
      	  NEW(observed[i],nAlpha(AB)+2,double);
      	  for(n=1; n<=NumSeqsCMSA(CMA[1]); n++){ 
		r=ResidueCMSA(1,n,i,CMA[1]); observed[i][r]+=1.0;
      	  }
	} // nullfreq[0]=nullfreq[1]; // indicates that arrays should be freed.
	NullFreq[0]=nullfreq;
	Obs[0]=observed;
	WtFrq[0]=Freq[0]=0;
	BackGrnd[0]=0; 
	// BackGrnd[0]=1; 
	FractSeqAln[0]=0; 
	NumSqMain[0]=0; 
	Hist[0] = FALSE; 
	SuperAln[0]=FALSE;  // Freq[] == 0 if no SuperAln, so don't need this.
      } else { // NumAnalysis == 1
	assert(!"Not implemented for NumAnalysis == 1!");
      }
    }  // end of dealing with regular family alignment.
  } // end of anal <= NumberMCMA loop
//************************** moved here from chn_rtf.cc ****************************
     if((ModeLPR <= '9' && ModeLPR >= '0') && NumberMCMA == 2){	// --> query + FG + BG...
	GetMergedIntegerSWT(cma, MCMA[1],MCMA[2]);	// swt[0] = original Main set
     }
     if(NumAnalysis <= 1){ // then show consensus alignment only.
		assert(!"This option not yet implemented");
     } else {
// ******************* KEY SECTION RIGHT HERE... *********************************
	Int4 a,anal;
	for(a=2,anal=1; anal <= NumberMCMA; anal++,a++){
	   if(Obs[anal]){
		WtNumSeq[a]=GetWtNumSq(&WtNumSq[a],Obs[anal],length);
		// WtNumSeq = used for residue highlighting in rtf! (see rtf_typ.cc)
		// WtNumSq = printed in rtf file only (not used for calculations).
	   } else {
		WtNumSq[a]=0; WtNumSeq[a]=0;
	   }
	   // fprintf(stderr,"WtNumSq[%d] = %g\n",a,WtNumSq[a]);
	}
	if(BG_CMA) WtNumSeq[a]=GetWtNumSq(&WtNumSq[a],ObsBG,length);
        // NOTE: above not weighting sequences; for this need to call Henikoff funct.
// ******************* END KEY SECTION... *********************************
     }
//************************** moved here from chn_rtf.cc ****************************
}

void	chn_typ::CheckWeightedMatrix(e_type qE, double  **WtFreq, a_type A)
{
	Int4	i,r;
	fprintf(stderr,"\nWeighted counts for PSI-BLAST\n");
	for(i=1; i <= LenSeq(qE); i++){
	    if(i % 50 == 1 && i < LenSeq(qE)){
		fprintf(stderr,"\npos: ");
		for(r=StartAlpha; r <= nAlpha(A); r++){
		   fprintf(stderr,"  %c ",AlphaChar(r,A));
		} fprintf(stderr,"\n");
	    }
	    if(WtFreq[i]){
	      fprintf(stderr,"%3d:",i);
	      for(r=StartAlpha; r <= nAlpha(A); r++){
		if(WtFreq[i][r] == 0.0) fprintf(stderr,"   ."); 
		else fprintf(stderr," %3d",(Int4) floor((100.0*WtFreq[i][r])+0.5));
	      } fprintf(stderr,"\n");
	    } else fprintf(stderr,"%3d: (nil)\n",i);
	} fprintf(stderr,"\n");
}

double	*chn_typ::GetWtNumSq(double *RtnWtNSq,double **Obsvd, Int4 length)
{
        h_type HG=Histogram("WtNumSq",-1,LenSeq(keyE)+OffSetSeq(keyE)+1,1.0);
        double WtNSq = 0.0;
	double *WtNSeq;

	NEW(WtNSeq,length+5,double);
	for(Int4 i=1; i <= length; i++){
          double wt_nsq=0.0,d=0.0;
          for(Int4 r=StartAlpha; r <= nAlpha(AB); r++){
		wt_nsq += Obsvd[i][r];	// This is most likely from WtCnts in swt_typ.
		if(r == ResSeq(i,keyE)) d=Obsvd[i][r];
	  } 
	  WtNSq += wt_nsq;
	  WtNSeq[i]=wt_nsq;
	  // d = 100.0*d/wt_nsq;	// d == query residues...
	  d = ceil(wt_nsq);
	  if(d > 0.0) IncdMHist(i+OffSetSeq(keyE),(Int4)d,HG);
	} WtNSq = WtNSq/(double) length;
	// PutHist(stderr,60,HG); 
	NilHist(HG);
	*RtnWtNSq=WtNSq;
	return WtNSeq;
}

void	chn_typ::GetSwtBG(Int4 anal)
//************************ Motif based input alignment ***********************
{
	assert(BG_CMA && anal == NumberMCMA);
	Int4	length = LengthCMSA(1,CMA[1]);
	Int4	i,m,r;
	// If user has input a cma background frequency file...
	// then go back and replace the nullfreq's at aligned positions.
	Int4	col,blk,j;
	char	*operation=0;

	// MCMA[NumberMCMA+1]=BG_CMA; // ??? attempt to fix alignment name problem...
	if(extra_files) fprintf(stderr,"Using input background frequency cma file\n");
	// 1. First check to make sure that the query sequence contains the motifs.
	Int4	BlkScore;
	Int4    *mtfpos=AlignSeqCMSA(stderr,&BlkScore,keyE,BG_CMA);

#if 1	//************************** smatrix... **********************************
	Int4 **gapscore=NULL;	// Don't use this right now...
	e_type	alnE;
#if 0
	alnE=FakeSeqCMSA(18,cma);
	AlnSeqSW(11,1,alnE,keyE,AB);
	PutSeq(stderr,alnE,AB);
	PutSeq(stderr,keyE,AB);
#else
	// alnE=keyE;
	alnE=qE;
#endif
	assert(LenSeq(alnE) == LenSeq(keyE));

	fprintf(stderr,"/*****************************************************/\n");
	operation=GapAlignSeqCMSA(stderr,insert,extend,&BlkScore,alnE,gapscore,BG_CMA);

std::cerr << "tmp traceback"; std::cerr << std::endl; 
std::cerr << operation; std::cerr << std::endl;

  // Traceback returned is reversed and ends at start; needs to be reversed

#endif	//************************** smatrix... **********************************

	double	pernats = SitesGSS(SitesCMSA(BG_CMA))->PerNats();
	double	AdjScore = 1.0/exp((double) BlkScore/pernats);
	Int4	netleng = length - TotalLenCMSA(BG_CMA);
	if(netleng > 0) AdjScore = AdjScore*bico(netleng+nBlksCMSA(BG_CMA),nBlksCMSA(BG_CMA));
	fprintf(stderr,"Adjusted probability = %.3g\n",AdjScore);
	if(FALSE && AdjScore > 0.05){	// about 0.02 significance without adjusting for length...
		print_error("Input cma file fails to find a significant match to query");
	}

    	NEWP(ObsBG,length+3,double);
    	NEWP(WtFrqBG,length+3,double);
    	NEWP(NullFreqBG,length+3,double);
	FractHighlightBG = (double) (0.20*TotalLenCMSA(BG_CMA))/(double)length;
	fprintf(stderr,"FractHighlightBG = %f\n",FractHighlightBG);
	// highlight half of the columns in the motif regions...

	NEW(FractSeqAlnBG,length+3,double);
	NEW(ObsBG[0],nAlpha(AB)+3,double);
      	NullFreqBG[0]=tFreqSeqSet(TrueDataCMSA(MCMA[anal]));
	for(r=0; r <=nAlpha(AB); r++) ObsBG[0][r] = NullFreqBG[0][r]*NumSeqsCMSA(BG_CMA);

	//******** First get default values for non-motif positions.. **********
	//******** First get default values for non-motif positions.. **********
	//******** First get default values for non-motif positions.. **********
    	for(i=1; i <= length; i++) { 
		// assert(ResSeq(i,keyE) == ResSeq(i,fakeSq1));
      	  	NullFreqBG[i]=NullFreqBG[0];
      	  	WtFrqBG[i]=NullFreqBG[0];
		ObsBG[i]=ObsBG[0];
		FractSeqAlnBG[i]=0.0;
	} // NullFreqBG[0] == NULL indicates that arrays should not be freed.

	//****************** Get values for motif positions.. ********************
	//****************** Get values for motif positions.. ********************
	//****************** Get values for motif positions.. ********************
	// PutSeq(stderr,keyE,AB); PutSeq(stderr,fakeSq1,AB);
	double **TmpObs,**TmpWtFreq,**col_freq;
	st_type S=SitesCMSA(BG_CMA);
	Int4 o,lenM,site,end;
	char	state=0;
	if(operation){
	  for(m=0,j=0,o=1; operation[o]; o++){
	    state=operation[o];
	    switch(state){
		case 'M': 	// match at start of motif block...
		  j++;
		case 'D':	// deletion at beginning of motif..
	  	  if(m > 0){ free(col_freq); free(TmpObs); free(TmpWtFreq); }
		  m++;
	  	  lenM = SiteLen(m,S); col=1;
	  	  col_freq=ColResFreqsCMSA(m, BG_CMA);
	  	  TmpWtFreq=ColResFreqsCMSA(m, &TmpObs, BG_CMA);
		  if(state=='M'){
		    if(col_freq[col] != NULL) NullFreqBG[j]=col_freq[col];
		    WtFrqBG[j] = TmpWtFreq[col];
		    ObsBG[j] = TmpObs[col];
		    FractSeqAlnBG[j]=1.0;
		  }
		 break;
		case 'm':	// match within motif.
		  j++; col++;
		  {
		    if(col_freq[col] != NULL) NullFreqBG[j]=col_freq[col];
		    WtFrqBG[j] = TmpWtFreq[col];
		    ObsBG[j] = TmpObs[col];
		    FractSeqAlnBG[j]=1.0;
		  }
		 break;
		case 'd':	// deletion within motif...
		 col++;
		 break;
		case 'I':	// insert within motif...
		  j++;
		 break;
		case 'i':	// insert between motif blocks...
		  j++;
		 break;
		case 'E':	// end of operational string...
		 assert(j==length);
		 break;
		default: print_error("bug in operation string!"); break;
	    }
	  }
	} else {
	 for(Int4 m=1; m <= nBlksCMSA(BG_CMA); m++){
	  TmpWtFreq=ColResFreqsCMSA(m, &TmpObs, BG_CMA);
	  lenM = SiteLen(m,S);
	  site = mtfpos[m];
	  end = site+lenM-1;
	  if(m==1) StartBG_CMA=site;
	  else if(m == nBlksCMSA(BG_CMA)) EndBG_CMA=end;
	  assert(site > 0 && end <= length);
	  // WARNING: I'm assuming that there are no indels in sequence 1 of BG_CMA!!!
	  col_freq=ColResFreqsCMSA(m, BG_CMA);
	  for(j=site,col=1; j <= end && col <= lenM; j++,col++){
		  if(col_freq[col] != NULL){
		   if(extra_files) fprintf(stderr,"%c%d (%.3f->%.3f)\n",
				  AlphaChar(ResSeq(j,keyE),AB), j,
				  NullFreqBG[col][ResSeq(j,keyE)],
				  col_freq[col][ResSeq(j,keyE)]);
		  }
		  if(col_freq[col] != NULL){
			  NullFreqBG[j]=col_freq[col];
		  }
		  WtFrqBG[j] = TmpWtFreq[col];
		  ObsBG[j] = TmpObs[col];
		  FractSeqAlnBG[j]=1.0;
	  } free(col_freq); free(TmpObs); free(TmpWtFreq);
	  fprintf(stderr,"\n");
	 }
	}
	fprintf(stderr,"Background frequency cma file is okay\n");
}

