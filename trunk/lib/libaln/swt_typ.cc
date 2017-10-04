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

#include "swt_typ.h"

#define TRY_SPEED_UP	1

swt_typ::swt_typ(cma_typ orthologs,cma_typ paralogs,BooLean mode, BooLean UseGlobalSqWts)
{ Init(orthologs,paralogs,mode,0,UseGlobalSqWts,0); }
// { Init(orthologs,paralogs,mode,(e_type)0,UseGlobalSqWts,(char*)0); }

swt_typ::swt_typ(cma_typ orthologs,cma_typ paralogs,BooLean mode,
	BooLean UseGlobalSqWts,char *wts_file_name)
{ Init(orthologs,paralogs,mode,0,UseGlobalSqWts,wts_file_name); }

swt_typ::swt_typ(cma_typ orthologs,cma_typ paralogs,BooLean mode,
	e_type cons_seq,BooLean UseGlobalSqWts)
{ Init(orthologs,paralogs,mode,cons_seq,UseGlobalSqWts,0); }

swt_typ::swt_typ(cma_typ orthologs,cma_typ paralogs,BooLean mode,
	e_type cons_seq,BooLean UseGlobalSqWts,char *wts_file_name)
{ Init(orthologs,paralogs,mode,cons_seq,UseGlobalSqWts,wts_file_name); }

swt_typ::swt_typ(cma_typ paralogs,BooLean UseGlobalSqWts)
{ Init(0,paralogs,0,0,UseGlobalSqWts,0); }

void	swt_typ::Init(hsw_typ HSW)
{
	assert(HSW); rtn_hsw=0; hsw=HSW; own_hsw=FALSE; WtsFileName=0; ocma=0; 
	Verbose=FALSE;
        Weight = hsw->Weight;       // position specific sequence weights
        WtCnts=hsw->WtCnts;       // Weighted counts at each position.
        WtFreq=hsw->WtFreq;
	Length=hsw->Length;
        NWtSq=hsw->NWtSq;
        mcma=hsw->mcma;
        AB=hsw->AB;
        fract_seq_aln=hsw->fract_seq_aln;
	compute_wts=FALSE;
}

void	swt_typ::Init(cma_typ orthologs,cma_typ paralogs,BooLean mode,
	e_type cons_seq,BooLean UseGlobalSqWts,char *wts_file_name)
{
	rtn_hsw=0; hsw=0; own_hsw=TRUE;
	if(UseGlobalSqWts) StartAlpha=0; else StartAlpha=1;
	if(wts_file_name) WtsFileName=wts_file_name; else WtsFileName=0;
	ocma=orthologs; use_ocma_pseudo=mode;
	if(ocma){
	   assert(nBlksCMSA(ocma) == 1); 
	   if(cons_seq){
	     assert(LengthCMSA(1,ocma) == LenSeq(cons_seq));
	     csq=CopySeq(cons_seq); Length=LenSeq(csq);
	   } else { csq=MkConsensusCMSA(ocma); Length=LenSeq(csq); }
	} else { csq=0; Length=LengthCMSA(1,paralogs); }
	if(WtsFileName){
		print_error("WtsFileName option for swt_typ not yet fully implemented");
		FILE *fp = open_file(WtsFileName,".swt","r");
		FReadHenikoffWeights(fp); fclose(fp);
		compute_wts=FALSE;
	} else { init(paralogs); compute_wts=TRUE; }
}

void	swt_typ::init(cma_typ paralogs)
{
	Int4	s,r,sq;
	mcma=paralogs; assert(nBlksCMSA(mcma) == 1);
	AB=AlphabetCMSA(mcma);
	Verbose=TRUE;
#if 0  // DEBUG.
WriteCMSA("junkZ.cma",mcma);
PutSeq(stdout,csq,AB);
#endif
	if(LengthCMSA(1,mcma) != Length) print_error("swt input error");
	assert(Length > 10);
	// rm_high_info=(Length/10)+1; 	// remove 1/10th of high info columns
	rm_high_info=(Length/20)+1; 	// remove 1/20th of high info columns
	keep_many_partic=(Length/5)+1;  // retain 1/5th of rest of columns
	// keep_many_partic=(Length/4)+1;  // retain 1/4th of rest of columns
	NWtSq=NumSeqsCMSA(mcma);
	NEW(Participant,Length+2,set_typ);
	NEW(IS,Length+2,set_typ);
	NEWP(RawCnts,Length+2,set_typ);
	NEW(nResTyp,Length+2,unsigned char);
	NEWP(Observed,Length+2,Int4);

	if(hsw==0){
	  NEWP(Weight,Length+2,double); NEWP(WtCnts,Length+2,double);
	  NEWP(WtFreq,Length+2,double); NEW(fract_seq_aln,Length+2,double);
	}

	MPWtFreq=0;
	MargProb=0;
	for(s=1; s <= Length; s++){

	     if(hsw==0){
		NEW(WtCnts[s],nAlpha(AB)+2,double); NEW(WtFreq[s],nAlpha(AB)+2,double);
		NEW(Weight[s],NWtSq+2,double);
	     }

		Participant[s]=MakeSet(NWtSq+1); // creates empty set.
		IS[s]=MakeSet(NWtSq+1); // creates empty set.
		NEW(RawCnts[s],nAlpha(AB)+2,set_typ);
		for(r=StartAlpha; r <= nAlpha(AB); r++) RawCnts[s][r]=MakeSet(NWtSq+1);
		// if(StartAlpha==1) RawCnts[s][0]=MakeSet(NWtSq+1);
		// ^ trying to fix memory fault...
		NEW(Observed[s],nAlpha(AB)+2,Int4);
	}
	is=MakeSet(NWtSq+1); // creates empty set.
	Int4 start;
	if(ocma) { start=2; ofreq=ColResFreqsCMSA(1,0,&oObs,ocma); }
	else { start=1; ofreq=0; oObs=0; }
#if TRY_SPEED_UP   // Speed up AFN(12/1/01)...
    {
	NEWP(PtrWtSq,NWtSq+3,unsigned char);	// pointer to fakeSeqs
	Int4 s0;
	for(sq=1; sq <= NWtSq; sq++){
           Int4 s0; s0=SitePos(1,sq,1,SitesCMSA(mcma)); // start of fake sequence.
           unsigned char *fsq = SeqSeqSet(sq,DataCMSA(mcma));	// fake seq
	   PtrWtSq[sq]=fsq+s0-1;
	   if(sq >= start){
	     for(s=1; s <= Length; s++,s0++){
		r=fsq[s0];
		if(r != 0){ AddSet(sq,Participant[s]); AddSet(sq,RawCnts[s][r]); }
	     }
	   }
	}
    }
#else
	for(s=1; s <= Length; s++){
	   for(sq=start; sq <= NWtSq; sq++){
		r=ResidueCMSA(1,sq,s,mcma);
		if(r != 0){ AddSet(sq,Participant[s]); AddSet(sq,RawCnts[s][r]); }
	   }
	}
#endif
	NEW(numParticipants,Length+2,Int4);
	NEW(info,Length+2,double); NEW(order,Length+2,Int4);
	NEW(info_order,Length+2,Int4);
	compute_mp=TRUE;
	freq=tFreqSeqSet(TrueDataCMSA(mcma));
	double I;
	UpperLimitInfo=0.0;
        for(r=StartAlpha; r <= nAlpha(AB); r++){
            double q = freq[r];
            if(q > 0.0) I = log(1.0/q); else I=0.0;
	    if(UpperLimitInfo < I) UpperLimitInfo=I;
        } assert(UpperLimitInfo > 0.0);
}

void	swt_typ::SumSeqWts(BooLean *skip, Int4 S,double minfract,Int4 to_use)
// b. Sum sequence weights: separated out for speed...
// This doesn't seem to speed up the program...and chnages weighting.  Need to check.
{
	register Int4	s,r,*obs_s;
	register double fraction,*wt_frq_s,sum;
	register double expected_wt=0.0;
	register unsigned char   nResTyp_s;
	double		nCol=0.0,numused=0;

	print_error("SumSeqWts( ) has serious problems that need to be fixed!");
        for(Int4 i=Length; i > 0; i--){
	     s=order[i];
	     if(skip[s]) continue; //  skip over highly conserved columns
	     nResTyp_s=nResTyp[s];
	     if(nResTyp_s <= 1) continue; 
	     fraction=(double)numParticipants[s]/(double)numParticipants[S];
	     if(fraction < minfract && numused > to_use) break;
	     numused++;
	     obs_s=Observed[s]; wt_frq_s=WtFreq[s];
	     if(StartAlpha == 1){
	       for(sum=0.0,r=nAlpha(AB); r >= 1; r--) if(obs_s[r] > 0) sum+=wt_frq_s[r];
	       for(r=nAlpha(AB); r >= 1; r--){
	         if(obs_s[r] > 0){
	           expected_wt+= fraction*(wt_frq_s[r]/sum)*(1./(double)(obs_s[r]*nResTyp[s]));
	           // expected_wt+= fraction*(wt_frq_s[r]/sum)*(double)nResTyp[s]/(double)(Observed[s][r]);
	         }
	       }
	     } else {	// StartAlpha==0;
	       for(sum=0.0,r=nAlpha(AB); r; r--) if(obs_s[r] > 0) sum+=wt_frq_s[r];
	       for(r=nAlpha(AB); r; r--){
	         if(obs_s[r] > 0){
	           expected_wt+= fraction*(wt_frq_s[r]/sum)*(1./(double)(obs_s[r]*nResTyp[s]));
	           // expected_wt+= fraction*(wt_frq_s[r]/sum)*(double)nResTyp[s]/(double)(Observed[s][r]);
	         }
	       }
	     }
	     nCol+=fraction;
	} 
	if(nCol > 0) Weight[S][0]=expected_wt/nCol; // Raw expected (avg) Weight.
	Weight[S][0] *= 2.0;  // assume ave wt == 0.5 maximum weight because
}

void	swt_typ::HenikoffWeights( )
{
        Int4    i,n,s,S,r;
	Int4	numused,to_use=keep_many_partic-rm_high_info;
        double  total,sum,minfract,fraction,*nCol,*pseudo;
	BooLean *skip;
	time_t	time1=time(NULL);

to_use = 20;
rm_high_info=10;
minfract=0.75;

	NEW(pseudo,nAlpha(AB)+2,double);
	NEW(nCol,NWtSq+2,double); NEW(skip,Length+2,BooLean);
      for(S=1;S<=Length; S++){
	// 1. Get participant information:
	SortColumns(S);
	fract_seq_aln[S] = (double)numParticipants[S]/(double)NWtSq;
	if(ocma && use_ocma_pseudo) {	// base pseudocounts on ocma alignment
		for(r=0; r<=nAlpha(AB); r++){
		   pseudo[r]=ofreq[S][r]/(double)(numParticipants[S] + 1.0);
		   // ADD ONE count for each residue...
		   // pseudo[r]=ofreq[S][r]/(double)NWtSq;
		   // pseudo[r]=(double)nAlpha(AB)*ofreq[S][r]/(double)NWtSq;
		}
	} else {
		for(r=0; r<=nAlpha(AB); r++) pseudo[r]=1.0/(double)NWtSq; 
		// due to ave wt=0.5:
		//for(r=0; r<=nAlpha(AB); r++) pseudo[r]=2.0/(double)NWtSq;
		// Ave wt=0.5 and set pseudo to 2 r residues for robustness --> 
		//for(r=0; r<=nAlpha(AB); r++) pseudo[r]=4.0/(double)NWtSq;
		// for(r=0; r<=nAlpha(AB); r++) pseudo[r]=(double)nAlpha(AB)/((double)NWtSq);
		// ADD ONE/NumberSeqs pseudocounts for each residue...
		// for(r=0; r<=nAlpha(AB); r++) pseudo[r]=freq[r]/(double)NWtSq; 
		// for(r=0; r<=nAlpha(AB); r++) pseudo[r]=2.0/((double)numParticipants[S]);
		// for(r=0; r<=nAlpha(AB); r++) pseudo[r]=1.0/((double)numParticipants[S]);
	}
	// 2. Initialize data structures:
	for(r=0; r <=nAlpha(AB); r++){ WtFreq[S][r]=pseudo[r]; }
        for(n=0;n<=NWtSq;n++){ Weight[S][n]=0.0; nCol[n]=0.0; }
	numused=0;
	if(numParticipants[S] == 0){		// case where no aligned sequences.
#if 0
	} else if(!ocma && numParticipants[S] == 1){	// case where only one sequence.
	  if(MemberSet(1,IS[S])){		// assume query sequence aligned
		Weight[S][1] = 1.0; 
		r = ResidueCMSA(1,1,S,mcma); 
		WtFreq[S][r] += 1.0; 
	  } else assert(MemberSet(1,IS[S]));
#else
	} else if(!ocma && numParticipants[S] == 1 && MemberSet(1,IS[S])){
		Weight[S][1] = 1.0; 
		r = ResidueCMSA(1,1,S,mcma); 
		WtFreq[S][r] += 1.0; 
#endif
	} else {
	  // 3. Skip over certain columns:
	  for(i=1; i<=Length; i++) { skip[i]=FALSE;  }
	  for(i=1; i<=rm_high_info; i++) { s=info_order[i]; skip[s]=TRUE; }
	  // 4. Sum sequence weights:
          for(i=1; i<=Length; i++){
	     s=order[i];
	     if(nResTyp[s] <= 1) continue; 
	     if(skip[s]) continue; //  skip over highly conserved columns
	     fraction=(double)numParticipants[s]/(double)numParticipants[S];
assert(fraction <= 1.0);
	     if(fraction < minfract && numused > to_use) break;
	     numused++; 
             for(n=1;n<=NWtSq;n++){
		if(!MemberSet(n,IS[s])) continue;
#if TRY_SPEED_UP   // Speed up AFN(12/1/01)...
		r=PtrWtSq[n][s];
#else
		r=ResidueCMSA(1,n,s,mcma);
#endif
		Weight[S][n] += fraction/(double)(Observed[s][r]*nResTyp[s]);
		// Weight[S][n] += 1./(double)(Observed[s][r]*nResTyp[s]);
                // Weight[S][n] += (double)nResTyp[s]/(double)(Observed[s][r]);
		nCol[n]+=fraction;
		// nCol[n]+=1.0;
             }
          }
	  // 5. Compute as avg. weights over number of aligned cols.
	  for(n=1;n<=NWtSq;n++) { 
	   if(nCol[n] > 0.0){ Weight[S][n] = Weight[S][n]/nCol[n]; }
	  }
	  // 6. Normalize weights and compute and normalize WtFreq...
          for(total=0.0,n=1; n <= NWtSq; n++) total+=Weight[S][n];
#if 0	//**************** Fix problem with pseudocounts here. *********************
	  // Now have actual total number of weighted sequences at each position....
	  // ave wt = total/NWtSq;
	  double add_to_total=0.0;
	  for(r=0; r <=nAlpha(AB); r++){
		WtFreq[S][r]=(total/NWtSq); 
		add_to_total+=(total/NWtSq); 
		
	  } total += add_to_total;
	  // add ave. weighted sequence as pseudocount for each residue type.
	  for(r=0; r <=nAlpha(AB); r++){
		WtFreq[S][r]=WtFreq[S][r]/total;
	  }
#endif	//********************************************************
	  if(total > 0.0) for(n=1; n <= NWtSq; n++){
#if TRY_SPEED_UP   // Speed up AFN(12/1/01)...
		r=PtrWtSq[n][S];
#else
                r = ResidueCMSA(1,n,S,mcma);
#endif
                if(r>0) WtFreq[S][r] += Weight[S][n]/total;
		// NOTE: WtFreq was initialized above.
          }
	}
	// 6. Normalize WtFreq...
        for(sum=0.0,r=StartAlpha; r <=nAlpha(AB); r++){ sum+=WtFreq[S][r]; }
assert(sum > 0.0);
        for(r=StartAlpha; r <=nAlpha(AB); r++){ WtFreq[S][r] /= sum; }
      } 
      // 8. COMPUTE EFFECTIVE NUMBER OF SEQUENCES HERE...
      // 8. Compute expected weight based on weighted frequencies.
      for(S=1; S <= Length; S++){
	  SortColumns(S);
	  // if(numParticipants[S] == 0) ;
	  // a. Skip over certain columns:
	  for(i=1; i<=Length; i++) { skip[i]=FALSE; }
	  for(i=1; i<=rm_high_info; i++) { s=info_order[i]; skip[s]=TRUE; }
#if 1
	  // b. Sum sequence weights:
          double expected_wt=0.0; numused=0; nCol[0]=0.0;
          for(i=1; i<=Length; i++){
	     s=order[i];
	     if(nResTyp[s] <= 1) continue; 
	     if(skip[s]) continue; //  skip over highly conserved columns
	     fraction=(double)numParticipants[s]/(double)numParticipants[S];
	     if(fraction < minfract && numused > to_use) break;
	     numused++;
	     for(sum=0.0,r=StartAlpha; r<=nAlpha(AB); r++) if(Observed[s][r]>0) sum+=WtFreq[s][r];
	     for(r=StartAlpha; r <=nAlpha(AB); r++){
	       if(Observed[s][r] > 0){
	         double f=WtFreq[s][r]/sum;
	         expected_wt+= fraction*f*(1./(double)(Observed[s][r]*nResTyp[s]));
	         // expected_wt+= fraction*f*(double)nResTyp[s]/(double)(Observed[s][r]);
	       }
	     } nCol[0]+=fraction;
	     // nCol[0]+=1.0;
	  } 
	  if(nCol[0] > 0) Weight[S][0]=expected_wt/nCol[0]; // Raw expected (avg) Weight.
	  Weight[S][0] *= 2.0;  // assume ave wt == 0.5 maximum weight because
				// can't be < 0.0 and assume normal distribution.
#else
	  SumSeqWts(skip,S,minfract,to_use);
#endif
          // 10. Compute weighted observed counts...
	  for(n=1; n <= NWtSq; n++){
#if 1	// WARNING: ATTEMPTING TO FIX A SERIOUS WEIGHTING BUG.
		// assert(Weight[S][0] > 0 && Weight[S][n] < 10000);
		if(Weight[S][0] <= 0){ Weight[S][n] = 0; continue; }
#endif
		Weight[S][n] = Weight[S][n]/Weight[S][0]; 
		if(Weight[S][n] > 1.0) Weight[S][n]=1.0;  // truncate weight at 1.0
	  }
      } // END: COMPUTE EFFECTIVE NUMBER OF SEQUENCES.
      free(skip); free(nCol);
      h_type HG = Histogram("Weighted Observed Counts",0,NWtSq,2.0);
      for(s=1; s <= Length; s++){
          for(r=0; r<=nAlpha(AB); r++) WtCnts[s][r]=0.0;
          for(n=1; n<=NWtSq; n++){
#if TRY_SPEED_UP   // Speed up AFN(12/1/01)...
		r=PtrWtSq[n][s];
#else
                r = ResidueCMSA(1,n,s,mcma);
#endif
                if(r > 0) WtCnts[s][r] += Weight[s][n]; // position specific weight
#if 0	// correct counts at gap positions...
		else {
		  // fprintf(stderr,"adjusting for gaps...\n");
		  for(Int4 r0=StartAlpha; r0 <= nAlpha(AB); r0++){
			WtCnts[s][r0] += freq[r0]*Weight[s][n]; // add partial counts.
		  }
		}
#endif
          }
          // for(r=0; r<=nAlpha(AB); r++)
          for(r=StartAlpha; r<=nAlpha(AB); r++)
	  {
		if(WtCnts[s][r]>0.0) IncdHist(100*WtCnts[s][r]/CardSet(RawCnts[s][r]),HG); 
	  } 
#if 1	//**************** Fix problem with pseudocounts here. *********************
	  // Use 1 pseudocount for weights...
	  double total_wt_cnts=0.0;
          for(r=StartAlpha; r<=nAlpha(AB); r++){ total_wt_cnts+=WtCnts[s][r]+0.01; }
	  for(r=StartAlpha; r<=nAlpha(AB); r++){ 
		WtFreq[s][r] = (WtCnts[s][r]+0.01)/total_wt_cnts;
	  }
#else
#endif	//********************************************************
      } if(Verbose) PutHist(stderr,60,HG); 
	NilHist(HG);
      compute_wts=FALSE; free(pseudo);
      double runtime=difftime(time(NULL),time1);
      if(Verbose) fprintf(stderr,"\ttime: %0.1f seconds (%0.2f minutes)\n",runtime,runtime/60.0);
}

double  swt_typ::relative_entropy(Int4 s)
{
        double  I,total=(double)numParticipants[s];
        Int4    r,*observed=Observed[s];

        for(I=0.0,r=StartAlpha; r <= nAlpha(AB); r++){
            double p = (double)observed[r]/total;
            double q = freq[r];
            if(p > 0.0 && q > 0.0) I += p*log(p/q);
        } return I;
}

BooLean	swt_typ::SortColumns(Int4 S)
// Put the number of participants that position S shares with other positions 
{
	Int4	r,s,i,total;
	keytyp	key;

	if(S < 1 || S > Length) return FALSE;
	minPart=NWtSq; maxPart=0;
	mininfo=DBL_MAX; maxinfo=0.0;
	dh_type	dH=dheap(Length+2,4);
	dh_type	iH=dheap(Length+2,4);
	for(i=0,s=1; s <= Length; s++,i++){
	  IntersectSet1(Participant[S],Participant[s],IS[s]);
	  nResTyp[s]=0;
	  for(total=0,r=StartAlpha; r <= nAlpha(AB); r++){
		Int4 obs=CardInterSet(IS[s],RawCnts[s][r]);
		if(obs > 0) nResTyp[s]++;
		Observed[s][r]=obs; total += obs;
	  } numParticipants[s]=total;
	  if(minPart > total) minPart = total;
	  if(maxPart < total) maxPart = total;
	  info[s]=relative_entropy(s);
	  key = (keytyp)-numParticipants[s] -1. + (info[s]/UpperLimitInfo);
	  // more informative positions below less informative ones.
	  insrtHeap(s,key,dH);
	  if(mininfo > info[s]) mininfo=info[s];
	  if(maxinfo < info[s]) maxinfo=info[s];
	  insrtHeap(s,(keytyp)info[s],iH);
	}
	for(i=1; (s=delminHeap(dH)) != NULL; i++) order[i]=s;
	for(i=1; (s=delminHeap(iH)) != NULL; i++) info_order[i]=s;
	Nildheap(dH); Nildheap(iH);
	return TRUE;
}

double  **LocalMarginalProb(e_type sE, char *operation, Int4 Start, Int4 Oper_len, 
        double *score, Int4 length, double **WtFreq, double *fract_seq_aln, double *bfreq, 
		Int4 PerNats, Int4 gapopen, Int4 gapextend, a_type AB, unsigned char StartAlpha) 
{
        double          **TMAT;
        Int4            i,im1,j,jm1,seq_len=LenSeq(sE);
        unsigned char   *seq=SeqPtr(sE);
        e_type          rE=CopySeq(sE);


        rE = ReverseSeq(rE);
        unsigned char   *rseq=XSeqPtr(rE);
        double **SMAT,**BMAT,**SDEL,**BDEL,**SINS,**BINS;
        NEWP(SMAT,length+2,double); NEWP(BMAT,length+2,double); NEWP(TMAT,length+2,double);
        NEWP(SDEL,length+2,double); NEWP(BDEL,length+2,double);
        NEWP(SINS,length+2,double); NEWP(BINS,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(SMAT[i],seq_len+1,double); NEW(BMAT[i],seq_len+1,double); 
                NEW(TMAT[i],seq_len+1,double);
                NEW(SDEL[i],seq_len+1,double); NEW(BDEL[i],seq_len+1,double);
                NEW(SINS[i],seq_len+1,double); NEW(BINS[i],seq_len+1,double);
        }

        double m2d,m2i,d2d,i2i;

        m2d = m2i = exp((double) -gapopen/(double) PerNats);
        d2d = i2i = exp((double) -gapextend/(double) PerNats);

        double **pmat_emit_prob,**rmat_emit_prob;

        NEWP(pmat_emit_prob,length+2,double);
        NEWP(rmat_emit_prob,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(pmat_emit_prob[i],nAlpha(AB)+2,double);
        }
        double *pmat_emit_prob_i,*WtFreq_i;

        for(i=1;i<=length;i++){
                pmat_emit_prob_i = pmat_emit_prob[i];
		WtFreq_i = WtFreq[i];
		if(fract_seq_aln[i] > 0.0) {
                	for(j=StartAlpha;j<=nAlpha(AB);j++){
                        	pmat_emit_prob_i[j] = WtFreq_i[j]/bfreq[j];
                	}
		} else {
			for(j=StartAlpha;j<=nAlpha(AB);j++){
				pmat_emit_prob_i[j] = 1.0;
			}
		}
		pmat_emit_prob_i[0]=0.5; // set 'X' residues to this odds score.
        }

        for(i=1;i<=length;i++){ rmat_emit_prob[i] = pmat_emit_prob[length-i+1]; }

        SDEL[1][1] = 0; SMAT[1][1] = pmat_emit_prob[1][seq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                SINS[j][1] = 0; SDEL[j][1] = 0;
                SMAT[j][1] = pmat_emit_prob[j][seq[1]];
        }
        SINS[1][1] = 0; SMAT[1][2] = pmat_emit_prob[1][seq[2]];
        SINS[1][2] = pmat_emit_prob[1][seq[1]] * m2i; 
        SDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                SDEL[1][i] = 0; SMAT[1][i] = pmat_emit_prob[1][seq[i]];
                SINS[1][i] = (SINS[1][im1] * i2i) + (SMAT[1][im1] * m2i); 
        }

        BDEL[1][1] = 0; BMAT[1][1] = rmat_emit_prob[1][rseq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                BINS[j][1] = 0; BDEL[j][1] = 0;
                BMAT[j][1] = rmat_emit_prob[j][rseq[1]];
        }                        
        BINS[1][1] = 0; BMAT[1][2] = rmat_emit_prob[1][rseq[2]];
        BINS[1][2] = rmat_emit_prob[1][rseq[1]] * m2i;
        BDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                BDEL[1][i] = 0; BMAT[1][i] = rmat_emit_prob[1][rseq[i]];
                BINS[1][i] = (BINS[1][im1] * i2i) + (BMAT[1][im1] * m2i);
        }

        double *smatj,*sdelj,*sinsj,*bmatj,*bdelj,*binsj;
        double *smatjm1,*sdeljm1,*sinsjm1,*bmatjm1,*bdeljm1,*binsjm1;
        double *pmat_emit_probj,*rmat_emit_probj,pmat_emit_probjres,rmat_emit_probjres;
        unsigned char res,rres;
        for(jm1=1,j=2;j<=length;j++,jm1++){
                smatj=SMAT[j];sdelj=SDEL[j];sinsj=SINS[j];
                bmatj=BMAT[j];bdelj=BDEL[j];binsj=BINS[j];
                smatjm1=SMAT[jm1];sdeljm1=SDEL[jm1];sinsjm1=SINS[jm1];
                bmatjm1=BMAT[jm1];bdeljm1=BDEL[jm1];binsjm1=BINS[jm1];
                pmat_emit_probj=pmat_emit_prob[j];
                rmat_emit_probj=rmat_emit_prob[j];
                for(im1=1,i=2;i<=seq_len;im1++,i++){
                        res=seq[i]; rres=rseq[i];
                        pmat_emit_probjres=pmat_emit_probj[res]; 
                        rmat_emit_probjres=rmat_emit_probj[rres];
                        smatj[i] += (smatjm1[im1] * pmat_emit_probjres);
                        smatj[i] += (sdeljm1[i] * pmat_emit_probjres); 
                        smatj[i] += (sinsjm1[im1] * pmat_emit_probjres);
			smatj[i] += pmat_emit_probjres;
                        sdelj[i] += (smatjm1[im1] * m2d);
                        sdelj[i] += (sdeljm1[i] * d2d);
                        sinsj[i] += (smatj[im1] * m2i);
                        sinsj[i] += (sinsj[im1] * i2i);
                        bmatj[i] += (bmatjm1[im1] * rmat_emit_probjres);
                        bmatj[i] += (bdeljm1[i] * rmat_emit_probjres);
                        bmatj[i] += (binsjm1[im1] * rmat_emit_probjres);
			bmatj[i] += rmat_emit_probjres;
                        bdelj[i] += (bmatjm1[im1] * m2d);
                        bdelj[i] += (bdeljm1[i] * d2d);
                        binsj[i] += (bmatj[im1] * m2i);
                        binsj[i] += (binsj[im1] * i2i);
                }
        }
        double *TMAT_j,*SMAT_j,*BMAT_l,*pmat_emit_prob_j;
        for(j=1;j<=length;j++){
                TMAT_j = TMAT[j]; SMAT_j = SMAT[j]; BMAT_l = BMAT[length-j+1];
                pmat_emit_prob_j = pmat_emit_prob[j];
                for(i=1;i<=seq_len;i++){
                        res=seq[i];
                        TMAT_j[i]=BMAT_l[seq_len-i+1]*(SMAT_j[i]/pmat_emit_prob_j[res]);
                }
        }

	double all = 0.;
	for(j=1;j<=length;j++){
		smatj=SMAT[j];
		for(i=1;i<=seq_len;i++){	
			all += smatj[i]; 
		}
	}

        for(j=1;j<=length;j++){
                TMAT_j = TMAT[j];
                for(i=1;i<=seq_len;i++){
                        TMAT_j[i] /= all;
//printf("TMAT[%d][%d] = %lf\n",j,i,TMAT_j[i]);
                }
        }

        if(operation !=NULL){
                NEW(score,Oper_len,double);
                i=1;j=1; 
                Int4 posProf=1, posSeq=Start;
                while(operation[i] != 'E'){
                        switch(operation[i++]){
                                case 'M':
                                case 'm': score[j++] = TMAT[posProf++][posSeq++]; break; 
                                case 'D': 
                                case 'd': j++; posProf++; break;
                                case 'I':
                                case 'i': j++; posSeq++; break; 
                                default : print_error("MarginalProb(): error in the operation array");
                        }
                }
        }

        for(i=0;i<=length+1;i++) { 
                free(SMAT[i]);free(SDEL[i]);free(SINS[i]);
                free(BMAT[i]);free(BDEL[i]);free(BINS[i]);
        }
        free(SMAT);free(SDEL);free(SINS);free(BMAT);free(BDEL);free(BINS);
        for(i=0;i<=length+1;i++) { free(pmat_emit_prob[i]); } 
        free(pmat_emit_prob);free(rmat_emit_prob); NilSeq(rE);
        if(!(all > 0.0)){
          for(i=0;i<=length+1;i++) free(TMAT[i]); free(TMAT);
          return 0;
        } free(TMAT[0]); TMAT[0]=0; // AFN: fixes memory leak.
        free(TMAT[length+1]); TMAT[length+1]=0; // AFN: fixes memory leak.
        return TMAT;
}

void	swt_typ::ComputeMargProb(Int4 pernats,Int4 gap_open, Int4 gap_extend)
{
	Int4		s;
        unsigned char	r;
	BooLean		success;

	assert(csq);
	if(!compute_mp) return;
	if(compute_wts) HenikoffWeights( ); compute_mp=FALSE;
	NEWP(MPWtFreq,Length+2,double);
	NEW(MargProb,Length+2,double);
	for(s=1; s <= Length; s++){
		NEW(MPWtFreq[s],nAlpha(AB)+2,double);
#if 0
		for(r=StartAlpha; r <= nAlpha(AB); r++){
		    MPWtFreq[s][r]=WtFreq[s][r];
		}
#endif
	}
        {
           Int4                 i,j;
	   // ss_type	data=TrueDataCMSA(ocma);
	   // e_type qE = SeqSeqSet(1,data);
	   // assert(LenSeq(csq) == Length);
	   // PutSeq(stdout,csq,AB); 
	   double  **mp=LocalMarginalProb(csq,NULL,0,0,NULL,
        		Length,WtFreq,fract_seq_aln,freq,
			pernats, gap_open,gap_extend,AB,StartAlpha);
	   if(mp == 0) success=FALSE; else success=TRUE;
// IncdHist(1.0,HG); IncdHist(0.0,HG); 
           h_type HG = Histogram("Marg Prob Weighted Freq",0,1,0.02);
           if(success){
	    for(j=1; j <= LenSeq(csq); j++){	// query sequence
	     MargProb[j]=mp[j][j];
             for(i=1; i <= Length; i++){		// main profile
		if(mp[i][j] >= 0.001) IncdHist(mp[i][j],HG); 
		if(mp[i][j] < 0.0  || mp[i][j] > 1.0) success=FALSE;
		for(r=StartAlpha; r <= nAlpha(AB); r++){
		   MPWtFreq[j][r] += mp[i][j]*WtFreq[i][r];
// assert(MPWtFreq[j][r] <= 1.0);
		   // if(i==j) MPWtFreq[j][r] += WtFreq[i][r];
		}
		if(!success) break;
		// sum probability that residue at j is aligned with profile at i.
             }
	     if(!success) break;
	    }
           }
	   if(success){
	     for(i=1; i <= LenSeq(csq); i++){
		double total=0.0;
		for(r=0; r <= nAlpha(AB); r++) total += MPWtFreq[i][r];
		for(r=0; r <= nAlpha(AB); r++){ MPWtFreq[i][r] /= total; }
	     }
	   }
	   if(mp){ // NEW: DEALLOCATE MEMORY.
	     for(i=1; i <= LenSeq(csq); i++) free(mp[i]); free(mp);
	   }
	   if(!success){ 
		for(s=1; s <= Length; s++) free(MPWtFreq[s]); 
		free(MPWtFreq); MPWtFreq=0;
	   } else if(Verbose) PutHist(stderr,60,HG); NilHist(HG);
        }
}

void	swt_typ::Test(FILE *fp)
{
	if(!own_hsw) print_error("swt_typ: hsw obtained from elsewhere");
	Int4	s,time1;
	time1=time(NULL);
	for(s=1; s <= Length; s++){
	     SortColumns(s);
	  // numParticipants[s]=CardSet(IS[s]);
	} 
	fprintf(fp,"\ttime: %d seconds (%0.2f minutes)\n",
                        time(NULL)-time1,(float)(time(NULL)-time1)/60.0);
}

void	swt_typ::Put(FILE *fp,Int4 S)
// Put the number of participants that position S shares with other positions 
{
	Int4	r,s,i,j;

	if(S < 1 || S > Length) return;
	if(!own_hsw) print_error("swt_typ: hsw obtained from elsewhere");
	SortColumns(S);
	fprintf(fp,"Pos %d:\n",S);
	for(i=0,j=1; j <= Length; j++,i++){
	  if(i % 50 == 0 && j < Length){
	    fprintf(fp,"\npos (cnts): ");
	    for(r=StartAlpha; r <= nAlpha(AB); r++) fprintf(fp,"   %c",AlphaChar(r,AB));
	    fprintf(fp,"  nats typs\n");
	  } 
	  s=order[j];
	  // s=info_order[j];
	  fprintf(fp,"%3d (%4d): ",s,numParticipants[s]);
	  for(r=StartAlpha; r <= nAlpha(AB); r++) fprintf(fp,"%4d",Observed[s][r]);
	  fprintf(fp,"  %3.1f %3d\n",info[s],nResTyp[s]);
	} fprintf(fp,"\n");
}

void	swt_typ::Put(FILE *fp)
{
	if(!own_hsw) print_error("swt_typ: hsw obtained from elsewhere");
	Int4	r,s,i;
	FillSet(is); 
	for(i=0,s=1; s <= Length; s++,i++){
	  if(i % 50 == 0 && s < Length){
	    fprintf(fp,"\npos (cnts): ");
	    for(r=StartAlpha; r <= nAlpha(AB); r++) fprintf(fp,"   %c",AlphaChar(r,AB));
	    fprintf(fp,"\n");
	  }
	  IntersectSet3(is,Participant[s]);
	  Int4 card=CardSet(Participant[s]);
	  fprintf(fp,"%3d (%4d): ",s,card);
	  for(r=StartAlpha; r <= nAlpha(AB); r++) fprintf(fp,"%4d",CardSet(RawCnts[s][r]));
	  fprintf(fp,"\n");
	} fprintf(fp,"\n");
	fprintf(fp,"\n%d sequences span the entire alignment\n",CardSet(is));
}

void	swt_typ::Free( )
{
   if(own_hsw){		// else don't use anything...
	for(Int4 s=1; s <= Length; s++){
	   if(MPWtFreq) free(MPWtFreq[s]);
	   if(ocma){ free(ofreq[s]); free(oObs[s]); }
	   NilSet(Participant[s]); NilSet(IS[s]);
// fprintf(stderr,"swt->Free 1d(%d)\n",s);
	   // Memory corruption occurring here when 
	   for(Int4 r=StartAlpha; r <= nAlpha(AB); r++){
// fprintf(stderr,"swt->Free 1d(%c%d = %d)\n",AlphaChar(r,AB),s,CardSet(RawCnts[s][r]));
		NilSet(RawCnts[s][r]);
	   }
	   free(RawCnts[s]); free(Observed[s]);
	   // if(own_hsw){ free(Weight[s]); free(WtCnts[s]); free(WtFreq[s]); }
	   free(Weight[s]); free(WtCnts[s]); free(WtFreq[s]);
	} free(RawCnts); free(Participant); free(IS); free(info);
// fprintf(stderr,"swt->Free 2\n");
	if(PtrWtSq) free(PtrWtSq);
	if(MPWtFreq) free(MPWtFreq);
	if(MargProb) free(MargProb);
	if(ocma){ free(ofreq); free(oObs); }
	free(Weight); free(WtCnts); free(WtFreq); free(fract_seq_aln);
	if(rtn_hsw) free(rtn_hsw);
	NilSet(is); 
	if(csq) NilSeq(csq); 
	free(Observed); free(numParticipants); 
	free(order); free(info_order); free(nResTyp);
   }
}

UInt4	*swt_typ::GetIntegerWts(unsigned char ***RtnSqWt)
#if 0	//**********************************************************************
	Use integer sequence weights of 1..100 (char) so that 
	can add and subtract these counts when add or subtract the sequences...
	Only need to store query residues versus non-query residues...
	Need to use Henikoff sequence weights (either position-specific or not).
#endif	//**********************************************************************
// Both chn_pps or chn_see can get consistent weighting of sequences...
// need to merge FG & BG sequences in chn_see to get original alignment weighting.
{

	// print_error("WARNING: not tested");
	double  wt; // use swt **Weight;
	if(compute_wts) HenikoffWeights( );
	UInt4   *AveSqIWt,num_seqs=NWtSq;
	unsigned char **SqWt;
	Int4	r,sq,s,set,num_sets=2;
	UInt4 iwt;

	NEWP(SqWt,Length+3,unsigned char);
	// fprintf(stderr,"num_seqs = %d\n",num_seqs);
	for(s=1; s <= Length; s++){
	   NEW(SqWt[s],num_seqs+3,unsigned char);
#if 0	// add residue pseudocounts here???
#endif
	}
	h_type HG=Histogram("seq. wts",0,100,6.25);
	NEW(AveSqIWt,num_seqs +3, UInt4);
	for(sq=1; sq <= num_seqs; sq++){
	     double ave=0.0,total=0.0;
	     UInt4  TotalIWt=0;
	     unsigned char *seq=GetAlnResInSiteCMSA(1,sq,mcma);
	     for(s=1; s <= Length; s++){
		r = seq[s];
		if(r != 0){
		    wt = Weight[s][sq];
		    iwt = (UInt4) floor((100.0*wt) + 0.5);
		    if(iwt > 100) iwt=100;
		    SqWt[s][sq] = (unsigned char) iwt;
	            ave+= (double) iwt; total+=1.0;
		}
	     }
	     iwt = (UInt4) floor(ave/total);
	     if(iwt < 2) iwt = 2;	// minimum weight = 0.02.
	     AveSqIWt[sq]=iwt;
	     for(s=1; s <= Length; s++){
		if(seq[s] != 0) SqWt[s][sq] = iwt;
		else SqWt[s][sq] = 0;
	     }
	     IncdHist(ave/total, HG);
	     // if(ave==0.0) PutSeq(stdout,TrueSeqCMSA(sq,mcma),AB);
	     if(ave==0.0 && Verbose) PutSeq(stdout,FakeSeqCMSA(sq,mcma),AB);
		// ave may be zero if weights are very small.
	}
	if(Verbose) PutHist(stdout,60,HG); 
	NilHist(HG);
	*RtnSqWt = SqWt;
	return AveSqIWt;
}

void	swt_typ::FReadHenikoffWeights(FILE *fp)
{
	Int4	s,xLength,xNWtSq,xLenAB,rtn;
	rtn=fread(&xLength,sizeof(Int4),1,fp);
	if(rtn != 1 || xLength != Length) print_error("FReadHenikoffWeights( ) error");
	rtn=fread(&xNWtSq,sizeof(Int4),1,fp);
	if(rtn != 1 || xNWtSq!= NWtSq) print_error("FReadHenikoffWeights( ) error");
	rtn=fread(&xLenAB,sizeof(Int4),1,fp);
	if(rtn != 1 || xLenAB!= nAlpha(AB)) print_error("FReadHenikoffWeights( ) error");
        for(s=1; s <= Length; s++){
            rtn=fread(WtCnts[s],sizeof(double),nAlpha(AB)+1,fp);
	    if(rtn != nAlpha(AB)+1) print_error("FReadHenikoffWeights( ) error");
            rtn=fread(WtFreq[s],sizeof(double),nAlpha(AB)+1,fp);
	    if(rtn != nAlpha(AB)+1) print_error("FWriteHenikoffWeights( ) error");
            rtn=fread(Weight[s],sizeof(double),NWtSq+1,fp);
	    if(rtn != NWtSq+1) print_error("FWriteHenikoffWeights( ) error");
	}
        rtn=fread(fract_seq_aln,sizeof(double),Length+1,fp);
	if(rtn != Length+1) print_error("FWriteHenikoffWeights( ) error");
}

void	swt_typ::FWriteHenikoffWeights(FILE *fp)
{
      Int4 s,n,rtn;
      rtn=fwrite(&Length,sizeof(Int4),1,fp);
      rtn=fwrite(&NWtSq,sizeof(Int4),1,fp);
      rtn=fwrite(&nAlpha(AB),sizeof(Int4),1,fp);
      for(s=1; s <= Length; s++){
          rtn=fwrite(WtCnts[s],sizeof(double),nAlpha(AB)+1,fp);
	  if(rtn != nAlpha(AB)+1) print_error("FWriteHenikoffWeights( ) error");
          rtn=fwrite(WtFreq[s],sizeof(double),nAlpha(AB)+1,fp);
	  if(rtn != nAlpha(AB)+1) print_error("FWriteHenikoffWeights( ) error");
          rtn=fwrite(Weight[s],sizeof(double),NWtSq+1,fp);
	  if(rtn != NWtSq+1) print_error("FWriteHenikoffWeights( ) error");
      }
      rtn=fwrite(fract_seq_aln,sizeof(double),Length+1,fp);
      if(rtn != Length+1) print_error("FWriteHenikoffWeights( ) error");
}


