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

// Gibbs Sampling algorithms for finding multiple sites in multiple sequences 
#include "ppg_typ.h"

/*********************** propagation routines *****************************/
BooLean ppg_typ::Update(char c)
{
   Int4         k,n,t,N;

   switch(c){
     case 'p': /** update for propagation **/
        this->update_propagate_gibbs( ); return TRUE;
     case 'c': /** check data structure **/
	if(!SamePerSeqCMSA(G->cmsa)) print_error("ordering error");
        return TRUE;
     case 'i': /** initialize gibbs **/
      InitCMSA(G->cmsa); this->update_propagate_gibbs( ); return TRUE;
     case 's': /** shift occurred **/
      return FALSE;
     default: return FALSE;  /** command ignored **/
  }
}

void	ppg_typ::update_propagate_gibbs()
/*** find ends of matrix ***/
{
	Int4	n,t,blocked;

	for(blocked=1, t=1; t<=nBlksCMSA(G->cmsa); t++) {
	    for(n=1; n <= NumSeqsCMSA(G->cmsa); n++){ G->strt[t][n]=blocked; }
	    blocked += LengthCMSA(t,G->cmsa);
	}
	for(blocked=-1, t=nBlksCMSA(G->cmsa); t > 0; t--) {
	    blocked += LengthCMSA(t,G->cmsa);
	    for(n = 1; n <= NumSeqsCMSA(G->cmsa); n++)
		{ G->end[t][n] = LenSeqCMSA(n,G->cmsa) - blocked; }
	}
	for(t=nBlksCMSA(G->cmsa); t > 0; t--) {
	    for(n=1; n <= NumSeqsCMSA(G->cmsa); n++) {
		if(G->strt[t][n] > G->end[t][n]){
		  fprintf(stderr,
			"t=%d;n=%d;strt=%d;end=%d;seqlen=%d;modelen=%d\n",
			t,n,G->strt[t][n],G->end[t][n],
			LenSeqCMSA(n,G->cmsa), LengthCMSA(t,G->cmsa));
		  print_error("update_propagate_gibbs( ): seq. too short");
		}
	    }
	}
}

static Int4 pseudo_aln_score_cmsa(register Int4 Length, register unsigned char *seq1,
        register unsigned char *seq2, register char **R)
{
        // register Int4 k,score=0;
        // for(k = 0; k < Length; k++){ score += R[seq1[k]][seq2[k]]; }
        register Int4 score;
        for(score=0; Length > 0; Length--,seq1++,seq2++){
                // Skip over gap residues...
                if(*seq1 && *seq2) score += R[*seq1][*seq2];
        } return score;
}

Int4    PseudoAlnScoreCMSA_3(Int4 sq1, Int4 sq2, cma_typ cma, Int4 Blk=0)
// use cma alignment to obtain pseudo pairwise scores for two aligned sequences.
{
        // assert(!"THIS ROUTINE NEEDS TO BE DEBUGGED!\n");
        // NOTE: OFF BY ONE!!!???
        Int4	s1,s2,blk,N=NumSeqsCMSA(cma),score;
        assert(sq1 >0 && sq1 <= N && sq2 > 0 && sq2 <= N);
        a_type  A=AlphabetCMSA(cma);
        unsigned char *seq1 = SeqSeqSet(sq1,DataCMSA(cma));
        unsigned char *seq2 = SeqSeqSet(sq2,DataCMSA(cma));
	if(Blk < 1){
	   for(score=0,blk=1; blk <= nBlksCMSA(cma); blk++){
        	s1 = SitePos(blk,sq1,1,SitesCMSA(cma));
        	s2 = SitePos(blk,sq2,1,SitesCMSA(cma));
		score += pseudo_aln_score_cmsa(LengthCMSA(blk,cma),seq1+s1,seq2+s2,AlphaR(A));
	  }
	} else {
        	s1 = SitePos(Blk,sq1,1,SitesCMSA(cma));
        	s2 = SitePos(Blk,sq2,1,SitesCMSA(cma));
		score = pseudo_aln_score_cmsa(LengthCMSA(Blk,cma),seq1+s1,seq2+s2,AlphaR(A));
	} return score;
}

#if 0
Int4	ppg_typ::PropagateSimilarSeqs(Int4 percent_ident,double MaxFrctn, double &LLR)
{
	Int4	i,j,x,end,k,t,s,n,N=NumSeqsCMSA(G->cmsa);
	Int4	sq,M,ntyp = nBlksCMSA(G->cmsa),moves=0;
	double	d,dd,prb,lpr,bst_lpr=TotLikeGibbs(G);
	BooLean	*moved = new BooLean [nBlksCMSA(G->cmsa)+2];
	
        Int4  card,sq,size=(Int4) ceil((double)N*0.10);  // ten percent of the sequences.
        set_typ Set=MakeSet(N+5); ClearSet(Set);
        for(n=1;(card=CardSet(Set)) < N || n <= 12; n++){
               if(n >= 12){
                 for(sq = 1; sq <= N; sq++){
                   if(!MemberSet(sq,Set)){ ppg->Propagate(sq,&prb,moved); AddSet(sq,Set); break; }
                 } n=12;
               } else { ppg->PropagateClusters(Set,size,&prb); }

	   if((lpr=TotLikeGibbs(G)) > bst_lpr){
		   bst_lpr=lpr; SaveBestCMSA(G->cmsa);
		   fprintf(stderr,"%d. %d seqs; map = %.2f (%.2f K)\n",i,x,lpr,G->temp0);
	   } else if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); this->Update('p'); }

	   if((lpr=TotLikeGibbs(G)) > bst_lpr){
                if((nruns >= 0 || stage > 1) && map < best){
                  map=best; SaveBestCMSA(G->cmsa);
                  if(G->mod_temp && G->temp0 > 0){
                        fprintf(stderr,"map = %.2f @ %.2f K (%.2f degrees F)\n",
                                                        map,G->temp0,((9.*G->temp0)/5.)-459);
                  }
                }
           } else {
                if(map > 0 && TotLike < 0.8*map){  // then start over...
                     if(AlignedCMSA(G->cmsa)){  // if cmsa->best != null
                        InitMAPCMSA(G->cmsa); ppg->Update('p');
                     }
                 }
           } if(CardSet(Set) == N) break;
        } NilSet(Set); free(moved); LLR=bst_lpr;
}
#endif

Int4	ppg_typ::PropagateRepSets(Int4 percent_ident,double MaxFrctn, double &LLR)
{
	Int4	NumSets,sq,M,moves=0,i,j,x,end,k,t,s,N=NumSeqsCMSA(G->cmsa);
	double	d,dd,prb,lpr,bst_lpr=TotLikeGibbs(G);
	BooLean	*moved = new BooLean [nBlksCMSA(G->cmsa)+2];
	
	set_typ *set=RtnTargetSizeSetsCMSA(NumSets, percent_ident,G->cmsa,MaxFrctn);
	dh_type dH=0; // dheap(N+5,4); // ClearSet(Set);
	for(i=1; i <= NumSets; i++){
	   x=CardSet(set[i]);
	   if(x==1){ Int4 *Sq=ListSet(set[i]); this->Propagate(Sq[0],&prb,moved);  free(Sq); }
	   else {
		dH=dheap(N+5,4); 
        	for(sq=1; sq <= N; sq++){
		   if(MemberSet(sq,set[i])){
	   		d=GetTotalProbCMSA(sq, G->cmsa);
	   		insrtHeap(sq,(keytyp)d,dH); 	// align worst first...
		  }
		} this->MultiPropagate(set[i],&dd,dH); Nildheap(dH); 
	   }
	   if((lpr=TotLikeGibbs(G)) > bst_lpr){
		   bst_lpr=lpr; SaveBestCMSA(G->cmsa);
		   fprintf(stderr,"%d. %d seqs; map = %.2f (%.2f K)\n",i,x,lpr,G->temp0);
	   } else if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); this->Update('p'); }
	   NilSet(set[i]);
	} free(set); free(moved); LLR=bst_lpr;
	// if(AlignedCMSA(G->cmsa)){ InitMAPCMSA(G->cmsa); this->Update('p'); }
}

Int4	ppg_typ::PropagateClusters(set_typ Set, Int4 size, double *L)
{
	Int4	i,j,x,sq,end,k,t,s,N=NumSeqsCMSA(G->cmsa);
	Int4	Sq,M,ntyp = nBlksCMSA(G->cmsa),moves = 0;
	keytyp	key;
	
	assert(SetN(Set) > N);
	assert(size < N);
#if 1
	// 1. Find a sequence Sq not yet realigned...
	M = N-CardSet(Set);
	assert(M > 0);
	dh_type dH=dheap(N+5,4); // ClearSet(Set);
	x=random_integer(M)+1; // random 0 <= integer < max 
	for(Sq=j=0,sq =1; sq <= N; sq++){
	   if(!MemberSet(sq,Set)){ j++; if(j==x){ Sq=sq; break; } }
	}
// fprintf(stderr,"Sq=%d; M=%d; j=%d; x=%d\n",Sq,M,j,x);
	assert(Sq > 0 && Sq <= N);

h_type	*HG;	NEW(HG,ntyp+3,h_type);
#if 0	// turn on to debug and evaluate...
	Int4	SelfScore= PseudoAlnScoreCMSA_3(Sq,Sq,G->cmsa);
	double	inc=(double) SelfScore/25.0;
// HG[0]=Histogram("scores",-SelfScore,SelfScore,5);
HG[0]=Histogram("scores",-SelfScore,SelfScore,inc);
// MinValHist(HG[0]); MaxValHist(HG[0]);
// MinimumHist(H); MaximumHist(H); NumBinsHist(H);
// Int4    *RtnHist(Int4 line_leng,h_type H);

char str[30]; 
for(t=1; t <= ntyp; t++){
	sprintf(str,"Scores %d",t); HG[t]=Histogram(str,-500,500,5); 
}
HG[ntyp+1]=Histogram("prob.",-500,5000,1);
#endif

	// Find other sequences most related to Sq.
	for(sq =1; sq <= N; sq++){
	    key=(keytyp) PseudoAlnScoreCMSA_3(Sq,sq,G->cmsa);
	    if(key > 0) insrtHeap(sq,(keytyp) -key,dH);
if(HG[0]){
 IncdHist(key,HG[0]);
 for(t=1; t <= ntyp; t++)
  { key=(keytyp) PseudoAlnScoreCMSA_3(Sq,sq,G->cmsa,t); IncdHist(key,HG[t]); }
}
	}
	// Sample the best matches to the current alignment in first...
	set_typ SetX=CopySet(Set); ClearSet(SetX); AddSet(Sq,SetX);
	dh_type dH2=dheap(N+5,4); 
        for(i=1; i <= size && !emptyHeap(dH); i++){
	   sq=delminHeap(dH); 
	   key=(keytyp) GetTotalProbCMSA(sq, G->cmsa);
if(HG[0]) IncdHist(key,HG[ntyp+1]);
	   insrtHeap(sq,(keytyp)-key,dH2); AddSet(sq,Set); AddSet(sq,SetX);
	} Nildheap(dH); dH=dH2;
	// Do the sampling...
if(HG[0]){
  PutHist(stderr,60,HG[0]); 
  PutHistArray(stderr,HG[0]);
#if 1
        double mean=MeanHist(HG[0]);
        double  min=MinimumHist(HG[0]);
        double  *xvalT, *yvalT,*ProbT=0, P=0.0;
        Int4    TotalT;
        Int4 NumBinsT = RtnHistArray(&xvalT,&yvalT,&TotalT,HG[0]);
        if(TotalT > 0){
          NEW(ProbT,NumBinsT+5,double);
          double tiny=0.0000000001;
          for(i=0; i<=NumBinsT+1; i++){
              // ProbT[i] = yvalT[i]/(double) TotalT;
              ProbT[i] = ((double) yvalT[i] + tiny)/((double)TotalT + (NumBinsT+2)*tiny);
              // if(ProbT[i] > 0.0) fprintf(stdout,"%.2f: prob[%d]=%.3f\n",xvalT[i],i,ProbT[i]);
          } // fprintf(stdout,"total prob = %.3f\n",P);
        } free(ProbT);
#endif
 NilHist(HG[0]);
 for(t=1; t <= ntyp; t++){ PutHist(stderr,60,HG[t]); NilHist(HG[t]); }
 PutHist(stderr,60,HG[ntyp+1]); NilHist(HG[ntyp+1]);
}
free(HG);
	this->MultiPropagate(SetX,L,dH); Nildheap(dH); NilSet(SetX);
#else
	dh_type dH=dheap(N+5,4); ClearSet(Set);
	for(sq =1; sq <= N; sq++){
		key=(keytyp) GetTotalProbCMSA(sq, G->cmsa);
		// prob += GetProbCMSA(b,n,cma);
		insrtHeap(sq,(keytyp) -key,dH);
	}
        for(i=1; i <= size; i++){ sq=delminHeap(dH); AddSet(sq,Set); } Nildheap(dH);
	this->MultiPropagate(Set,L);
#endif
}

Int4	ppg_typ::PropagateRandom(set_typ Set, Int4 size, double *L)
{
	Int4	i,sq,end,k,t,s,N=NumSeqsCMSA(G->cmsa);
	Int4	ntyp = nBlksCMSA(G->cmsa),moves = 0;
	keytyp	key;
	assert(SetN(Set) > N);
	dh_type dH=dheap(N+5,4); ClearSet(Set);
	for(sq =1; sq <= N; sq++){ insrtHeap(sq,(keytyp) Random(),dH); }
        for(i=1; i <= size; i++){ sq=delminHeap(dH); AddSet(sq,Set); } Nildheap(dH);
	this->MultiPropagate(Set,L);
}

Int4	ppg_typ::MultiPropagate(set_typ Set, double *L, dh_type dH)
// *L = likelihood; * 
{
	Int4	i,sq,end,k,t,s,N=NumSeqsCMSA(G->cmsa);
	Int4	ntyp = nBlksCMSA(G->cmsa),moves = 0;
        double  p,q,sum;
	// BooLean	*moved; NEW(moved,ntyp+3,BooLean);
	assert(SetN(Set) > N);

	// 1. Remove sequences in Set from the alignment.
	if(dH){
	  while(!emptyHeap(dH)){
		sq=delminHeap(dH); assert(sq <= N); assert(MemberSet(sq,Set));
		for(t = 1; t <= ntyp; t++)
              	  { PosSiteCMSA(t,sq,G->pos,G->cmsa); RmSiteCMSA(t,sq,G->pos[1], G->cmsa); }
	  }
	} else {
	  for(sq=1; sq <= N; sq++){
	    if(!MemberSet(sq,Set)) continue;
	    for(t = 1; t <= ntyp; t++)
	      { PosSiteCMSA(t,sq,G->pos,G->cmsa); RmSiteCMSA(t,sq,G->pos[1], G->cmsa); }
	  }
	}

	// 2. Add sequences back in one-at-a-time...
	for(sq=1; sq <= N; sq++){
	  if(!MemberSet(sq,Set)) continue;
	  // 2.a. compute conditional probabilities.
          for(p=0.0,t=1; t<=ntyp; t++) {
             sum=this->CondProbPropagate(t,sq); 	// computes matrix...
             q=G->matrix[t][G->end[t][sq]]; 
	     p+=log(q)+log(sum);
          } *L=p;
	  // 2.b. Sample backwards through the matrix to add sites back in.
	  for(s=G->end[ntyp][sq], t = ntyp; t > 0;  t--) {
	     if(G->mod_temp == 2) s=this->BestSitePropagate(s,t,sq);
	     else s=this->ChooseSitePropagate(s,t,sq); 
	     AddSiteCMSA(t,sq,s, G->cmsa);
	     if(t > 1) s -= LengthCMSA(t-1,G->cmsa);
	  }
	} return moves;
}


Int4	ppg_typ::Propagate(Int4 sq, double *L, BooLean *moved)
// *L = likelihood; * 
{
	Int4	end,k,t,s,*oldsite=G->oldsite,ntyp = nBlksCMSA(G->cmsa),moves = 0;
        double  p,q,sum;

	// 1: remove sites from sequence n.
	for(t = 1; t <= ntyp; t++) {
	   PosSiteCMSA(t,sq,G->pos,G->cmsa);
	   s = G->pos[1]; RmSiteCMSA(t,sq,s, G->cmsa);
	   oldsite[t]=s;
	}

	// 2: compute conditional probabilities.
        for(p=0.0,t=1; t<=ntyp; t++) {
             sum=this->CondProbPropagate(t,sq); 
             q=G->matrix[t][G->end[t][sq]]; p+=log(q)+log(sum);
        } *L = p;

	// 3: sample backwards through the matrix adding sampled 
	//	sites back into sequence.
	for(s=G->end[ntyp][sq], t = ntyp; t > 0;  t--) {
	     if(G->mod_temp == 2) s=this->BestSitePropagate(s,t,sq);
	     else s=this->ChooseSitePropagate(s,t,sq); 
	     if(s != oldsite[t]) { moved[t]=TRUE; moves++; } else moved[t]=FALSE;
	     AddSiteCMSA(t,sq,s, G->cmsa);
	     if(t > 1) s -= LengthCMSA(t-1,G->cmsa);
	}
	return moves;
}

Int4     ppg_typ::ChooseSitePropagate(register Int4 end, Int4 t, Int4 n)
// sample a t site in sequence n of S. return site location
// WARNING: Assumes that BLOCKED sites have zero probability. 
{
        register double	rand_no, cum_prob;
        register Int4	site;
        register double	*P=G->matrix[t];

        if(end <= G->end[t][n]){
          do {
		rand_no = (double) Random()/(double) RANDOM_MAX;
          } while(rand_no == 0.0);
          rand_no *= P[end];	/** cumulative probability **/
          for(site=G->strt[t][n],cum_prob = 0.0; site <= end; site++){
           if((cum_prob += (double) (P[site]-P[site-1])) >= rand_no)
			return site;
          }
	} else print_error("ChooseSitePropagate( ): end error");
        fprintf(stderr,
		"site = %d; n = %d; start = %d; end = %d (=%d?)\n", 
		site,n,G->strt[t][n],G->end[t][n],end);
        fprintf(stderr,
		"block %d; total prob = %g; rand_no = %g; cum_prob = %g\n",
                t,P[end],rand_no,cum_prob);
	for(site=G->strt[t][n]; site <=  G->end[t][n]; site++){
		fprintf(stderr," s = %d; P[s] = %g\n",site,P[site]);
	}
	// PutMSA(stderr, G->cmsa);
	PutSeqAlnCMSA(stderr, G->cmsa);
	assert(cum_prob < rand_no);
        print_error("ChooseSiteGibbs( ) - this should not happen!?");
}

Int4	ppg_typ::BestSitePropagate(register Int4 end, Int4 t, Int4 n)
// Pick the best type t site in sequence n of S. return site location
//  WARNING: Assumes that BLOCKED sites have zero probability. 
{
        register double	max_prob;
        register Int4	max_site,site;
        register double	*P = G->matrix[t];

        if(end <= G->end[t][n]){
           for(max_site=site=G->strt[t][n],max_prob=0.; site<=end; site++){
		if(max_prob < (double) (P[site]-P[site-1])){
			max_prob=P[site]-P[site-1]; max_site=site;
		}
           }
           return max_site;
	} print_error("BestSiteGibbs( ) - this should not happen!?");
}

double	ppg_typ::LogLikePropagate()
{
	Int4	end,k,t,s,ntyp,n;
	double	p,q,sum;

	ntyp = nBlksCMSA(G->cmsa);
	for(p=0,n = 1; n <= NumSeqsCMSA(G->cmsa); n++){
	   for(t=1; t<=ntyp; t++) {
	     sum=this->CondProbPropagate(t,n); end=G->end[t][n];
	     q=G->matrix[t][end]; p+=log(q)+log(sum);
	   }
	}
	return (1.4427*p);
}

double	ppg_typ::CondProbPropagate(Int4 t, Int4 n)
/*********************************************************************

       A          B             C

    [P(A1)][       0     ][      0      ]
    [P(A2)][       0     ][      0      ]
       :       :      :      :
       :       :      :      :
    [P(Ai)][P(Bi,A<i-end)][      0      ]
       :       :      :      :
       :       :      :      :
    [P(Ai)][P(Bi,A<i-end)][P(Ci,Bj,Ak)  ]

	P(A,B,C) = P(A) * P(B|A) * P(C|B).

   Note: P(M0,i) = 1.0 & THERE ARE ROUNDING ERRORS but this is dealt with.

   end[t][n] = end for sequence n and model t.

 *********************************************************************/
{
        register unsigned char    *seq;
        register Int4    i,end,len;
        register double  sum,*P = G->matrix[t];
        register double  *Pm1 = G->matrix[t-1];
	double		normL;
	fm_type		M=ModelCMSA(t,G->cmsa);

	/*** before convergence use seg'ed sequences ***/
        if(G->stage < 1) seq=XSeqPtrCMSA(n,G->cmsa);
	/*** but after convergence use 'raw' sequences ***/
	else seq = SeqPtrCMSA(n,G->cmsa);
	/*** set prob = 0 for i = 1..start-1 ****/
        for(i=0; i < G->strt[t][n]; i++) P[i] = 0.0;

	/*** compute P(Mi) ****/
        end = G->end[t][n];
  if(G->mod_temp && G->temp0 > 0){  // WARNING: NOT SURE THAT THIS IS CORRECT
  // if(G->mod_temp){
    normL=NormLikelihoodFModel(M);
    do {
        for(sum=0.0, i=G->strt[t][n]; i<= end; i++){
           P[i] = pow((LikelihoodFModel(seq,i,M)/normL),G->temp);
	   sum += P[i];		
        }
	if(sum >= DBL_MAX){ 
	    fprintf(stderr,
		"!!WARNING: overflow in CondProbPropagate( ) (norm=%g)!!\n", normL);
	    fprintf(stderr,"Sampling Temperature = %.2f K\n",G->temp0);
	    fprintf(stderr,"Adjusting likelihoods to compensate...\n");
		normL *= 100000.0;	/** increase normL **/
	} else if(sum == 0.0){
	    fprintf(stderr,
		"!!WARNING: underflow in CondProbPropagate( ) (norm=%g)!!\n", normL);
	    fprintf(stderr,"Sampling Temperature = %.2f K\n",G->temp0);
	    fprintf(stderr,"Adjusting likelihoods to compensate...\n");
		normL /= 100000.0;	/** decrease normL **/
	} else break;
    } while(TRUE);
  } else {
        for(sum=0.0, i=G->strt[t][n]; i<= end; i++){
           P[i] = (double)LikelihoodFModel(seq, i, M);
	   sum += P[i];		
        }
  }
	if(sum >= DBL_MAX) print_error("sum > DBL_MAX");
	// normalize probabilities.
        for(i=G->strt[t][n]; i<= end; i++) { P[i] /= sum; }

	// compute Markovian conditional probability.
	if(t > 1) len = LengthCMSA(t-1,G->cmsa); else len=0;

        for(i=G->strt[t][n]; i<= end; i++){ P[i] *= Pm1[i-len]; }
#if 0
/**** add gap penalty right here based on distance --v ****/
???????
        for(i=G->strt[t][n]; i<= end; i++){ P[i] *= Pm1[i-len]*gapfunct[i-len]; }
???????
#endif

/*** treat gaps as another residue (insertions or deletions) so
     that in addition to summing probabilities over all residue positions
     also need to sum over gap positions as well.  ************/

	// convert to cumulative probabilities.
        for(P[0] = 0.0,i=G->strt[t][n]; i<= end; i++){ 
		P[i] += P[i-1]; 
#if 0
	if(P[i] >= 1.0) fprintf(stderr,"?P[%d][%d][%d] = %g (%g + %g)\n",
			t,n,i,P[i],P[i]-P[i-1],P[i-1]);
#endif
	}
	return sum;
}

