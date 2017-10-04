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

#include "psm_typ.h"

psm_typ::psm_typ(Int4 maxrpts,cma_typ cma)
// -G200..1300,20..100/0..50:500,35..250
{ init(maxrpts,"-P200..1300,20..100/0..50:500,35..250",'D', 693, cma); }
// 693 pernats == 1000 perbits.

psm_typ::psm_typ(Int4 maxrpts,char *pssm_arg, cma_typ cma)
// input as "-P100..1100,20..120:500,40..400"
// or as -P150..1400,15..70/5..70:500,40..400  (Best settings??)
{ init(maxrpts,pssm_arg,'D',693,cma); }		// 693 pernats == 1000 perbits.

psm_typ::psm_typ(Int4 maxrpts,char *pssm_arg,cma_typ cma,Int4 pernats)
{ init(maxrpts,pssm_arg,'D', pernats, cma); }

void	psm_typ::GetProb(Int4 left_flank,Int4 right_flank, cma_typ cma)
{
    Int4        j,t,k;
    double	predicted_accuracy,**tprob;

    if(nblk == 1){
	tprob=ComputeDSC(NULL,&predicted_accuracy,nblk,left_flank,right_flank,cma);
	prob_h=tprob[1]; prob_s=tprob[2]; prob_c=tprob[3]; free(tprob);
// for(Int4 i=0;i<=length;i++) fprintf(stderr,"prob(%d) = %.2f %.2f %.2f\n",i,prob[1][i],prob[2][i],prob[3][i]);
    } else {
      NEW(prob_h,length+4,double);
      NEW(prob_s,length+4,double);
      NEW(prob_c,length+4,double);
      for(k=1,t=1; t <= nblk; t++) {
	tprob=ComputeDSC(NULL,&predicted_accuracy,t,left_flank,right_flank,cma);
        for(j=1; j <= lenblk[t]; j++) {
          prob_h[k] = tprob[1][j]; 	// Helix
          prob_s[k] = tprob[2][j];	// Stand
          prob_c[k] = tprob[3][j]; 	// coil probability
	  k++;
        } free(tprob[1]); free(tprob[2]); free(tprob[3]); free(tprob);
      }
    }
#if 1	// New bayesian update using observed structures...
	// 1. Scan cma for structural information and use to update probabilities.
	Int4	pos[4],s,sq,ts,b;
	char	ss;
	BooLean	flag=0;
	double	numObs=1.0;	// one pseudocount for now...
	gss_typ *gss=gssCMSA(cma);
	for(sq=1; sq <= NumSeqsCMSA(cma); sq++){
	    e_type E = gss->TrueSeq(sq);
	    if(IsStructSeq(E)){
		flag=1;
		numObs+=1.0;
		for(k=1,b=1; b<=nblk; b++){
		   if(PosSiteCMSA(b,sq,pos,cma)){
		     for(s=1; s<=lenblk[b]; s++){
			if(!IsDeletedCMSA(sq,s,cma)){
		           ts=gss->TruePos(sq,pos[1]+s-1);
			   ss=StructSeq(ts,E);
			   switch(ss){
			     case 'H': prob_h[k] +=0.333; break;
			     case 'h': prob_h[k] +=0.222; prob_c[k] +=0.111; break;
			     case 'S': prob_s[k] +=0.333; break;
			     // case 's': prob_s[k] +=0.167; prob_c[k] +=0.083; break;
			     case 's': prob_s[k] +=0.222; prob_c[k] +=0.111; break;
			     case 'C': prob_c[k] +=0.333; break;
			     case 'c': prob_c[k] +=0.333; break;
			     case 'u': break;	// do nothing 
			     default: print_error("psm_typ::GetProb() input error");
				break;
			   }
			} k++;
		     }
		   }
		}
	    }
	}
	if(flag){	// then normalize probabilities...
	   double	denom;
	   for(k=1,b=1; b<=nblk; b++){
		for(s=1; s<=lenblk[b]; s++){
		   denom= prob_h[k] + prob_s[k] + prob_c[k];
		   prob_h[k] = prob_h[k]/denom;
		   prob_s[k] = prob_s[k]/denom;
		   prob_c[k] = prob_c[k]/denom;
		   k++;
		}
	   }
	}
#endif
}

void	psm_typ::init(Int4 maxrpts, char *psm_arg, char mode, Int4 pernats, cma_typ cma)
{
	Int4 	left_flank=0,right_flank=0;

	target_pvalue=0.05; 
	length=TotalLenCMSA(cma);
	nblk=nBlksCMSA(cma); PerNats = pernats; MaxRpts = maxrpts;
	NEW(lenblk,nblk+3,unsigned short);
	for(Int4 b=1; b<=nblk; b++) lenblk[b]=(unsigned short)LengthCMSA(b,cma);

        // psg_typ psg(cma,pernats);   
	psg = new psg_typ(cma,'h','h',0,pernats,'w');
	M = psg->BlcksMatrices(maxrpts);
	a_type	A = AlphabetCMSA(cma);
	if(mode = 'D'){ GetProb(left_flank,right_flank,cma); }
	else {		// Want else prob=ComputeGOR(cma); // that calls gor_typ g;
	   print_error("GOR not currently implemented");
		// gor_typ g; double ***inf = g.FanoInfo(cma); prob = g.AlnProbab(cma, inf); 
	} 
	idp = new idp_typ(psm_arg,'S',prob_c,length,maxrpts,nAlpha(A),AlphaSurface(A),
			nblk,lenblk,pernats);
	// idp->Put(stderr);
	consensus = ConsensusSeqCMSA(cma); 
	Put(stderr);
}

void	psm_typ::Put(FILE *fp)
{
	char	c;
	Int4	i,j,b;
	fprintf(fp,"Insertion & Deletion penalties:\n");
        fprintf(fp," POS:   m2i  i2i  m2d  d2d  i2m  m2m  d2m  (coil strand helix) SS\n");
	for(b=1,j=0,i=1;i<=length;i++){
            fprintf(fp,"%c %-3d: %4d %4d %4d %4d %4d %4d %4d  (%.2f  %.2f  %.2f)",consensus[i],i,
		-idp->MatToIns(i),-idp->InsToIns(i),-idp->MatToDel(i),
		-idp->DelToDel(i),-idp->InsToMat(i),-idp->MatToMat(i),
		-idp->DelToMat(i),prob_c[i],prob_s[i],prob_h[i]);
            if(prob_c[i] > prob_s[i] && prob_c[i] > prob_h[i]) c='C';
            else if(prob_s[i] > prob_c[i] && prob_s[i] > prob_h[i]) c='S';
            else if(prob_h[i] > prob_c[i] && prob_h[i] > prob_s[i]) c='H';
            else c=' ';
            fprintf(fp,"   %c\n",c);
	    j++;
	    if(j==lenblk[b]){ assert(b <= nblk); b++; j=0; fprintf(fp,"\n"); }
        }
	// idp->Put(fp);
	// fprintf(fp,"Scoring matrix:\n");
	// for(Int4 b=1; b <= nblk; b++) PutSMatrix(fp,M[b]);
}

#if 0
void	psm_typ::PutAlign(FILE *fp, e_type E1, char *operation, Int4 operation_length, 
		Int4 start, Int4 rpts)
{
        Int4    i, x;
        Int4    tmp = 0, e1 = start, e2 = start-1; 
        Int4    j = 1, s = 1, t = 1, k = start;
        e_type  E2 = ConsensusSMatrix(M, nblk, rpts);
        Int4    rpt_len = LenSeq(E2)/rpts;
        Int4    *lengths;
	a_type 	A = SMatrixA(M[1]);

        NEW(lengths,rpts + 1,Int4);
        unsigned char *ptr1 = SeqPtr(E1);
        unsigned char *ptr2 = SeqPtr(E2);
        char *seq1, *seq2, *beetw;

        NEW(seq1,operation_length+2,char);
        NEW(seq2,operation_length+2,char);
        NEW(beetw,operation_length+2,char);
        for(i=1; operation[i]!='E'; i++){
	  switch(operation[i]){
	     case 'M': case 'm':
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = AlphaChar(ptr2[j],A);
                        if(ptr1[k] == ptr2[j]) beetw[i] = ':';
                        else if(ValSMatrix(s,ptr1[k],M[t]) > 0) beetw[i] = '.';
                        else beetw[i] = ' ';
                        j++;k++;
                        if(s>=LenSMatrix(M[t])) {s=1;t++;}
                        else s++;
		break;
	     case 'D': case 'd':
                        seq1[i] = '-'; seq2[i] = AlphaChar(ptr2[j],A); beetw[i] = ' ';
                        j++;
                        if(s>=LenSMatrix(M[t])) {s=1;t++;} else s++;
		break;
	     case 'I': case 'i':
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = '-'; beetw[i] = ' ';
                        k++;
		break;
	     default: 
		 std::cerr << operation; std::cerr << std::endl; std::cerr << operation[i]; std::cerr << std::endl;
			 print_error("PutAlign(): error in operation array");
		break;
	   }
        } fprintf(fp,"\n\n");
        Int4 br = 0; j=1; 
        for(j=1,i=1; i < operation_length-1; ){
                while(br < rpt_len){
                        if (seq2[i++] != '-') br++; else lengths[j] += 1;
                } br = 0; j++;
        } Int4 total=0;
        for(i=1; i<=rpts; i++){
                x=1;
                for(j=1; j<=rpt_len + lengths[i]; j++){
                        if (seq2[total+j] == '-') x++; else break;
                }
                for(j=x; j<=rpt_len + lengths[i]; j++) fprintf(fp,"%c",seq2[total+j]);
                fprintf(fp,"\n");
                for(j=x; j<=rpt_len + lengths[i]; j++)
                        fprintf(fp,"%c",beetw[total+j]);
                fprintf(fp,"\n");
                for(tmp=0,j=x; j<=rpt_len + lengths[i]; j++){
                        fprintf(fp,"%c",seq1[total+j]);
                        if (seq1[total+j]!='-') tmp++;
                }
                e2+=(x+tmp-1);e1=e2-tmp+1;
                fprintf(fp," %d-%d",e1,e2);
                fprintf(fp,"\n\n\n");
                total += rpt_len + lengths[i];
        } free(seq1); free(seq2); free(beetw); free(lengths); NilSeq(E2);
}
#endif

double	psm_typ::RelMap(cma_typ cma) { return RelMap(0,cma); }

double	psm_typ::RelMap(Int4 rm_seq, cma_typ cma)
{
    double      indel_penalty=0.0;
    double      interblk_penalty=0.0;
    double	map = UnGappedRelMapCMSA(cma);
    gss_typ	*gss=gssCMSA(cma);
    Int4	start[3];

    // assert(nBlksCMSA(cma) == 1);
    for(Int4 s = 1; s <= NumSeqsCMSA(cma); s++){
      if(s != rm_seq){
      	for(Int4 b=1; b <= nblk; b++){
	    PosSiteCMSA(b,s,start,cma);
	    indel_penalty -= (double)gss->GapPenalty(s,start[1],b,idp)/(double)PerNats; 
#if 1	// between block penalties...
	    if(b < nblk){
		Int4 end = idp->EndBlk(b);
		Int4 gap = ResBetweenCMSA(b,s,cma);
	    	interblk_penalty -= gap*((double)idp->InsToIns(end)/(double)PerNats);
// fprintf(stderr,"InsToIns(%d) = %d; gap = %d\n",end,idp->InsToIns(end),gap);
	    }
#endif
	}
      }
    }
#if 0
    fprintf(stderr,"map = %g; penalty=%g; interblk=%g\n",
			map,indel_penalty,interblk_penalty);
#endif
    return (map-indel_penalty-interblk_penalty);
}

char	*psm_typ::FindBestRpts(e_type E, Int4 start_rpts, Int4 *Rpts, 
	Int4 *Score, Int4 *Start, Int4 *OperLen)
// Find the 
{
	char	*low_operation=0,*high_operation,*best_operation;
	Int4	low_score,high_score;
	Int4	low_start,high_start;
	Int4	low_oper_len,high_oper_len;
	Int4	high_rpts,low_rpts;
	Int4	score,threshold,maxscore,oldscore;
	double	pval,oldpval;

        double *Pvalues = ThresholdScore(E,target_pvalue,&threshold,&maxscore);
	// PutSeqInfo(stderr,E);
        fprintf(stderr,"ThresholdScore = %d; maxscore = %d; target_pvalue = %g\n",
		threshold,maxscore,target_pvalue);

	assert(start_rpts > 0); 
	high_rpts=start_rpts; 
        high_operation=Align(E,high_rpts,&high_score,&high_start,&high_oper_len);
#if 0	// Fix for multiple blocks
	score = high_score; 
	if(score <= maxscore){
                if(score >= threshold) pval = Pvalues[score]; else pval=1.0;
        } else pval = 0.0;
#endif
	// 1. Find lower limit...
	for(low_rpts=high_rpts-1; low_rpts >= 0; low_rpts--){
	  if(low_rpts > 0){
	     low_operation=Align(E,low_rpts,&low_score,&low_start,&low_oper_len);
          } else { low_score=0; low_operation=0; }
          score = high_score-low_score;

	  if(score <= maxscore){
                if(score >= threshold) pval = Pvalues[score]; else pval=1.0;
          } else pval = 0.0;
	  if(pval <= target_pvalue){	// then high repeats are okay.
		break;
	  } else {
                fprintf(stderr,"Rpts(low) = %d; Score = %d (%.3f)\n",high_rpts,score,pval);
		free(high_operation);
		high_rpts=low_rpts; high_operation=low_operation; 
		high_score=low_score; high_start=low_start; 
		high_oper_len=low_oper_len;
		if(high_rpts==0) score = 0;
	  }
	}
	if(high_rpts < start_rpts || high_rpts == MaxRpts){
	  if(low_operation) free(low_operation);
	  *Rpts = high_rpts; *Score= high_score; *OperLen=high_oper_len; 
          fprintf(stderr,"Final*: Rpts = %d; Score = %d (%.3f); Tot = %d ",
				high_rpts,score,pval,*Score);
	  // fprintf(stderr,"\n");
	  PutSeqInfo(stderr,E);
	  *Start=high_start; best_operation=high_operation;
	} else {			// 2. Find upper limit...
	  while(high_rpts < MaxRpts){
             if(high_rpts > 1) 
		fprintf(stderr,"Rpts(high) = %d; Score = %d (%.3f)\n",high_rpts,score,pval);
	     if(low_operation) free(low_operation);
	     low_rpts=high_rpts; low_operation=high_operation; 
	     low_score=high_score; low_start=high_start; low_oper_len=high_oper_len;
	     oldpval = pval; oldscore=score;
	     high_rpts++;
             high_operation=Align(E,high_rpts,&high_score,&high_start,&high_oper_len);
             score = high_score-low_score;
	     if(score <= maxscore){
                if(score >= threshold) pval = Pvalues[score]; else pval=1.0;
             } else pval = 0.0;
	     if(pval > target_pvalue){	// then high repeats are bad.
		free(high_operation);
		score = oldscore; pval=oldpval; break;
	     } else if(high_rpts == MaxRpts){ // then high repeats are good on exit.
	        if(low_operation) free(low_operation);
	        low_rpts=high_rpts; low_operation=high_operation; 
	        low_score=high_score; low_start=high_start; low_oper_len=high_oper_len;
		break;
	     }
	  }  
	  *Rpts = low_rpts; *Score= low_score; *OperLen=low_oper_len; *Start=low_start;
          fprintf(stderr,"Final: Rpts = %d; Score = %d (%.3f); Tot = %d ",
					low_rpts,score,pval,*Score);
	  // fprintf(stderr,"\n");
	  PutSeqInfo(stderr,E);
	  best_operation=low_operation;
	} free(Pvalues);
	return best_operation;
}

double	*psm_typ::ThresholdScore(e_type E, double TargetPvalue, Int4 *TargetScore,
	Int4 *MaxScore)
// compute P-values via Monte Carlo simulation
{
        char	*operation;
        Int4	Score,Start,Oper_len,tmp_score,hpsz,item,N,factor=25;
	Int4	max,s;
	double	pval,*Pvalues;

	assert(TargetPvalue > 0.0 && TargetPvalue < 0.99);

        e_type	sE = CopySeq(E);
	hpsz = (factor*(Int4) ceil(1.0/TargetPvalue)) + factor;
	dh_type dH = dheap(hpsz+3,4);
        h_type H = Histogram("scores",-1000,5000,100);
	max=-99999;
        for(item=1;item<=hpsz;item++){
              sE = ShuffleSeq(sE);
	      operation=Align(sE,1,&Score,&Start,&Oper_len);
              free(operation);
              IncdHist(Score,H);
	      insrtHeap(item,(keytyp)-Score,dH);
	      if(max < Score) max = Score;
        } NilSeq(sE);
	if(max < 0) max = 0; 
	*MaxScore = max;
	NEW(Pvalues,max+5,double);
	for(pval=0.0,Score=max,N=0; !emptyHeap(H); ){
	      tmp_score=-minkeyHeap(dH);
	      if(tmp_score < 0) break; // fixes bug???
	      item=delminHeap(dH); N++;
	      pval = (double) N/ (double) hpsz;
	      for(s=max; s>=tmp_score; s--) Pvalues[s]=pval; max=tmp_score;
	      if(pval >= TargetPvalue && Score > tmp_score) break;
	      else Score=tmp_score;
	} PutHist(stderr,60,H);
	NilHist(H); Nildheap(dH);
	*TargetScore=Score;
	return Pvalues;
}

//============== Aleksandar's Complex alignment routines ================

static void	complex_affine_aln_psm(Int4 *matj, Int4 *matjm1, 
	   Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
           Int4  delo, Int4  dele, Int4  inso, Int4 ie, 
	   Int4 seq_len, unsigned char *seq, Int4 *scorej)
{
   register Int4	i,im1,s,Mjm1_im1,Djm1_im1,Ijm1_im1;

   for(im1=0,i=1;i<=seq_len;im1++,i++){
	Ijm1_im1 = insjm1[im1]; Djm1_im1 = deljm1[im1]; Mjm1_im1 = matjm1[im1];

        if(Mjm1_im1 >= Djm1_im1 && Mjm1_im1 >= Ijm1_im1)
		matj[i] = scorej[seq[i]]+Mjm1_im1;			// m2m
        else if(Djm1_im1 >= Ijm1_im1) matj[i]=scorej[seq[i]]+Djm1_im1;	// d2m
        else matj[i] = scorej[seq[i]]+Ijm1_im1;				// i2m
        if((s=delo+matjm1[i]) > deljm1[i]) delj[i]=s+dele;		// m2d
	else delj[i]=deljm1[i]+dele;					// d2d
        if((s=inso+matj[im1]) > insj[im1]) insj[i]=s+ie;		// m2i
	else insj[i]=insj[im1]+ie;					// i2i
   }
}

static void	complex_affine_aln_psm(Int4 *matj, Int4 *matjm1, 
	   Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
           Int4  delo, Int4  dele, Int4  inso, Int4 ie, 
	   Int4 seq_len, unsigned char *seq, Int4 *scorej, Int4 *grease_pen)
// about 80% of compute time is spent here...
{
   register Int4	i,im1,s,Mjm1_im1,Djm1_im1,Ijm1_im1;

   for(im1=0,i=1;i<=seq_len;im1++,i++){
	Ijm1_im1 = insjm1[im1]; Mjm1_im1 = matjm1[im1]; Djm1_im1 = deljm1[im1]; 

        if(Mjm1_im1>=Djm1_im1 && Mjm1_im1>=Ijm1_im1)
		matj[i] = scorej[seq[i]]+Mjm1_im1;			// m2m
        else if(Djm1_im1 >= Ijm1_im1) matj[i]=scorej[seq[i]]+Djm1_im1;  // d2m
        else matj[i] = scorej[seq[i]]+Ijm1_im1;				// i2m
        if((s=delo+matjm1[i]) > deljm1[i]) delj[i]=s+dele;		// m2d
	else delj[i]=deljm1[i]+dele;					// d2d
        if((s=inso+matj[im1]) > insj[im1]) insj[i]=s+grease_pen[i]+ie;  // m2i
	else insj[i]=insj[im1]+grease_pen[i]+ie;			// i2i
   }
}

static void	complex_gap_align_psm(Int4 *matj, register Int4 *matjm1, 
	   Int4 *insj, Int4 *delj, register Int4 *deljm1,
           Int4  delpen, Int4  inso, Int4 ie, Int4 seq_len, 
	   unsigned char *seq, Int4 *scorej, Int4 *gpen_s,Int4 *grease_pen)
{
   register Int4 pen,max,i0,g,s;

   for(Int4 im1=0,i=1;i<=seq_len;im1++,i++){

        for(max=SHRT_MIN,i0=im1,g=0;g < i;g++,i0--){
            if((pen=gpen_s[g]) <= SHRT_MIN) break;
            if(max < (s=pen+matjm1[i0])) max=s;
            if(max < (s=pen+deljm1[i0])) max=s;
        } matj[i] = scorej[seq[i]]+max;

        for(max=SHRT_MIN,i0=i,g=0;g <= i;g++,i0--){
            if((pen=gpen_s[g]) <= SHRT_MIN) break;
            if(max < (s=pen+matjm1[i0])) max=s;
            if(max < (s=pen+deljm1[i0])) max=s;
        } delj[i] = delpen+max;

	if(grease_pen){
           if((s=inso+matj[im1]) > insj[im1]) insj[i]=s+ie + grease_pen[i];
	   else insj[i]=insj[im1]+ie + grease_pen[i];
	} else {
           if((s=inso+matj[im1]) > insj[im1]) insj[i]=s+ie; else insj[i]=insj[im1]+ie;
	}
   }
}

static char *get_traceback_operation_psm(Int4 seq_len, Int4 prof_len,
        Int4 nblks, Int4 *block_lengths, Int4 *start_prof, Int4 **MAT, 
        Int4 **DEL, Int4 **INS, Int4 *oper_len, Int4 *grease_pen, 
        Int4 *alignscore, Int4 *start, idp_typ *idp)
{
        Int4    l,s,i,i0,j,k,bl_nmbr,best_i,inse,end;
        Int4    *gapfunct,gapend,i_gap,gapstart;
        Int4    mxm, gap,g;
        char    state,*operation, *back_operation;
        Int4 	*io=idp->MatToIns();
        Int4 	*ie=idp->InsToIns();
        Int4 	*od=idp->MatToDel();
        Int4 	*de=idp->DelToDel();
        Int4 	**gpen=idp->GapFunct();

        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
        j=prof_len; i=seq_len;

        for(state='E',k=0;state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                        back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
                        back_operation[k++]='E';
                        bl_nmbr=nblks; s=start_prof[bl_nmbr];
                        mxm=MAT[j][i];state='M';i0=i;
                        for(i=1;i<=seq_len;i++){
                           if(MAT[j][i]>mxm){ state='M'; mxm=MAT[j][i];i0=i;}
                           if(DEL[j][i]>mxm){ state='D'; mxm=DEL[j][i];i0=i;}
                        } i=i0;*alignscore=mxm;
                        break;
                case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        if(j==s && gpen[bl_nmbr-1]!=0){
                           bl_nmbr--; s=start_prof[bl_nmbr];
                           j--; back_operation[k++]='M'; 
                           state='G'; gapstart=i-1;  
                        } else { // else use affine gaps.
			   if(j==s) {back_operation[k++]='M'; bl_nmbr--; s=start_prof[bl_nmbr];}
			   else back_operation[k++]='m';
                           j--; i--; 
                           if(MAT[j][i]>INS[j][i] && MAT[j][i]>DEL[j][i]) state='M';
                           else if(INS[j][i]>DEL[j][i]) state='I';
                           else state='D'; 
                        } break;
                case 'D': // previously sampled deletion
                        if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                        if(j==s && gpen[bl_nmbr-1]!=0){
                            j--; back_operation[k++]='D'; 
                            bl_nmbr--; s=start_prof[bl_nmbr];
                            state='G'; gapstart=i; 
                        } else { 
			     if(j==s) {back_operation[k++]='D'; bl_nmbr--; s=start_prof[bl_nmbr];}
			     else back_operation[k++]='d';
			     j--;
                             if(DEL[j][i]+de[j]>=MAT[j][i]+od[j]+de[j]) state='D';
                             else state='M';
                        } break;
                case 'I': // previously sampled insertion.
                          i--; end = start_prof[bl_nmbr]+block_lengths[bl_nmbr]-1;
                          if(grease_pen){
                                if(j==end) {inse = ie[j];}
                                else { inse = ie[j] + grease_pen[i]; }
                          } else inse = ie[j];
			  if(j==end) {back_operation[k++]='i'; *alignscore -= ie[j];}
			  else back_operation[k++]='I'; 
                          if (INS[j][i]+inse > MAT[j][i]+io[j]+inse) state='I';
                          else state='M'; 
                          break;                   
                case 'G':  // use gap function (between blocks).
                       	  gapfunct=gpen[bl_nmbr];
                      	  gapend = i_gap = gapstart;
                      	  mxm=gapfunct[0]+MAT[j][i_gap];
                      	  for(g=0; g <= gapend; g++,i_gap--){
                       		if(gapfunct[g] <= SHRT_MIN) break; // fixes bug???
                       		gap=gapfunct[g]+MAT[j][i_gap];
                       		if(gapfunct[g]+MAT[j][i_gap]>=mxm) {best_i=i_gap; state='M';mxm=gap;}
                       		gap=gapfunct[g]+DEL[j][i_gap];
                       		if(gapfunct[g]+DEL[j][i_gap]>=mxm) {best_i=i_gap; state='D';mxm=gap;}
                          } g = (gapstart-best_i);
                          for(l=0;l < g; l++) { back_operation[k++]='i';} i =best_i; 
                  	  break;
                default:  print_error("this should not happen");
           }
        } *oper_len = k;
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
        return operation;
}

char	*psm_typ::AlignHMM(e_type E, Int4 Rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len)
{
        Int4    i,j,jm1,n,r,s,prof_len;
        Int4    alph_len,total;
        Int4    *block_lengths, *start_prof;
        Int4    **MAT,**DEL,**INS,**score;
        Int4    nblks = Rpts*idp->nBlks();
	Int4	seq_len=LenSeq(E);
	unsigned char *seq=XSeqPtr(E);
         
        assert(Rpts <= idp->MaxRpts());
        MEW(block_lengths, nblks+1, Int4);
        MEW(start_prof, nblks+2, Int4);

        for(n=0,i=1;i<=nblks;i++){
                start_prof[i]=n+1;
                block_lengths[i] = LenSMatrix(M[i]);
                n+=block_lengths[i];
        } start_prof[nblks+1]=0;
        prof_len=n; 

        Int4    *grease_pen=0;
        grease_pen=idp->SeqPenalties(seq,seq_len); // returns 0 if ixal==0 && ixah==0;
        alph_len = nAlpha(SMatrixA(M[1]));
        NEWP(MAT,prof_len+1,Int4); NEWP(DEL,prof_len+1,Int4);
        NEWP(INS,prof_len+1,Int4);
        for(i=0;i<=prof_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
                MEW(INS[i],seq_len+1,Int4);
        }
        NEWP(score,prof_len+1,Int4);
        for(i=0;i<=prof_len;i++) NEW(score[i],alph_len+1,Int4);
        for(total=0,i=1;i<=nblks;i++){
           for(s=1;s<=block_lengths[i];s++){
                for(r=0;r<=alph_len;r++){
                   score[total+s][r]=ValSMatrix(s,r,M[i]);
                }
           } total+=block_lengths[i];
        }

        Int4 *io=idp->MatToIns();
        Int4 *ie=idp->InsToIns();
        Int4 *od=idp->MatToDel();
        Int4 *de=idp->DelToDel();
        Int4 **gapen=idp->GapFunct();

        DEL[0][0]=INS[0][0]=MAT[0][0]=0;
        DEL[1][0]=INS[1][0]=MAT[1][0]=od[1]+de[1];
        for(s=2,jm1=1,j=2;j<=prof_len;jm1++,j++) {
                if(j!=start_prof[s]){ MAT[j][0]=INS[j][0]=DEL[j][0]=DEL[jm1][0]+de[j];
                } else {
                   if(gapen[s-1]){
                        MAT[j][0]=DEL[jm1][0]+od[j]+de[j]+gapen[s-1][0]; 
                        DEL[j][0]=DEL[jm1][0]+od[j]+de[j]+gapen[s-1][0];
                        INS[j][0]=DEL[jm1][0]+od[j]+de[j]+gapen[s-1][0];
                   } else {
                        MAT[j][0]=DEL[jm1][0]+od[j]+de[j]; 
                        DEL[j][0]=DEL[jm1][0]+od[j]+de[j];
                        INS[j][0]=DEL[jm1][0]+od[j]+de[j];
                   } s++; 
                }
        } for(i=1;i<=seq_len;i++) {DEL[0][i] = od[1];  MAT[0][i] = 0; INS[0][i] = 0;}

        register Int4 J,Jm1;
        for(s=1,J=1,Jm1=0;J<=prof_len;J++,Jm1++){
           Int4 End = start_prof[s]+block_lengths[s]-1;
           if(J!=start_prof[s+1] || gapen[s]==0){
              if((grease_pen && J!=End)){
                complex_affine_aln_psm(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J],grease_pen);
              } else { // don't use grease penalities between blocks...
                complex_affine_aln_psm(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J]);
              }
           } else { // impose gap penalties between blocks when at first position...
              complex_gap_align_psm(MAT[J],MAT[Jm1],INS[J],DEL[J],DEL[Jm1],
                od[J]+de[J],io[J],ie[J],seq_len,seq,score[J],gapen[s],grease_pen);
           }
           if(J==start_prof[s+1]) s++;
        }
        char *operation=get_traceback_operation_psm(seq_len,prof_len,nblks,
                block_lengths,start_prof,MAT,DEL,INS,Oper_len,grease_pen, 
                Score, Start, idp);
        free(block_lengths);free(start_prof);
        for(i=0;i<=prof_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        free(MAT);free(DEL);free(INS);
        for(i=0;i<=prof_len;i++) free(score[i]);
        free(score);
        if(grease_pen) free(grease_pen);
        return operation;
}

//============== End Aleksandar's Complex alignment routines ================

char	*psm_typ::Align(e_type E, Int4 rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len) 
{
	char	*operation;
	operation=ComplexAlnSMatrix(idp,LenSeq(E),XSeqPtr(E),rpts,M,Score,Start,Oper_len); 
	// printf("Score=%d\n",*Score);
	return operation;
}

char 	*psm_typ::Align(FILE *fp, e_type E, Int4 rpts, Int4 *Score, Int4 *Start, 
	Int4 *Oper_len) 
{
	char *operation=Align(E, rpts, Score, Start, Oper_len);
	fprintf(fp,"score=%d\n",*Score/40); 
	PutAlign(fp,E,operation,*Oper_len,*Start,rpts);
	return operation;
}

void	psm_typ::ReComputeSMX(cma_typ newcma)
{
	assert(length==TotalLenCMSA(newcma));
	assert(nblk==nBlksCMSA(newcma));
	for(Int4 i=1;i<=nblk; i++) NilSMatrix(M[i]); free(M);;
	delete psg;
        psg = new psg_typ(newcma,PerNats);
	M = psg->BlcksMatrices(MaxRpts);
}

char	*psm_typ::SampleAlign(e_type E, Int4 rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len,
	Int4 wt) 
{
	char	*operation;
	smx_typ *smx = psg->SampleBlcksMatrices(MaxRpts,wt);
#if 0	// Debug...
	if(SeqI(E) == 4971){
	    for(Int4 i=1;i<=nblk; i++) PutSMatrix(stderr, smx[i]);
	    // put_cmaseq_smatrixSW(stderr, operation, Int4 n2, unsigned char *seq2,
            // UInt4 offset, Int4 J, Int4 nmod, smx_typ *M)
	}
#endif
	operation=ComplexAlnSMatrix(idp,LenSeq(E),XSeqPtr(E),rpts,smx,Score,Start,Oper_len); 
	for(Int4 i=1;i<=nblk; i++) NilSMatrix(smx[i]); free(smx);
	return operation;
}

void	psm_typ::Free( )
{
	Int4	i;
	for(i=1;i<=nblk; i++) NilSMatrix(M[i]); free(M);
	free(consensus);
	free(prob_h); free(prob_s); free(prob_c);
	free(lenblk);
	delete idp; // prob_c is copied freed by idp;
	delete psg;
}

void psm_typ::PutAlign(FILE *fp, e_type E1, char *operation, Int4 operation_length, 
               Int4 start, Int4 rpts)
{
        Int4    h, i, total, x, *bl_score, *rpt_score, *space;
        Int4    tmp,e1,e2,j,s,r,t,k,pos;
	Int4	*bl_len;
        e_type  E2 = ConsensusSMatrix(M, nblk, rpts);
        Int4    rpt_len = LenSeq(E2)/rpts;
        Int4    *lengths;
	a_type 	A = SMatrixA(M[1]);
	Int4 	*io = idp->MatToIns();
	Int4	*ie = idp->InsToIns();
	Int4	*od = idp->MatToDel();
	Int4    *de = idp->DelToDel();

	Int4    *grease_pen = idp->SeqPenalties(SeqPtr(E1),LenSeq(E1));
	operation_length-=2;

	NEW(bl_len,nblk+2,Int4);
	for(i=1;i<=nblk;i++) bl_len[i] = idp->Length(i);
        NEW(lengths,nblk*rpts+2,Int4);
	NEW(rpt_score,nblk*rpts+2,Int4);
	NEW(bl_score,nblk*rpts+2,Int4);
	NEW(space,rpts+2,Int4);

        unsigned char *ptr1 = SeqPtr(E1);
        unsigned char *ptr2 = SeqPtr(E2);

        char *seq1, *seq2, *beetw;

        NEW(seq1,operation_length+2,char);
        NEW(seq2,operation_length+2,char);
        NEW(beetw,operation_length+2,char);
        for(k=start,pos=s=j=t=1,i=1; i<=operation_length; i++){
                if (operation[i] == 'M' || operation[i] == 'm'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = AlphaChar(ptr2[j],A);
                        if(ptr1[k] == ptr2[j]) beetw[i] = ':';
                        else if(ValSMatrix(s,ptr1[k],M[t]) > 0) beetw[i] = '.';
                        else beetw[i] = ' ';
			bl_score[t] += ValSMatrix(s,ptr1[k],M[t]);
                        j++;k++;
                        if(s>=LenSMatrix(M[t])) {
				if((t%nblk) == 0) {pos = 1;} else {pos++;}
				s=1; t++;
			} else {s++;pos++;}
                } else if (operation[i] == 'D' || operation[i] == 'd'){
                        seq1[i] = '-'; seq2[i] = AlphaChar(ptr2[j],A); beetw[i] = ' ';
			if (operation[i-1] != 'D' && operation[i-1] != 'd') bl_score[t] += (od[pos] + de[pos]);
                        else bl_score[t] += de[pos];
                        j++;
                        if(s>=LenSMatrix(M[t])) {
				if((t%nblk) == 0) {pos = 1;} else {pos++;}
				s=1;t++;
			} else {s++;pos++;} 
                } else if (operation[i] == 'I'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = '-'; beetw[i] = ' ';
			if (operation[i-1] != 'I'){
				if(grease_pen != NULL) bl_score[t] += (io[pos-1]+ie[pos-1]+grease_pen[k]);
				else bl_score[t] += (io[pos-1]+ie[pos-1]);
			} else {
                                if(grease_pen != NULL) bl_score[t] += ie[pos-1]+grease_pen[k];
                                else bl_score[t] += ie[pos-1];
                        } k++;
                } else if (operation[i] == 'i'){
                        seq1[i] = AlphaChar(ptr1[k],A); seq2[i] = '-'; beetw[i] = ' ';
                        k++;
                } else print_error("PutAlign(): error in operation array");
        }
        fprintf(fp,"\n");
        Int4 br = 0; j=1; i=1; h=1;
        while (i <= operation_length){
                while(br < bl_len[h]){
                        if (seq2[i++] != '-') br++;
                        else lengths[j] += 1;
                }
//fprintf(fp,"lengths[%d]=%d\n",j,lengths[j]);
                br = 0; j++;
		if(h == nblk) h=1; else h++;
        }
	h=1;
	for(total=0,i=1; i<=nblk*rpts; i++){
		x=0;
		for(j=1; j<=bl_len[h] + lengths[i]; j++){
			if (seq2[total+j] == '-') x++;
			else break;
		} space[i] = x;
		total += bl_len[h]+lengths[i];
		if(h == nblk) h=1;
                else h++;
//fprintf(fp,"space[%d]=%d\n",i,space[i]);
	}
	h=1;
	fprintf(fp,"repeat 1:\n");
        for(e1=start,e2=start-1,total=0,i=1; i<=nblk*rpts; i++){
		r=bl_score[i];
		x=space[i]+1;
		fprintf(fp,"blk%2d: ",h);
                for(j=x; j<=bl_len[h] + lengths[i]; j++)
                        fprintf(fp,"%c",seq2[total+j]);
                fprintf(fp,"\n");
		fprintf(fp,"       ");
                for(j=x; j<=bl_len[h] + lengths[i]; j++)
                        fprintf(fp,"%c",beetw[total+j]);
                fprintf(fp,"\n");
		fprintf(fp,"%5d  ",r);
                for(tmp=0,j=x; j<=bl_len[h] + lengths[i]; j++){
                        fprintf(fp,"%c",seq1[total+j]);
                        if (seq1[total+j]!='-') tmp++;
                }
		e1=e2+x; e2+=(x+tmp-1);
		if(i==nblk*rpts) space[nblk*rpts+1] = LenSeq(E1)-e2;
                fprintf(fp," %d-%d (%d)",e1,e2,space[i+1]);
                fprintf(fp,"\n\n");
                total += bl_len[h] + lengths[i];
                if(h == nblk && i!=nblk*rpts) {h=1; fprintf(fp,"repeat %2d:\n",i/nblk+1);}
                else h++;
        }

        free(seq1); free(seq2); free(beetw); 
	free(lengths); free(bl_score); free(space);
	NilSeq(E2);
}

