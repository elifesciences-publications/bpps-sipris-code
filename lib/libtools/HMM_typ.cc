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

#include "HMM_typ.h"

HMM_typ::HMM_typ(Int4 maxrpts,cma_typ cma,double *exp_rpt_gap)
// -G200..1300,20..100/0..50:500,35..250
{ init(maxrpts,"-P200..1300,20..100:500,35..250",200, cma,exp_rpt_gap); }

HMM_typ::HMM_typ(Int4 maxrpts,char *pssm_arg,cma_typ cma,Int4 pernats,double *exp_rpt_gap)
{ init(maxrpts,pssm_arg,pernats,cma,exp_rpt_gap); }

void	HMM_typ::GetProb(Int4 left_flank,Int4 right_flank, cma_typ cma)
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
	// New bayesian update using observed structures...
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
			     default: print_error("HMM_typ::GetProb() input error");
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
}

void	HMM_typ::init(Int4 maxrpts, char *hmm_arg, Int4 pernats, cma_typ cma, double *exp_rpt_gap)
{
	Int4 	n,i,t,left_flank=0,right_flank=0;
	double	*erg=0;

	if(exp_rpt_gap==0){
		NEW(erg,maxrpts+3,double);
		for(i=0;i<=maxrpts;i++){ erg[i] = 1.; }
	}
	target_pvalue=0.05; 
	length=TotalLenCMSA(cma);
	nblk=nBlksCMSA(cma); PerNats = pernats; MaxRpts = maxrpts;
	NEW(lenblk,maxrpts*nblk+3,unsigned short);
	NEW(start_block,maxrpts*nblk+2,Int4);
	NEW(end_block,maxrpts*nblk+2,Int4);
        for(n=0,i=1;i<=maxrpts*nblk;i++){
                start_block[i]=n+1;
		if((t=i%nblk)==0) t=nblk;		// AFN: Fixes one bug!
                lenblk[i] = (unsigned short)LengthCMSA(t,cma);
                n+=lenblk[i];
		end_block[i] = n;
fprintf(stderr,"start_block[%d] = %d; end_block[%d] = %d\n",
		i,start_block[i],i,end_block[i]);
        } 
	psg = new psg_typ(cma,'h','h',0,pernats,'w');
	smx = psg->BlcksMatrices(maxrpts);
	AB=AlphabetCMSA(cma);
	a_type	A = AlphabetCMSA(cma);
	GetProb(left_flank,right_flank,cma);
	double *exp_gap;
	NEW(exp_gap,nblk+3,double);
	for(i=0;i<=nblk;i++) exp_gap[i] = ExpectedGapLengthCMSA(cma,i);
	if(exp_rpt_gap){
	  tpb = new tpb_typ(hmm_arg,'S',prob_c,length,maxrpts,nAlpha(A),nblk,lenblk,
		exp_gap,exp_rpt_gap,pernats);
	} else {
	  tpb = new tpb_typ(hmm_arg,'S',prob_c,length,maxrpts,nAlpha(A),nblk,lenblk,
		exp_gap,erg,pernats); free(erg);
	}
	// tpb->Put(stderr);
	free(exp_gap);
	consensus = ConsensusSeqCMSA(cma); consensus[0]=' ';
	Put(stderr);

        NEWP(mat_emit_prob,maxrpts*length+3,Int4);
        for(i=0;i<=maxrpts*length;i++) NEW(mat_emit_prob[i],nAlpha(AB)+3,Int4);
        for(n=0,i=1;i<=maxrpts*nblk;i++){
           for(Int4 s=1;s<=lenblk[i];s++){
                for(Int4 r=0;r<=nAlpha(AB);r++){
                   mat_emit_prob[n+s][r]=ValSMatrix(s,r,smx[i]);
                }
           } n+=lenblk[i];
	}
}

void	HMM_typ::Put(FILE *fp)
{
	char	c;
	const char head[]=" pos:   m2i  m2d  m2m  i2i  i2m  d2d  d2m  (coil strand helix) SS\n";
	Int4	i,j,b;
	// for(b=1; b <= nblk; b++) PutSMatrix(fp,smx[b]); fprintf(fp,"\n\n");
	fprintf(fp,"Insertion & Deletion penalties:\n");
        fprintf(fp,"blk 1:\n%s",head);
	for(b=1,j=0,i=0;i<=length;i++){
            fprintf(fp,"%c %-3d: %4d %4d %4d %4d %4d %4d %4d  (%.2f  %.2f  %.2f)",consensus[i],i,
		-tpb->MatToIns(i), -tpb->MatToDel(i), -tpb->MatToMat(i),
		-tpb->InsToIns(i), -tpb->InsToMat(i), -tpb->DelToDel(i),
		-tpb->DelToMat(i),prob_c[i],prob_s[i],prob_h[i]);
            if(prob_c[i] > prob_s[i] && prob_c[i] > prob_h[i]) c='C';
            else if(prob_s[i] > prob_c[i] && prob_s[i] > prob_h[i]) c='S';
            else if(prob_h[i] > prob_c[i] && prob_h[i] > prob_s[i]) c='H';
            else c=' ';
            fprintf(fp,"   %c\n",c);
	    if(j==lenblk[b]){ 
		assert(b <= nblk); b++; j=0; 
		if(b <= nblk) fprintf(fp,"\nblk %d:\n%s",b,head); else fprintf(fp,"\n"); 
	    } j++;
        }
}

char	*HMM_typ::FindBestRpts(e_type E, Int4 start_rpts, Int4 *Rpts, 
	Int4 *Score, Int4 *Start, Int4 *OperLen)
// Find the optimum number of repeats...
{
	char	*best_operation=0,*low_operation=0,*high_operation;
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

double  *HMM_typ::ThresholdScore(e_type E, double TargetPvalue, Int4 *TargetScore,
        Int4 *MaxScore)
// compute P-values via Monte Carlo simulation
{
        char    *operation;
        Int4    Score,Start,Oper_len,tmp_score,hpsz,item,N,factor=25;
        Int4    max,s;
        double  pval,*Pvalues;

        assert(TargetPvalue > 0.0 && TargetPvalue < 0.999999);

        e_type  sE = CopySeq(E);
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

char	*HMM_typ::Align(e_type E, Int4 rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len) 
{ return Viterbi(E,rpts,Score,Start,Oper_len); }

char 	*HMM_typ::Align(FILE *fp, e_type E, Int4 rpts, Int4 *Score, Int4 *Start, 
	Int4 *Oper_len) 
{
	char *operation=Viterbi(E, rpts, Score, Start, Oper_len);
	fprintf(fp,"score=%d\n",*Score/40); 
	PutAlign(fp,E,operation,*Oper_len,*Start,rpts);
	return operation;
}

void	HMM_typ::ReComputeSMX(cma_typ newcma)
{
	assert(length==TotalLenCMSA(newcma));
	assert(nblk==nBlksCMSA(newcma));
	for(Int4 i=1;i<=nblk; i++){
		assert(lenblk[i] == LengthCMSA(i,newcma));
		NilSMatrix(smx[i]); 
	} free(smx);;
	delete psg;
        // psg = new psg_typ(newcma,PerNats);
	psg = new psg_typ(newcma,'h','h',0,PerNats,'w');
	smx = psg->BlcksMatrices(MaxRpts);
}

Int4    *HMM_typ::InsEmitProb(unsigned char *seq, Int4 len)
{
	float	*surface=AlphaSurface(AB);
        Int4	*ins_emit;
	
        NEW(ins_emit,len+2,Int4);
        for(Int4 i=1;i<=len;i++){
           ins_emit[i] = (Int4) floor((PerNats*surface[seq[i]])+0.5);
        } return ins_emit;
}

char	*HMM_typ::SampleAlign(e_type E, Int4 Rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len,
	Int4 wt) 
{
	char	*operation=0;
	Int4	n,i;
        smx_typ *tmpM=smx; smx=psg->SampleBlcksMatrices(MaxRpts,wt);
        for(n=0,i=1;i<=MaxRpts*nblk;i++){
           for(Int4 s=1;s<=lenblk[i];s++){
                for(Int4 r=0;r<=nAlpha(AB);r++){
                   mat_emit_prob[n+s][r]=ValSMatrix(s,r,smx[i]);
                }
           } n+=lenblk[i];
	}
	operation=Viterbi(E,Rpts,Score,Start,Oper_len);
	// reset back to normal...
        for(i=1;i<=nblk; i++) NilSMatrix(smx[i]); free(smx); smx=tmpM;
        for(n=0,i=1;i<=MaxRpts*nblk;i++){
           for(Int4 s=1;s<=lenblk[i];s++){
                for(Int4 r=0;r<=nAlpha(AB);r++){
                   mat_emit_prob[n+s][r]=ValSMatrix(s,r,smx[i]);
                }
           } n+=lenblk[i];
	}
	return operation;
}

void	HMM_typ::Free( )
{
	Int4	i;
	for(i=1;i<=nblk; i++) NilSMatrix(smx[i]); free(smx);
	free(consensus);
	free(prob_h); free(prob_s); free(prob_c);
	free(lenblk); free(start_block); free(end_block);
	delete tpb; // prob_c is copied & freed by tpb;
	delete psg;
        for(i=0;i<=MaxRpts*length;i++) free(mat_emit_prob[i]);
        free(mat_emit_prob);
}

void	HMM_typ::inner_loop(Int4 j, Int4 *matj, Int4 *matjm1, Int4 *insj, 
	   Int4 *insjm1, Int4 *delj, Int4 *deljm1,
           Int4 m2m, Int4 d2m, Int4 i2m, Int4  m2d, Int4  d2d, Int4  m2i, Int4 i2i, 
	   Int4 seq_len, unsigned char *seq, Int4 *scorej, Int4 *ins_emit)
// about 80% of compute time is spent here...
{
   register Int4	i,im1,M,I,D;

   for(im1=0,i=1;i<=seq_len;im1++,i++){

	// Find optimum path to match(j,i) state.
	M = matjm1[im1] + m2m; D = deljm1[im1] + d2m; I = insjm1[im1] + i2m;  // = to M transitions
        if(M >= D && M >= I){ matj[i]= M + scorej[seq[i]]; traceM[j][i]='m'; }	// m2m 
        else if(D >= I){ matj[i] = D+scorej[seq[i]]; traceM[j][i]='d'; } 	// d2m
        else { matj[i] = I + scorej[seq[i]]; traceM[j][i]='i'; }		// i2m

	// Find optimum path to delete(j,i) state.			 
	M = matjm1[i] + m2d; D = deljm1[i] + d2d;		// == to D transitions
        if(M > D){ delj[i]=M; traceD[j][i]='m'; } 		// m2d
	else { delj[i] = D; traceD[j][i]='d'; }			// d2d

	// Find optimum path to insert(j,i) state.
	M = matj[im1] + m2i; I = insj[im1] + i2i;		// == to I transitions
        if(M > I){ insj[i]= M + ins_emit[i]; traceI[j][i]='m'; }  	// m2i
	else { insj[i]= I + ins_emit[i]; traceI[j][i]='i'; }		// i2i
   }
}

char	*HMM_typ::Viterbi(e_type E, Int4 Rpts, Int4 *Score, Int4 *Start, Int4 *Oper_len)
{
        Int4    i,j,im1,jm1,seq_len=LenSeq(E);
        Int4    nblks = Rpts*tpb->nBlks();
	Int4	hmm_len=length*Rpts;
         
        assert(Rpts <= tpb->MaxRpts()); // assert(Rpts == 1);
        ins_emit_prob=InsEmitProb(SeqPtr(E),seq_len); 
	Int4 *ins_emit_null; NEW(ins_emit_null,seq_len+3,Int4);
        NEWP(MAT,hmm_len+2,Int4); NEWP(DEL,hmm_len+2,Int4); NEWP(INS,hmm_len+2,Int4);
	NEWP(traceM,hmm_len+2,char); NEWP(traceD,hmm_len+2,char); 
	NEWP(traceI,hmm_len+2,char);
        for(i=0;i<=hmm_len;i++){
                MEW(MAT[i],seq_len+1,Int4); MEW(DEL[i],seq_len+1,Int4);
                MEW(INS[i],seq_len+1,Int4);
		NEW(traceM[i],seq_len+2,char); NEW(traceD[i],seq_len+2,char); 
		NEW(traceI[i],seq_len+2,char);
        }
        DEL[0][0]=tpb->MatToDel(0);	// A kludge; these should be done at DEL[1][0] 
	INS[0][0]=tpb->MatToIns(0);	// set to zero in tpb_typ for now.
	MAT[0][0]=tpb->MatToMat(0);     // set to zero in tpb_typ for now.
	traceM[0][0]=traceI[0][0]=traceD[0][0]='B';
	Int4    *d2d=tpb->DelToDel();
        for(jm1=0,j=1;j<=hmm_len;jm1++,j++) {
		MAT[j][0]=-999999; traceM[j][0] = ' ';	// impossible;
		INS[j][0]=-999999; traceI[j][0] = ' ';	// impossible;
		DEL[j][0]=DEL[jm1][0] + d2d[jm1]; traceD[j][0] = 'd';
        }
	Int4 i2i0 = tpb->InsToIns(0);
	for(im1=0,i=1;i<=seq_len;im1++,i++) {
	   DEL[0][i] = -999999; traceD[0][i] = ' ';  // impossible
	   MAT[0][i] = -999999; traceM[0][i] = ' ';  // impossible
	   INS[0][i] = INS[0][im1] + i2i0; traceI[0][i] = 'i';
	}
	Int4 *ins_emit=ins_emit_prob;
        for(Int4 J=1,Jm1=0;J<=hmm_len;J++,Jm1++){
	   if(J==length) ins_emit=ins_emit_null; else ins_emit=ins_emit_prob;
           inner_loop(J,MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
		tpb->MatToMat(Jm1),tpb->DelToMat(Jm1),tpb->InsToMat(Jm1),
		tpb->MatToDel(Jm1),tpb->DelToDel(Jm1),
		tpb->MatToIns(J),tpb->InsToIns(J),
		seq_len,SeqPtr(E),mat_emit_prob[J],ins_emit);
        } j=hmm_len; i=seq_len;
	if(MAT[j][i] > DEL[j][i] && MAT[j][i] > INS[j][i]) *Score=MAT[j][i];
	else if(DEL[j][i] > INS[j][i]) *Score = DEL[j][i]; else *Score = INS[j][i];
        char *operation=get_traceback(seq_len,hmm_len,Rpts,Oper_len,Start);
        for(i=0;i<=hmm_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        for(i=0;i<=hmm_len;i++){free(traceM[i]); free(traceI[i]); free(traceD[i]); }
	free(traceM); free(traceI); free(traceD); free(MAT);free(DEL);free(INS);
        if(ins_emit_prob) free(ins_emit_prob); free(ins_emit_null);
        return operation;
}

char	*HMM_typ::get_traceback(Int4 seq_len, Int4 hmm_len, Int4 Rpts, 
		Int4 *oper_len, Int4 *start)
/******************************************************************************
  states:
  'i' 'm' and 'd' within blocks
  'M' and 'D' at start of block.
  'I' between blocks.
  'B' at start of block.  'E' at end of block.
 ******************************************************************************/
{
        Int4    i,j,k,bl_nmbr,tot_blks;
        char    state,*operation, *back_operation;

        assert(ins_emit_prob);
        NEW(back_operation,seq_len+hmm_len+3,char);
	tot_blks=nblk*Rpts;
        for(k=0,bl_nmbr=tot_blks,j=hmm_len,i=seq_len,state='E';state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                        back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
                        back_operation[k++]='E'; 
			if(MAT[j][i] > DEL[j][i] && MAT[j][i] > INS[j][i]) state = 'm';
			else if(DEL[j][i] > INS[j][i]){ state = 'd'; }
			else state = 'i';
                        break;
                case 'm': case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; }
			else {
                           if(j != start_block[bl_nmbr]) back_operation[k++]='m';
			   else { back_operation[k++]='M'; bl_nmbr--; }
			   state = traceM[j][i]; j--; i--; 
			} break;
                case 'd': case 'D': // previously sampled deletion
                        if(j <= 1){ back_operation[k++]='D'; state='B'; i++; }
			else {
                             if(j != start_block[bl_nmbr]) back_operation[k++]='d';
			     else { back_operation[k++]='D'; bl_nmbr--; }
			     state = traceD[j][i]; j--;
			} break;
                case 'i':  case 'I': // previously sampled insertion.
                        if(j==end_block[bl_nmbr]) { back_operation[k++]='i'; }
			else { back_operation[k++]='I'; }
			state = traceI[j][i]; i--;
                        break;                   
                default: 
		     fprintf(stderr,"seq_len=%d\n",seq_len);
		     fprintf(stderr,"state = %c; i=%d; j = %d; k = %d\n", state,i,j,k);
		     back_operation[k]=0; std::cerr << back_operation; std::cerr << std::endl;
		     print_error("this should not happen");
           }
        } *oper_len = k;
        NEW(operation,k+4,char);
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
#if 1	// Temporary fix...
	for(i=*oper_len-1; i > 0; i--){
	   if(!(operation[i] == 'i' || operation[i] == 'E')){
		operation[i+1] = 'E'; operation[i+2] = 0; 
		*oper_len = strlen(operation); break;
	   }
	} 
#endif 
        return operation;
}

static	char *probToprbHMM(double *prob, char *operation, Int4 operation_length)
{
	char    *prb;
	Int4	rnd;
        NEW(prb,operation_length+3,char);
	for(Int4 i=1,I=0; i <= operation_length; I++,i++){
	   if(I > 0 && prb[I] !='\n'){
	    switch(operation[i]){
	     case 'M': case 'D': prb[I] = '\n'; I++; break;
	    }
	   } prb[I]= ' ';
	   if(operation[i] == 'm' || operation[i] == 'M'){
		  rnd=(Int4) floor(10.0*prob[i]);
		  switch(rnd){
			case 0: prb[I]=' '; break;
			case 1: prb[I]='1'; break;
			case 2: prb[I]='2'; break;
			case 3: prb[I]='3'; break;
			case 4: prb[I]='4'; break;
			case 5: prb[I]='5'; break;
			case 6: prb[I]='6'; break;
			case 7: prb[I]='7'; break;
			case 8: prb[I]='8'; break;
			default: prb[I]='9'; break;
		  }
	   }
	}
 std::cerr << prb; std::cerr << std::endl;
	return prb;
}

void HMM_typ::PutAlign(FILE *fp, e_type E1, char *operation, Int4 operation_length, 
               Int4 start, Int4 rpts)
{ PutAlign(fp,E1,operation,operation_length,start,rpts,FALSE); }
// { PutAlign(fp,E1,operation,operation_length,start,rpts,TRUE); }

void HMM_typ::PutAlign(FILE *fp, e_type E1, char *operation, Int4 operation_length, 
               Int4 start, Int4 rpts, BooLean use_marg_prob)
{
        Int4    b,i,I,j,s,t,k,pos,*start_sq,*end_sq;
	char	state;

std::cerr << operation; std::cerr << std::endl;
        unsigned char *ptr1 = SeqPtr(E1);
        char *seq1, *seq2, *beetw;
        NEW(seq1,nblk*rpts+operation_length+3,char);
        NEW(seq2,nblk*rpts+operation_length+3,char);
        NEW(beetw,nblk*rpts+operation_length+3,char);
	NEW(start_sq,nblk*rpts+3,Int4); NEW(end_sq,nblk*rpts+3,Int4);
	char	*prb=0;
	if(use_marg_prob){
	  double *prob=MarginalProb(E1,operation, start, operation_length);
	  prb=probToprbHMM(prob,operation,operation_length);
	  free(prob);
	}
	start_sq[1]=start; end_sq[1]=start-1;
        for(b=1,I=0,k=start,pos=s=j=t=1,i=1; i < operation_length; i++){
	   state=operation[i];
	   if(I > 0 && seq1[I] !='\n'){
	    switch(state){
	     case 'M': case 'D':
		seq1[I] = seq2[I]=beetw[I] = '\n'; I++;
		b++; start_sq[b]=k; end_sq[b]=k-1;
	    }
	   }
	   switch(state){
	     case 'M': case 'm':
                 seq1[I] = AlphaChar(ptr1[k],AB); seq2[I] = consensus[pos];
                 if(seq1[I] == seq2[I]) beetw[I] = ':';
                 else if(ValSMatrix(s,ptr1[k],smx[t]) > 0) beetw[I] = '.';
                 else beetw[I] = ' ';
                 j++;k++; I++; end_sq[b]++;
		 if(s>=lenblk[t]) {
			if((t%nblk) == 0) {pos = 1;} else {pos++;}
			s=1; t++;
		 } else {s++;pos++;}
	     break;
	     case 'D': case 'd':
                 seq1[I] = '-'; seq2[I] = consensus[pos]; beetw[I] = ' ';
                 j++; I++;
                 if(s >= lenblk[t]) {
			if((t%nblk) == 0) {pos = 1;} else {pos++;}
			s=1;t++;
		 } else {s++;pos++;} 
	     break;
	     case 'I':
                 seq1[I] = tolower(AlphaChar(ptr1[k],AB));
		 seq2[I] = '-'; beetw[I] = ' ';
                 k++; I++; end_sq[b]++;
	     break;
	     case 'E': case 'B': break;
	     case 'i': k++; break;
	     default: print_error("PutAlign(): error in operation array");
	   }
	} seq1[I] = seq2[I]=beetw[I] = '\n'; I++;
	fprintf(fp,"\n\n");
	PutSeqInfo2(fp,E1);
	fprintf(fp,"\n");
	for(b=1,i=0; i < I; i=j+1,b++){
	    fprintf(fp,"\nblk %-3d ",b);
	    for(j=i; seq2[j] != '\n'; j++) fprintf(fp,"%c",seq2[j]);
	    fprintf(fp,"\n        ");
	    for(j=i; beetw[j] != '\n'; j++) fprintf(fp,"%c",beetw[j]);
	    fprintf(fp,"\nseq     ");
	    for(j=i; seq1[j] != '\n'; j++) fprintf(fp,"%c",seq1[j]);
	    fprintf(fp,"  %d-%d",start_sq[b],end_sq[b]);
	    if(b < nblk) fprintf(fp," (%d)",start_sq[b+1]-end_sq[b]-1);
	    else fprintf(fp," (%d)",LenSeq(E1)-end_sq[b]-1);
	    fprintf(fp,"\n");
	    if(prb){
	       fprintf(fp,"prob    ");
	       for(Int4 x=i; x < j; x++) fprintf(fp,"%c",prb[x]);
	       fprintf(fp,"\n");
	    }
	} fprintf(fp,"\n\n");
	// fprintf(fp,"seq1 %s\n",seq1);
        free(seq1); free(seq2); free(beetw); free(start_sq); free(end_sq);
	if(prb) { free(prb); }
}

Int4	HMM_typ::GetScore(e_type E, char *operation, Int4 Start, Int4 Oper_len)
{
        Int4    i,j,o,im1,jm1,seq_len=LenSeq(E),Score;
	Int4	hmm_len=length*MaxRpts;
	char	state;
	unsigned char *seq=SeqPtr(E);
        Int4	*ins_emit;
         
        // assert(Rpts <= tpb->MaxRpts()); // assert(Rpts == 1);
        ins_emit_prob=InsEmitProb(SeqPtr(E),seq_len); 
	Int4 *ins_emit_null; NEW(ins_emit_null,seq_len+3,Int4);
	Score=0; jm1=-1; j=0; im1=0; o=i=1;
	switch(operation[1]){
	      // case 'i': // BeginToIns; not used yet...
	      case 'M': 		
		if(Start == 1){ state='B'; }
		else { 		// compute score up to start position.
		   Score=tpb->MatToIns(0);	// i.e., BeginToIns
		   for( ; i < Start; i++,im1++) Score+=tpb->InsToIns(0);
		   state='i';
		} break;
	      case 'D': assert(Start==1); state='B'; break;
	      default: print_error("HMM_typ::Score(): error in operation array");
	}
	for( ; o < Oper_len; o++){
		// no ins_emit score between repeats...
		assert(o < (Oper_len) || j <= hmm_len && i <= seq_len);
	        if(j==length) ins_emit=ins_emit_null; else ins_emit=ins_emit_prob;
		switch(operation[o]){
		  case 'M':	// match at start of block
			j++; jm1++; 
			Score+=mat_emit_prob[j][seq[i]];
			switch(state){
		   	  case 'B': Score+=tpb->MatToMat(jm1); 	// i.e., BeginToMatch;
			  case 'i': Score+=tpb->InsToMat(jm1); break;
			  case 'm': Score+=tpb->MatToMat(jm1); break;
			  case 'd': Score+=tpb->DelToMat(jm1); break;
			  default: assert(!"HMM_typ::Score(): operation array error");
			} i++; im1++;
		   break;
		  case 'm':	// match within a block
			j++;jm1++;
			Score+=mat_emit_prob[j][seq[i]];
			switch(state){
			  case 'I': Score+=tpb->InsToMat(jm1); break;
			  case 'm': case 'M': Score+=tpb->MatToMat(jm1); break;
			  case 'd': case 'D': Score+=tpb->DelToMat(jm1); break;
			  default: assert(!"HMM_typ::Score(): operation array error");
			} i++; im1++;
		   break;
		  case 'I':	// insertion within a block
			Score+=ins_emit[i];
			switch(state){
			  case 'I': Score+=tpb->InsToIns(j); break;
			  case 'm': Score+=tpb->MatToIns(j); break;
			  case 'M': Score+=tpb->MatToIns(j); break;
			  default:  
			  fprintf(stderr,"operation = %s\nstate = %c; o = %d\n",
				operation,state,o);
			  assert(!"HMM_typ::Score(): operation array error");
			} i++; im1++;
		   break;
		  case 'i':	// insertion between blocks
			Score+=ins_emit[i];
			switch(state){
			  case 'i': Score+=tpb->InsToIns(j); break;
			  case 'm': Score+=tpb->MatToIns(j); break;
			  default: assert(!"HMM_typ::Score(): operation array error");
			} i++; im1++;
		   break;
		  case 'D':
			j++; jm1++;
			switch(state){
			  case 'B': Score+=tpb->MatToDel(jm1); Score+=tpb->DelToDel(jm1);
			   break;
			  case 'm': Score+=tpb->MatToDel(jm1); break;
			  default:  assert(!"HMM_typ::Score(): operation array error");
			} // j++; jm1++;
		   break;
		  case 'd':
			j++; jm1++;
			switch(state){
			  case 'B': Score+=tpb->MatToDel(jm1); Score+=tpb->DelToDel(jm1);
			    break;
			  case 'M': case 'm': Score+=tpb->MatToDel(jm1); break;
			  case 'D': case 'd': Score+=tpb->DelToDel(jm1); break;
			  default:  assert(!"HMM_typ::Score(): operation array error");
			} // j++; jm1++;
		   break;
		  case 'E':
        		free(ins_emit_prob); free(ins_emit_null);
        		return Score;
		   break;
		  default: assert(!"HMM_typ::Score(): error in operation array");
		} state=operation[o];
        } print_error("HMM_typ::Score(): error in operation array");
	return 0;
}

double	*HMM_typ::MarginalProb0(e_type E, char *operation, Int4 Start, Int4 Oper_len) 
{
        double          *score=0;
        Int4            i,im1,j,jm1,seq_len=LenSeq(E);
        unsigned char   *seq=XSeqPtr(E);
        e_type          rE=CopySeq(E);
        rE = ReverseSeq(rE);
        unsigned char   *rseq=XSeqPtr(rE);
        ins_emit_prob=InsEmitProb(seq,seq_len); 

        double **SMAT,**BMAT,**TMAT,**SDEL,**BDEL,**SINS,**BINS;
        NEWP(SMAT,length+2,double); NEWP(BMAT,length+2,double); NEWP(TMAT,length+2,double);
        NEWP(SDEL,length+2,double); NEWP(BDEL,length+2,double);
        NEWP(SINS,length+2,double); NEWP(BINS,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(SMAT[i],seq_len+1,double); NEW(BMAT[i],seq_len+1,double); NEW(TMAT[i],seq_len+1,double);
                NEW(SDEL[i],seq_len+1,double); NEW(BDEL[i],seq_len+1,double);
                NEW(SINS[i],seq_len+1,double); NEW(BINS[i],seq_len+1,double);
        }
        Int4    *Lm2m=tpb->MatToMat();
        Int4    *Ld2m=tpb->DelToMat();
        Int4    *Li2m=tpb->InsToMat();
        Int4    *Lm2i=tpb->MatToIns();
        Int4    *Li2i=tpb->InsToIns();
        Int4    *Lm2d=tpb->MatToDel();
        Int4    *Ld2d=tpb->DelToDel();
        double  *exp_gap=tpb->ExpBlkGap();
        double  *exp_rpt_gap=tpb->ExpRptGap();
        char    *tpb_arg=tpb->Arg();
        double          *rprob_c,*rexp_gap,*rexp_rpt_gap;
        unsigned short  *rlenblk;
        NEW(rprob_c,length+1,double);
        NEW(rlenblk,nblk+2,unsigned short);
        NEW(rexp_gap,nblk+2,double);
        NEW(rexp_rpt_gap,MaxRpts+3,double);
        for(i=1;i<=length;i++) { rprob_c[i] = prob_c[length-i+1]; }
        for(i=1;i<=nblk;i++) { rlenblk[i] = lenblk[nblk-i+1]; }
        for(i=1;i<=nblk;i++) { rexp_gap[i] = exp_gap[nblk-i+1]; }
        for(i=1;i<=nblk;i++) { rexp_rpt_gap[i] = exp_rpt_gap[MaxRpts-i+1]; }
        tpb_typ *rtpb;
        rtpb = new tpb_typ(tpb_arg,'S',rprob_c,length,MaxRpts,nAlpha(AB),nblk,rlenblk,
                rexp_gap,rexp_rpt_gap,PerNats);
        Int4 *Lrm2m=rtpb->MatToMat();
        Int4 *Lrd2m=rtpb->DelToMat();
        Int4 *Lri2m=rtpb->InsToMat();
        Int4 *Lrm2i=rtpb->MatToIns();
        Int4 *Lri2i=rtpb->InsToIns();
        Int4 *Lrm2d=rtpb->MatToDel();
        Int4 *Lrd2d=rtpb->DelToDel();

        double *m2m,*d2m,*i2m,*m2i,*i2i,*m2d,*d2d;
        double *rm2m,*rd2m,*ri2m,*rm2i,*ri2i,*rm2d,*rd2d;
        double *rins_emit_prob,**rmat_emit_prob;
        double *pins_emit_prob, **pmat_emit_prob;

        NEW(m2m,length+2,double);NEW(d2m,length+2,double);NEW(i2m,length+2,double);
        NEW(m2i,length+2,double);NEW(i2i,length+2,double);NEW(m2d,length+2,double);NEW(d2d,length+2,double);
        NEW(rm2m,length+2,double);NEW(rd2m,length+2,double);NEW(ri2m,length+2,double);
        NEW(rm2i,length+2,double);NEW(ri2i,length+2,double);NEW(rm2d,length+2,double);NEW(rd2d,length+2,double);
        NEW(pins_emit_prob,seq_len+2,double);NEWP(pmat_emit_prob,length+2,double);
        NEW(rins_emit_prob,seq_len+2,double);NEWP(rmat_emit_prob,length+2,double);
        for(i=1;i<=length;i++){
                NEW(pmat_emit_prob[i],nAlpha(AB)+1,double);
                NEW(rmat_emit_prob[i],nAlpha(AB)+1,double);
        }

        pmat_emit_prob = psg->TargetFreq();

        for(i=0;i<=length;i++){
                m2m[i] = exp((double) Lm2m[i]/PerNats);
                d2m[i] = exp((double) Ld2m[i]/PerNats);
                i2m[i] = exp((double) Li2m[i]/PerNats);
                m2i[i] = exp((double) Lm2i[i]/PerNats);
                i2i[i] = exp((double) Li2i[i]/PerNats);
                m2d[i] = exp((double) Lm2d[i]/PerNats);
                d2d[i] = exp((double) Ld2d[i]/PerNats);

                rm2m[i] = exp((double) Lrm2m[i]/PerNats);
                rd2m[i] = exp((double) Lrd2m[i]/PerNats);
                ri2m[i] = exp((double) Lri2m[i]/PerNats);
                rm2i[i] = exp((double) Lrm2i[i]/PerNats);
                ri2i[i] = exp((double) Lri2i[i]/PerNats);
                rm2d[i] = exp((double) Lrm2d[i]/PerNats);
                rd2d[i] = exp((double) Lrd2d[i]/PerNats);
        }
        for(i=1;i<=seq_len;i++){
                rins_emit_prob[i] = exp((double) ins_emit_prob[seq_len-i+1]/PerNats);
                pins_emit_prob[i] = exp((double) ins_emit_prob[i]/PerNats);
        }
        for(i=1;i<=length;i++){
        //      for(j=0;j<=nAlpha(AB);j++){
rmat_emit_prob[i] = pmat_emit_prob[length-i+1];//exp((double) mat_emit_prob[length-i+1][j]/PerNats) * blosum62freq[j];
        //              pmat_emit_prob[i][j] = exp((double) mat_emit_prob[i][j]/PerNats) * blosum62freq[j];
        //      }
        }
        SDEL[1][1] = m2d[0]; SMAT[1][1] = m2m[0] * pmat_emit_prob[1][seq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                SINS[j][1] = 0; SDEL[j][1] = SDEL[jm1][1] * d2d[jm1];
                SMAT[j][1] = SDEL[jm1][1] * d2m[jm1] * pmat_emit_prob[j][seq[1]];
        }
        SINS[1][1] = 0; SMAT[1][2] = m2i[0] * pmat_emit_prob[1][seq[2]];
        SINS[1][2] = m2m[0] * pmat_emit_prob[1][seq[1]] * m2i[1] * pins_emit_prob[2];
        SDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                SDEL[1][i] = 0; SMAT[1][i] = m2i[0] * pow(i2i[0],i-2) * pmat_emit_prob[1][seq[i]];
                SINS[1][i] = (SINS[1][im1] * i2i[1] * pins_emit_prob[i]) + (SMAT[1][im1] * m2i[1] * pins_emit_prob[i]);
        }
        BDEL[1][1] = rm2d[0]; BMAT[1][1] = rm2m[0] * rmat_emit_prob[1][rseq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                BINS[j][1] = 0; BDEL[j][1] = BDEL[jm1][1] * rd2d[jm1];
                BMAT[j][1] = BDEL[jm1][1] * rd2m[jm1] * rmat_emit_prob[j][rseq[1]];
        }                        
        BINS[1][1] = 0; BMAT[1][2] = rm2i[0] * rmat_emit_prob[1][rseq[2]];
        BINS[1][2] = rm2m[0] * rmat_emit_prob[1][rseq[1]] * rm2i[1] * rins_emit_prob[2];
        BDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                BDEL[1][i] = 0; BMAT[1][i] = m2i[0] * pow(i2i[0],i-2) * rmat_emit_prob[1][rseq[i]];
                BINS[1][i] = (BINS[1][im1] * ri2i[1] * rins_emit_prob[i]) + (BMAT[1][im1] * rm2i[1] * rins_emit_prob[i]);
        }
        double *smatj,*sdelj,*sinsj,*bmatj,*bdelj,*binsj;
        double *smatjm1,*sdeljm1,*sinsjm1,*bmatjm1,*bdeljm1,*binsjm1;
        double *pmat_emit_probj,*rmat_emit_probj,pmat_emit_probjres,rmat_emit_probjres;
        double pins_emit_probi,rins_emit_probi;
        double m2m_jm1,d2m_jm1,i2m_jm1,m2d_jm1,d2d_jm1,m2i_j,i2i_j;
        double rm2m_jm1,rd2m_jm1,ri2m_jm1,rm2d_jm1,rd2d_jm1,rm2i_j,ri2i_j;
        unsigned char res,rres;
        for(jm1=1,j=2;j<=length;j++,jm1++){
                smatj=SMAT[j];sdelj=SDEL[j];sinsj=SINS[j];
                bmatj=BMAT[j];bdelj=BDEL[j];binsj=BINS[j];
                smatjm1=SMAT[jm1];sdeljm1=SDEL[jm1];sinsjm1=SINS[jm1];
                bmatjm1=BMAT[jm1];bdeljm1=BDEL[jm1];binsjm1=BINS[jm1];
                pmat_emit_probj=pmat_emit_prob[j];
                rmat_emit_probj=rmat_emit_prob[j];
                m2m_jm1=m2m[jm1];d2m_jm1=d2m[jm1];i2m_jm1=i2m[jm1];
                m2d_jm1=m2d[jm1];d2d_jm1=d2d[jm1];m2i_j=m2i[j];i2i_j=i2i[j];
                rm2m_jm1=rm2m[jm1];rd2m_jm1=rd2m[jm1];ri2m_jm1=ri2m[jm1];
                rm2d_jm1=rm2d[jm1];rd2d_jm1=rd2d[jm1];rm2i_j=rm2i[j];ri2i_j=ri2i[j];
                for(im1=1,i=2;i<=seq_len;im1++,i++){
                        res=seq[i]; rres=rseq[i];
                        pmat_emit_probjres=pmat_emit_probj[res]; rmat_emit_probjres=rmat_emit_probj[rres];
                        pins_emit_probi=pins_emit_prob[i];rins_emit_probi=rins_emit_prob[i];
                        smatj[i] += (smatjm1[im1] * m2m_jm1 * pmat_emit_probjres);
                        smatj[i] += (sdeljm1[im1] * d2m_jm1 * pmat_emit_probjres); 
                        smatj[i] += (sinsjm1[im1] * i2m_jm1 * pmat_emit_probjres);
                        sdelj[i] += (smatjm1[i] * m2d_jm1);
                        sdelj[i] += (sdeljm1[i] * d2d_jm1);
                        sinsj[i] += (smatj[im1] * m2i_j * pins_emit_probi);
                        sinsj[i] += (sinsj[im1] * i2i_j * pins_emit_probi);
                        bmatj[i] += (bmatjm1[im1] * rm2m_jm1 * rmat_emit_probjres);
                        bmatj[i] += (bdeljm1[im1] * rd2m_jm1 * rmat_emit_probjres);
                        bmatj[i] += (binsjm1[im1] * ri2m_jm1 * rmat_emit_probjres);
                        bdelj[i] += (bmatjm1[i] * rm2d_jm1);
                        bdelj[i] += (bdeljm1[i] * rd2d_jm1);
                        binsj[i] += (bmatj[im1] * rm2i_j * rins_emit_probi);
                        binsj[i] += (binsj[im1] * ri2i_j * rins_emit_probi);
                }
        }
        double *tmatj,*bmatlenmjp1,sum=0.;
        for(j=1;j<=length;j++){
                tmatj=TMAT[j];smatj=SMAT[j];bmatlenmjp1=BMAT[length-j+1];pmat_emit_probj=pmat_emit_prob[j];
                for(i=1;i<=seq_len;i++){
                        res=seq[i];
                        tmatj[i]=bmatlenmjp1[seq_len-i+1]*(smatj[i]/pmat_emit_probj[res]);
//printf("s[%d][%d]=%g b[%d][%d]=%g t[%d][%d]=%g\n",j,i,SMAT[j][i],j,i,BMAT[length-j+1][seq_len-i+1],j,i,TMAT[j][i]);
                        sum+=tmatj[i];
                }

        }
        double last_cell=TMAT[length][seq_len];
        for(j=1;j<=length;j++){
//              tmatj=TMAT[j];
                for(i=1;i<=seq_len;i++){
                        TMAT[j][i]/=last_cell;;
printf("TMAT[%d][%d]=%g\n",j,i,TMAT[j][i]);
                }
        }

        NEW(score,Oper_len+2,double);
        j=1; 
        Int4 posProf=1, posSeq=Start;
        for(i=1;i<=Oper_len-2;i++){
                switch(operation[i]){
                        case 'M':
                        case 'm': score[j++] = TMAT[posProf++][posSeq++]; break; 
                        case 'D': 
                        case 'd': j++; posProf++; break;
                        case 'I':
                        case 'i': j++; posSeq++; break; 
                        default : print_error("MarginalProb(): error in the operation array");
                }
        }
        for(i=1;i<=Oper_len-2;i++) printf("score[%d]=%g\n",i,score[i]);
        for(i=0;i<=length+1;i++) { 
                free(SMAT[i]);free(SDEL[i]);free(SINS[i]);
                free(BMAT[i]);free(BDEL[i]);free(BINS[i]);free(TMAT[i]);
        }
        free(SMAT);free(SDEL);free(SINS);free(BMAT);free(BDEL);free(BINS);free(TMAT);
        free(rprob_c);free(rlenblk);free(rexp_gap);free(rexp_rpt_gap);
        free(m2m);free(d2m);free(i2m);free(m2i);free(i2i);free(m2d);free(d2d);
        free(rm2m);free(rd2m);free(ri2m);free(rm2i);free(ri2i);free(rm2d);free(rd2d);
        free(pins_emit_prob);free(rins_emit_prob);
        for(i=1;i<=length;i++) { free(pmat_emit_prob[i]); free(rmat_emit_prob[i]); }
        free(pmat_emit_prob);free(rmat_emit_prob);
        return score;
}

double	*HMM_typ::MarginalProb(e_type E, char *operation, Int4 Start, Int4 Oper_len) 
{
	double		*score=0;
	Int4		i,im1,j,jm1,seq_len=LenSeq(E);
        unsigned char   *seq=XSeqPtr(E);
        e_type          rE=CopySeq(E);
        rE = ReverseSeq(rE);
        unsigned char   *rseq=XSeqPtr(rE);
	ins_emit_prob=InsEmitProb(seq,seq_len); 
	double **SMAT,**BMAT,**TMAT,**SDEL,**BDEL,**SINS,**BINS;
        NEWP(SMAT,length+2,double); NEWP(BMAT,length+2,double); NEWP(TMAT,length+2,double);
        NEWP(SDEL,length+2,double); NEWP(BDEL,length+2,double);
        NEWP(SINS,length+2,double); NEWP(BINS,length+2,double);
        for(i=0;i<=length+1;i++){
                NEW(SMAT[i],seq_len+1,double); NEW(BMAT[i],seq_len+1,double); NEW(TMAT[i],seq_len+1,double);
                NEW(SDEL[i],seq_len+1,double); NEW(BDEL[i],seq_len+1,double);
                NEW(SINS[i],seq_len+1,double); NEW(BINS[i],seq_len+1,double);
        }
        Int4 		*Lm2m=tpb->MatToMat();
        Int4 		*Ld2m=tpb->DelToMat();
        Int4 		*Li2m=tpb->InsToMat();
        Int4 		*Lm2i=tpb->MatToIns();
        Int4 		*Li2i=tpb->InsToIns();
        Int4 		*Lm2d=tpb->MatToDel();
        Int4 		*Ld2d=tpb->DelToDel();
	double 		*exp_gap=tpb->ExpBlkGap();
	double 		*exp_rpt_gap=tpb->ExpRptGap();
	char 		*tpb_arg=tpb->Arg();
        double  	*rprob_c,*rexp_gap,*rexp_rpt_gap;
	unsigned short 	*rlenblk;

        NEW(rprob_c,length+1,double);
        NEW(rlenblk,nblk+2,unsigned short);
        NEW(rexp_gap,nblk+2,double);
        NEW(rexp_rpt_gap,MaxRpts+3,double);
        for(i=1;i<=length;i++) { rprob_c[i] = prob_c[length-i+1]; }
        for(i=1;i<=nblk;i++) { rlenblk[i] = lenblk[nblk-i+1]; }
        for(i=1;i<=nblk;i++) { rexp_gap[i] = exp_gap[nblk-i+1]; }
        for(i=1;i<=nblk;i++) { rexp_rpt_gap[i] = exp_rpt_gap[MaxRpts-i+1]; }
        tpb_typ *rtpb;
        rtpb = new tpb_typ(tpb_arg,'S',rprob_c,length,MaxRpts,nAlpha(AB),nblk,rlenblk,
                rexp_gap,rexp_rpt_gap,PerNats);
        Int4 *Lrm2m=rtpb->MatToMat();
        Int4 *Lrd2m=rtpb->DelToMat();
        Int4 *Lri2m=rtpb->InsToMat();
        Int4 *Lrm2i=rtpb->MatToIns();
        Int4 *Lri2i=rtpb->InsToIns();
        Int4 *Lrm2d=rtpb->MatToDel();
        Int4 *Lrd2d=rtpb->DelToDel();
	double *m2m,*d2m,*i2m,*m2i,*i2i,*m2d,*d2d;
	double *rm2m,*rd2m,*ri2m,*rm2i,*ri2i,*rm2d,*rd2d;
	double *rins_emit_prob,**rmat_emit_prob;
	double *pins_emit_prob, **pmat_emit_prob;	
	NEW(m2m,length+2,double);NEW(d2m,length+2,double);NEW(i2m,length+2,double);
	NEW(m2i,length+2,double);NEW(i2i,length+2,double);NEW(m2d,length+2,double);NEW(d2d,length+2,double);
	NEW(rm2m,length+2,double);NEW(rd2m,length+2,double);NEW(ri2m,length+2,double);
	NEW(rm2i,length+2,double);NEW(ri2i,length+2,double);NEW(rm2d,length+2,double);NEW(rd2d,length+2,double);
	NEW(pins_emit_prob,seq_len+2,double);
	NEW(rins_emit_prob,seq_len+2,double);NEWP(rmat_emit_prob,length+2,double);
	pmat_emit_prob = psg->TargetFreq();
	for(i=0;i<=length;i++){
		m2m[i] = exp((double) Lm2m[i]/PerNats);
		d2m[i] = exp((double) Ld2m[i]/PerNats);
		i2m[i] = exp((double) Li2m[i]/PerNats);
		m2i[i] = exp((double) Lm2i[i]/PerNats);
		i2i[i] = exp((double) Li2i[i]/PerNats);
		m2d[i] = exp((double) Lm2d[i]/PerNats);
		d2d[i] = exp((double) Ld2d[i]/PerNats);
		rm2m[i] = exp((double) Lrm2m[i]/PerNats);
		rd2m[i] = exp((double) Lrd2m[i]/PerNats);
		ri2m[i] = exp((double) Lri2m[i]/PerNats);
		rm2i[i] = exp((double) Lrm2i[i]/PerNats);
		ri2i[i] = exp((double) Lri2i[i]/PerNats);
		rm2d[i] = exp((double) Lrm2d[i]/PerNats);
		rd2d[i] = exp((double) Lrd2d[i]/PerNats);
	}
	for(i=1;i<=seq_len;i++){
		rins_emit_prob[i] = exp((double) ins_emit_prob[seq_len-i+1]/PerNats) * blosum62freq[rseq[i]];
		pins_emit_prob[i] = exp((double) ins_emit_prob[i]/PerNats) * blosum62freq[seq[i]];
	}
	for(i=1;i<=length;i++){
		rmat_emit_prob[i] = pmat_emit_prob[length-i+1];
	}

        SDEL[1][1] = m2d[0]; SMAT[1][1] = m2m[0] * pmat_emit_prob[1][seq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                SINS[j][1] = 0; SDEL[j][1] = SDEL[jm1][1] * d2d[jm1];
                SMAT[j][1] = SDEL[jm1][1] * d2m[jm1] * pmat_emit_prob[j][seq[1]];
        }
        SINS[1][1] = 0; SMAT[1][2] = m2i[0] * pmat_emit_prob[1][seq[2]];
        SINS[1][2] = m2m[0] * pmat_emit_prob[1][seq[1]] * m2i[1] * pins_emit_prob[2];
        SDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                SDEL[1][i] = 0; SMAT[1][i] = m2i[0] * pow(i2i[0],i-2) * pmat_emit_prob[1][seq[i]];
                SINS[1][i] = (SINS[1][im1] * i2i[1] * pins_emit_prob[i]) + (SMAT[1][im1] * m2i[1] * pins_emit_prob[i]);
        }


        BDEL[1][1] = rm2d[0]; BMAT[1][1] = rm2m[0] * rmat_emit_prob[1][rseq[1]];
        for(jm1=1,j=2;j<=length;jm1++,j++) {
                BINS[j][1] = 0; BDEL[j][1] = BDEL[jm1][1] * rd2d[jm1];
                BMAT[j][1] = BDEL[jm1][1] * rd2m[jm1] * rmat_emit_prob[j][rseq[1]];
        }                        
        BINS[1][1] = 0; BMAT[1][2] = rm2i[0] * rmat_emit_prob[1][rseq[2]];
        BINS[1][2] = rm2m[0] * rmat_emit_prob[1][rseq[1]] * rm2i[1] * rins_emit_prob[2];
        BDEL[1][2] = 0;
        for(im1=2,i=3;i<=seq_len;im1++,i++) {
                BDEL[1][i] = 0; BMAT[1][i] = rm2i[0] * pow(ri2i[0],i-2) * rmat_emit_prob[1][rseq[i]];
                BINS[1][i] = (BINS[1][im1] * ri2i[1] * rins_emit_prob[i]) + (BMAT[1][im1] * rm2i[1] * rins_emit_prob[i]);
        }



        double *smatj,*sdelj,*sinsj,*bmatj,*bdelj,*binsj;
        double *smatjm1,*sdeljm1,*sinsjm1,*bmatjm1,*bdeljm1,*binsjm1;
        double *pmat_emit_probj,*rmat_emit_probj,pmat_emit_probjres,rmat_emit_probjres;
	double pins_emit_probi,rins_emit_probi;
        double m2m_jm1,d2m_jm1,i2m_jm1,m2d_jm1,d2d_jm1,m2i_j,i2i_j;
        double rm2m_jm1,rd2m_jm1,ri2m_jm1,rm2d_jm1,rd2d_jm1,rm2i_j,ri2i_j;
	unsigned char res,rres;
        for(jm1=1,j=2;j<=length;j++,jm1++){
		smatj=SMAT[j];sdelj=SDEL[j];sinsj=SINS[j];
		bmatj=BMAT[j];bdelj=BDEL[j];binsj=BINS[j];
		smatjm1=SMAT[jm1];sdeljm1=SDEL[jm1];sinsjm1=SINS[jm1];
		bmatjm1=BMAT[jm1];bdeljm1=BDEL[jm1];binsjm1=BINS[jm1];
		pmat_emit_probj=pmat_emit_prob[j];
		rmat_emit_probj=rmat_emit_prob[j];
		m2m_jm1=m2m[jm1];d2m_jm1=d2m[jm1];i2m_jm1=i2m[jm1];
		m2d_jm1=m2d[jm1];d2d_jm1=d2d[jm1];m2i_j=m2i[j];i2i_j=i2i[j];
		rm2m_jm1=rm2m[jm1];rd2m_jm1=rd2m[jm1];ri2m_jm1=ri2m[jm1];
		rm2d_jm1=rm2d[jm1];rd2d_jm1=rd2d[jm1];rm2i_j=rm2i[j];ri2i_j=ri2i[j];
                for(im1=1,i=2;i<=seq_len;im1++,i++){
			res=seq[i]; rres=rseq[i];
			pmat_emit_probjres=pmat_emit_probj[res]; rmat_emit_probjres=rmat_emit_probj[rres];
			pins_emit_probi=pins_emit_prob[i];rins_emit_probi=rins_emit_prob[i];
                        smatj[i] += (smatjm1[im1] * m2m_jm1 * pmat_emit_probjres);
                        smatj[i] += (sdeljm1[im1] * d2m_jm1 * pmat_emit_probjres); 
                        smatj[i] += (sinsjm1[im1] * i2m_jm1 * pmat_emit_probjres);
                        sdelj[i] += (smatjm1[i] * m2d_jm1);
                        sdelj[i] += (sdeljm1[i] * d2d_jm1);
                        sinsj[i] += (smatj[im1] * m2i_j * pins_emit_probi);
                        sinsj[i] += (sinsj[im1] * i2i_j * pins_emit_probi);
                        bmatj[i] += (bmatjm1[im1] * rm2m_jm1 * rmat_emit_probjres);
                        bmatj[i] += (bdeljm1[im1] * rd2m_jm1 * rmat_emit_probjres);
                        bmatj[i] += (binsjm1[im1] * ri2m_jm1 * rmat_emit_probjres);
                        bdelj[i] += (bmatjm1[i] * rm2d_jm1);
                        bdelj[i] += (bdeljm1[i] * rd2d_jm1);
                        binsj[i] += (bmatj[im1] * rm2i_j * rins_emit_probi);
                        binsj[i] += (binsj[im1] * ri2i_j * rins_emit_probi);
                }
        }
	for(j=1;j<=length;j++){
		for(i=1;i<=seq_len;i++){
			res=seq[i];
			TMAT[j][i]=BMAT[length-j+1][seq_len-i+1]*(SMAT[j][i]/pmat_emit_prob[j][res]);
		}		
	}

double tdel=SDEL[length][seq_len];
double tins=SINS[length][seq_len];
printf("TMAT=%g  TDEL=%g  TINS=%g\n",TMAT[length][seq_len],tdel,tins);

	double last_cell=TMAT[length][seq_len]+tdel+tins;
	for(j=1;j<=length;j++){
		for(i=1;i<=seq_len;i++){
			TMAT[j][i]/=last_cell;;
printf("TMAT[%d][%d]=%g\n",j,i,TMAT[j][i]);
		}
	}

	NEW(score,Oper_len+2,double);
	j=1; 
	Int4 posProf=1, posSeq=Start;
	for(i=1;i<=Oper_len-2;i++){
		switch(operation[i]){
			case 'M':
			case 'm': score[j++] = TMAT[posProf++][posSeq++]; break; 
			case 'D': 
			case 'd': j++; posProf++; break;
			case 'I':
			case 'i': j++; posSeq++; break; 
			default : print_error("MarginalProb(): error in the operation array");
		}
	}
	for(i=1;i<=Oper_len-2;i++) printf("score[%d]=%g\n",i,score[i]);
	for(i=0;i<=length+1;i++) { 
		free(SMAT[i]);free(SDEL[i]);free(SINS[i]);
		free(BMAT[i]);free(BDEL[i]);free(BINS[i]);free(TMAT[i]);
	}
	free(SMAT);free(SDEL);free(SINS);free(BMAT);free(BDEL);free(BINS);free(TMAT);
        free(rprob_c);free(rlenblk);free(rexp_gap);free(rexp_rpt_gap);
	free(m2m);free(d2m);free(i2m);free(m2i);free(i2i);free(m2d);free(d2d);
        free(rm2m);free(rd2m);free(ri2m);free(rm2i);free(ri2i);free(rm2d);free(rd2d);
        free(pins_emit_prob);free(rins_emit_prob);
	for(i=1;i<=length;i++) { free(pmat_emit_prob[i]); free(rmat_emit_prob[i]); }
	free(pmat_emit_prob);free(rmat_emit_prob);

	return score;
}

//============== End Aleksandar's Complex alignment routines ================


