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

#include "smatrix.h"

/************************* dynamic programming **************************/
Int4	FastAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *X[2])
/*******************************************************************
 Performs a fast dynamic programming alignment for global score.
********************************************************************/
{
	Int4	r,m,s,lenM,start,end,j,**scrmtx;
	Int4	*Dm1,*Tmp,*D,Fm1,F,K;
	unsigned char *sq;

	D = X[0]; Dm1 = X[1];
	for(r=0; r<=len; r++) Dm1[r]=0;
        for(start=0,m=1; m<=nmod; m++) start+=LenSMatrix(M[m]);
	end = len-start+1; seq--;
	for(lenM=0,start=1,m=1; m<=nmod; m++) {
	   D[start-1] = Fm1 = INT4_MIN;
	   scrmtx = M[m]->score; K = M[m]->K;
	   for(sq=seq+start,r=start; r<=end; r++,sq++) {
		s=Dm1[r-lenM];
		for(j=K; j; j--) s+= scrmtx[j][sq[j]]; 
		F = D[r-1];
		if(F < Fm1) F = Fm1; else Fm1=F;
                if(s < F) s = F;
		D[r] = s; 
	   }
	   lenM = LenSMatrix(M[m]);
	   start+=lenM; end += lenM;
	   Tmp=Dm1; Dm1=D; D=Tmp; 
	}
	return s;
}

Int4	FastLocalAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *X[2])
/*******************************************************************
 Perform a marginally faster dynamic programming alignment for local score
********************************************************************/
{
	Int4	r,m,s,lenM,start,end,j,**scrmtx;
	Int4	*Dm1,*Tmp,*D,Fm1,F,K,max;
	unsigned char *sq;

	D = X[0]; Dm1 = X[1];
	for(r=0; r<=len; r++) Dm1[r]=0;
        for(start=0,m=1; m<=nmod; m++) start+=LenSMatrix(M[m]);
	end = len-start+1; seq--;
	for(lenM=0,start=1,m=1; m<=nmod; m++) {
	   D[start-1] = Fm1 = INT4_MIN;
	   scrmtx = M[m]->score; K = M[m]->K;
	   for(sq=seq+start,r=start; r<=end; r++,sq++) {
		for(s=max=0,j=K; j; j--) {
		    s += scrmtx[j][sq[j]]; 
		    if(s < 0) s=0; else if(s > max) max = s;
		} 
		s=Dm1[r-lenM]+max;
		F = D[r-1];
		if(F < Fm1) F = Fm1; else Fm1=F;
                if(s < F) s = F;
		D[r] = s; 
	   }
	   lenM = LenSMatrix(M[m]);
	   start+=lenM; end += lenM;
	   Tmp=Dm1; Dm1=D; D=Tmp; 
	}
	return s;
}

void	put_dp_matrixSMX(Int4 len, smx_typ M, unsigned char *profile, 
	unsigned char *seq, a_type A, Int4 startM,Int4 start, Int4 end)
{
	Int4	r,m,i,lenM;
	
	end = MINIMUM(Int4,len,end);
	start = MAXIMUM(Int4,start,0);
	end = MAXIMUM(Int4,start,end);
	start = MINIMUM(Int4,start,end);
	lenM=LenSMatrix(M);
	fprintf(stderr,"     |");
	for(i=startM,m=1; m<=lenM; m++,i++) fprintf(stderr,"  %c  |",
				AlphaChar(profile[i],A));
	fprintf(stderr,"\n     |");
	for(m=1; m<=lenM; m++) fprintf(stderr,"%3d  |",m);
	fprintf(stderr,"\n");
	for(r=start; r<=end; r++) {	
	   fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
	   for(i=startM,m=1; m<=lenM; m++,i++) {
		fprintf(stderr,"%5d|",ValSMatrix(m,seq[r],M));
	   }
	   fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
}

static void	put_dp_matrixSW(Int4 **D, Int4** T, Int4 len, smx_typ M,
	unsigned char *profile, unsigned char *seq, a_type A,
	Int4 startM,Int4 start, Int4 end)
{
	Int4	r,m,i,lenM;
	
	end = MINIMUM(Int4,len,end);
	start = MAXIMUM(Int4,start,0);
	end = MAXIMUM(Int4,start,end);
	start = MINIMUM(Int4,start,end);
	lenM=LenSMatrix(M);
	fprintf(stderr,"     |");
	for(i=startM,m=1; m<=lenM; m++,i++) fprintf(stderr,"  %c  |",
				AlphaChar(profile[i],A));
	fprintf(stderr,"\n     |");
	for(m=1; m<=lenM; m++) fprintf(stderr,"%3d  |",m);
	fprintf(stderr,"\n");
	for(r=start; r<=end; r++) {	
	   fprintf(stderr,"%c%4d|",AlphaChar(seq[r],A),r);
	   for(i=startM,m=1; m<=lenM; m++,i++) {
		if(D[i][r] == SHRT_MIN) fprintf(stderr," -inf ");
		else fprintf(stderr,"%5d",D[i][r]);
		if(T[i][r] < 0) fprintf(stderr,"^");
		else if(T[i][r] == 0) fprintf(stderr,"\\");
		else fprintf(stderr,"<");
	   }
	   fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
}

static void     put_dp_matrix0(double **D, Int4 len, Int4 prolen)
{
   Int4    r,m,v,start,end;

   for(start=0,end=14; start <= prolen; start+=15,end+=15) {
        fprintf(stderr,"     |");
        for(m=start; m<=end && m <= prolen; m++) { fprintf(stderr,"%3d  |",m); }
        fprintf(stderr,"\n");
        for(r=0; r<=len; r++) {
           fprintf(stderr,"%5d|",r);
           for(m=start; m<=end && m <= prolen; m++) {
                if(D[r][m] == 0) fprintf(stderr,"  -  ");
                else {
                        v = (Int4) log10(D[r][m]);
                        fprintf(stderr,"%5d",v);
                }
                fprintf(stderr,"|");
           }
           fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
   }
}

static void	put_dp_matrix(Int4 **D, BooLean** T, Int4 len, Int4 nmod, smx_typ *M)
{
	Int4	r,m,lenM,start,end;

	for(m=1; m<=nmod; m++) fprintf(stderr,"len[%d] = %d\n",
			m,LenSMatrix(M[m]));
	fprintf(stderr,"     |");
	for(m=0; m<=nmod; m++) {
		fprintf(stderr,"%4d  |",m);
	}
	fprintf(stderr,"\n");
	for(r=0; r<=len; r++) {	
	   fprintf(stderr,"%5d|",r);
	   for(m=0; m<=nmod; m++) {
		if(D[m][r] == INT4_MIN) fprintf(stderr," -inf ");
		else fprintf(stderr,"%6d",D[m][r]);
		if(T[m][r]) fprintf(stderr,"<");
		else fprintf(stderr,"|");
	   }
	   fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
}

static Int4 fast_aln_seq_overlap_smatrix(Int4 len, unsigned char *seq,
	Int4 nmod, smx_typ *M, char *overlap,
	Int4 (*score_smtrx)(register unsigned char *, register Int4, 
	register Int4 **))
// alignment with Overlaps
// len == len sequence...
{
	Int4	r,m,s,lenM,score,start,end,**D,**F,scr;
	Int4	**scrmtx,K;
	unsigned char *sq;
	double	p;
	Int4 begin=-overlap[0];

	// 0. Make sure that sequence is Int4 enough to accommodate motifs
	for(s=0,m=1; m<=nmod; m++){
		s+=(LenSMatrix(M[m])-overlap[m]);
	} if(s > len) return SHRT_MIN;
	
	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(F,nmod+3,Int4*);
	for(m=0; m<= nmod; m++) { 
	   MEW(D[m],len+overlap[0]+3,Int4); D[m]-=begin;
	   MEW(F[m],len+overlap[0]+3,Int4); F[m]-=begin;
	}
	/*** 2. Dynamic programming step. ***/
	for(r=begin; r<=len; r++) { D[0][r] = 0; }
	for(start=begin,m=1; m<=nmod; m++){ 
	   D[m][start]=INT4_MIN; F[m][start]=INT4_MIN; 
	   start+= (LenSMatrix(M[m])-overlap[m]);
	} 
	end = len-start+1-overlap[0]; score=SHRT_MIN;
	sq = seq-1;
	for(lenM=0,start=begin+1,m=1; m<=nmod; m++,start+=lenM) {
	   K = M[m]->K; scrmtx = M[m]->score; 
	   for(scr=INT4_MIN,r=start; r<= end; r++) {	
		s=D[m-1][r-lenM]+score_smtrx((sq+r),K,scrmtx); 
		/*** compare with site at previous seq. position ***/
                F[m][r] = MAXIMUM(Int4,D[m][r-1],F[m][r-1]);
                if(s < F[m][r]) s = F[m][r];
		D[m][r] = s; 
		if(s > scr) score=scr=s; 
	   }
	   lenM = LenSMatrix(M[m])-overlap[m];
	   end += lenM;
	}
	/*** 4. free allocated memory ***/
	r=overlap[0];
	for(m=0; m<=nmod; m++) { 
	  D[m]+=begin; F[m]+=begin; free(D[m]);free(F[m]);
	}
	free(D); free(F); 
	return score;
}

static Int4 aln_seq_overlap_smatrix(Int4 len, unsigned char *seq,
	Int4 nmod, smx_typ *M, Int4 *pos, char *overlap,
	Int4 (*score_smtrx)(register unsigned char *, register Int4,
        register Int4 **))
/*******************************************************************
 Perform a dynamic programming alignment of sequence E and m colinear 
 scoring matrices M and return the optimal score.  

 from calling environment: len = LenSeq(E); seq = XSeqPtr(E);

********************************************************************/
{
	Int4	r,m,s,lenM,score,max_r,start,end,**D,**F,scr;
	BooLean	**T,t;
	Int4	**scrmtx,K;
	unsigned char *sq;
	double	p;
	Int4	begin=-overlap[0];

	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(T,nmod+3,BooLean*);
	MEW(F,nmod+3,Int4*);
	for(m=0; m<= nmod; m++) { 
	  MEW(D[m],len+overlap[0]+3,Int4); D[m]-=begin;
	  MEW(F[m],len+overlap[0]+3,Int4); F[m]-=begin;
	  NEW(T[m],len+overlap[0]+3,BooLean); T[m]-=begin;
	}
	/*** 2. Dynamic programming step. ***/
	for(r=begin; r<=len; r++) { D[0][r] = 0; }
	for(start=begin,m=1; m<=nmod; m++){ 
	   D[m][start]=INT4_MIN; F[m][start]=INT4_MIN; 
	   start+= (LenSMatrix(M[m])-overlap[m]);
	}
	end = len-start+1-overlap[0]; 
	sq = seq-1;
	for(lenM=0,start=begin+1,m=1; m<=nmod; m++,start+=lenM) {
	   K = M[m]->K; scrmtx = M[m]->score; 
	   for(scr=INT4_MIN,r=start; r<= end; r++) {	
		s=D[m-1][r-lenM]+score_smtrx((sq+r),K,scrmtx); t=TRUE;
		/*** compare with site at previous seq. position ***/
                F[m][r] = MAXIMUM(Int4,D[m][r-1],F[m][r-1]);
                if(s < F[m][r]) { s = F[m][r]; t=FALSE; }
		T[m][r] = t; D[m][r] = s; 
		if(s > scr){ score=scr=s; max_r=r; }
#if 0
fprintf(stderr,"%d: r=%d; s=%d\n", m,r,s);
#endif
	   }
	   lenM = LenSMatrix(M[m])-overlap[m]; end += lenM;
	}
    if(pos!=NULL){
	/*** 3. Trace back step. ***/
	for(m=nmod,r=max_r; m > 0 && r > begin; ){
	   if(T[m][r]){
		pos[m]=r;	/** record optimum alignment positions **/
		m--;
// fprintf(stderr,"%d: r=%d\n",m,r);
		if(m <= 0) break;
		r-=LenSMatrix(M[m])-overlap[m];
	   } else r--;
	}
    }
	/*** 4. free allocated memory ***/
	for(m=0; m<=nmod; m++) {
	  D[m]+=begin; F[m]+=begin; T[m]+=begin;
	  free(D[m]);free(F[m]); free(T[m]);
	}
	free(D); free(T); free(F); 
// fprintf(stderr,"score = %d\n",score);
	return score;
}

Int4	AlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
	Int4 *pos)
{ 
	if(pos==NULL) return fast_aln_seq_smatrix(len,seq,nmod,M,score_smatrix);
	else return aln_seq_smatrix(len, seq, nmod, M, pos, score_smatrix); 
}

Int4	FastAlnSeqOverlapSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	char *overlap, smx_typ *M, Int4 *X[2])
/*******************************************************************
 Performs a fast dynamic programming alignment for global score.
********************************************************************/
{
	Int4	r,m,s,lenM,start,end,j,**scrmtx;
	Int4	*Dm1,*Tmp,*D,Fm1,F,K;
	unsigned char *sq;
	Int4	begin=-overlap[0];

	D = X[0]-begin; Dm1 = X[1]-begin;
	for(r=begin; r<=len; r++) Dm1[r]=0;
        for(start=begin,m=1; m<=nmod; m++) start+=LenSMatrix(M[m])-overlap[m];
	end = len-start+1-overlap[0]; seq--;
	for(lenM=0,start=begin+1,m=1; m<=nmod; m++) {
	   D[start-1] = Fm1 = INT4_MIN;
	   scrmtx = M[m]->score; K = M[m]->K;
	   for(sq=seq+start,r=start; r<=end; r++,sq++) {
		s=Dm1[r-lenM];
		for(j=K; j; j--) s+= scrmtx[j][sq[j]]; 
		F = D[r-1];
		if(F < Fm1) F = Fm1; else Fm1=F;
                if(s < F) s = F;
		D[r] = s; 
	   }
	   lenM = LenSMatrix(M[m])-overlap[m];
	   start+=lenM; end += lenM;
	   Tmp=Dm1; Dm1=D; D=Tmp; 
	}
	return s;
}

Int4	AlnSeqOverlapSMatrix(Int4 len, unsigned char *seq, Int4 nmod,
	char *overlap, smx_typ *M, Int4 *pos)
{
	if(pos==NULL)
	     return fast_aln_seq_overlap_smatrix(len,seq,nmod,M,
			overlap,score_smatrix);
	else return aln_seq_overlap_smatrix(len,seq,nmod,M,pos,
			overlap,score_smatrix);
}

////////// OVERLAP with GAPS /////////////////

Int4	gap_func_aln_seq_overlap_smatrix(Int4 len,unsigned char *seq,Int4 nmod, 
	smx_typ *M,Int4 *pos,Int4 **gapscore, char *overlap,
	Int4 (*score_smtrx)(unsigned char *,Int4,smx_typ))
/*******************************************************************
 Perform a dynamic programming alignment of sequence E and m colinear 
 scoring matrices M using the gapscores provided and allowing 
 overlaps; returns the optimal score and motif positions in *pos.  
********************************************************************/
{
	Int4	r,m,s,lenM,score,max_r,start,end,**D,scr;
	Int4	g0,g,s2,r2,s1,g_opt,endX,end2,gs;
	short	**T,t;
	double	p;
	Int4    begin=-overlap[0];

	// 0. Make sure that sequence is Int4 enough to accommodate motifs
	for(s=0,m=1; m<=nmod; m++){
	   s+=LenSMatrix(M[m])-overlap[m];
	} if(s > len) return SHRT_MIN;
	
	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(T,nmod+3,short*);
	for(m=0; m<= nmod; m++) { 
		MEW(D[m],len+overlap[0]+3,Int4); D[m]-=begin;
		MEW(T[m],len+overlap[0]+3,short); T[m]-=begin;
	}
	/*** 2. Initialization step. ***/
	for(r=begin; r<=len; r++) { D[0][r] = 0; }
	for(start=begin,m=1; m<=nmod; m++){ 	// set seq ends to -inf 
	   D[m][start]=INT4_MIN; 
	   start+=LenSMatrix(M[m])-overlap[m];
	}
	end = len-start+1-overlap[0]; max_r=begin; score = SHRT_MIN;
	/*** 3. Dynamic programming step. ***/
	for(start=begin+1,lenM=-overlap[0],m=1; m<=nmod; m++,start+=lenM) {
	  for(scr=INT4_MIN,r=start; r<= end; r++) {	
	     /*** compare with site at all previous sites using gapscores ***/
	     s = score_smtrx(seq,r,M[m]);
	     /*** NEW: reverse loop so can exit early ***/
/************************************************************************
          start        r 
           |                              Note: lenM == core length(m-1) 
           |-----------:===(m)====-----
    -------|-:.........:.|--------  
    -------:-:===(m-1)===:  g = -overlap[m-1]...             r2 = r-lenM
  ===(m-1)===:<-- g -->:.:  g = ...(r-(start+overlap[m-1]))  (r2 = start-lenM)
 :<-lenM ->:                                                             
 ************************************************************************/
	     end2 = r-start-overlap[m-1]; 
	     for(s2=INT4_MIN,r2=r-lenM, g=-overlap[m-1]; g <= end2; g++,r2--){  
		if(g >= 0) g0=g; else g0=0; // set gaps less than 0 to zero for now
		if((gs=gapscore[m-1][g0]) == SHRT_MIN) break; // within gap limit? 
		s1=D[m-1][r2]+ s + gs;
		if(s2 < s1){ s2 = s1; g_opt = g; }
#if 0
fprintf(stderr,"%d: r=%d; r2=%d;s=%d;s1=%d;s2=%d;g=%d;gs=%d\n",
		m,r,r2,s,s1,s2,g,gs);
#endif
	     }
	     D[m][r] = s2; T[m][r] = g_opt;
	     if(m == nmod && s2 > scr){ score=scr=s2; max_r=r; }
	  }
	  lenM = LenSMatrix(M[m])-overlap[m]; end += lenM;
	}
   if(pos != NULL){
	/*** 3. Trace back step. ***/
	for(m=nmod,r=max_r; m > 0 && r > begin; ){
	   	g = T[m][r];
		pos[m]=r;	// record optimum alignment positions
#if 0
fprintf(stderr,"%d: g=%d; r=%d\n",m,g,r);
#endif
		m--; if(m <= 0) break;
		r = r - g - LenSMatrix(M[m]); // don't subtract overlap!!
	}
   }
	/*** 4. free allocated memory ***/
	for(m=0; m<=nmod; m++) {
		D[m]+=begin; free(D[m]);T[m]+=begin; free(T[m]);
	}
	free(D); free(T);
#if 0
fprintf(stderr,"score = %d\n",score);
#endif
	return score;
}

Int4	GapFuncAlnSeqOverlapSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos, Int4 **gapscore,char *overlap)
{ return gap_func_aln_seq_overlap_smatrix(len,seq,nmod,M,pos,gapscore,
	overlap, ScoreSMatrix); }

////////// END OVERLAP with GAPS /////////////////
////////// CtermFix OVERLAP with GAPS /////////////////

Int4    GapFuncAlnSeqOverlapSMatrixCterm(Int4 len, unsigned char *seq, Int4 nmod,
        smx_typ *M, Int4 *pos, Int4 **gapscore,char *overlap)
{ return gap_func_aln_seq_overlap_smatrix_Cterm(len,seq,nmod,M,pos,gapscore,
        overlap, ScoreSMatrix); }

Int4	gap_func_aln_seq_overlap_smatrix_Cterm(Int4 len,unsigned char *seq,
	Int4 nmod, smx_typ *M,Int4 *pos,Int4 **gapscore, char *overlap,
	Int4 (*score_smtrx)(unsigned char *,Int4,smx_typ))
/*******************************************************************
 same as gap_func_aln_seq_overlap_smatrix( ) but requires that
 the last motif is aligned against the Cterminus.
********************************************************************/
{
	Int4	r,m,s,lenM,score,max_r,start,end,**D,scr;
	Int4	g0,g,s2,r2,s1,g_opt,endX,end2,gs;
	short	**T,t;
	double	p;
	Int4    begin=-overlap[0];

	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(T,nmod+3,short*);
	for(m=0; m<= nmod; m++) { 
		MEW(D[m],len+overlap[0]+3,Int4); D[m]-=begin;
		MEW(T[m],len+overlap[0]+3,short); T[m]-=begin;
	}
	/*** 2. Initialization step. ***/
	for(r=begin; r<=len; r++) { D[0][r] = 0; }
	for(start=begin,m=1; m<=nmod; m++){ 	// set seq ends to -inf 
	   D[m][start]=INT4_MIN; 
	   start+=LenSMatrix(M[m])-overlap[m];
	}
	end = len-start+1-overlap[0]; max_r=begin; score = SHRT_MIN;
	/*** 3. Dynamic programming step. ***/
	for(start=begin+1,lenM=-overlap[0],m=1; m<=nmod; m++,start+=lenM) {
	  if(m == nmod){
	     endX = end - overlap[m-1];
	     for(scr=INT4_MIN,r=start; r< end; r++) {	
	       if(r==endX) s = score_smtrx(seq,r,M[m]);
	       else s = -9999;
	       end2 = r-start-overlap[m-1]; 
	       for(s2=INT4_MIN,r2=r-lenM, g=-overlap[m-1]; g <= end2; g++,r2--){  
		if(g >= 0) g0=g; else g0=0; // set gaps less than 0 to zero for now
		if((gs=gapscore[m-1][g0]) == SHRT_MIN) break; // within gap limit? 
		s1=D[m-1][r2]+ s + gs;
		if(s2 < s1){ s2 = s1; g_opt = g; }
	       }
	       D[m][r] = s2; T[m][r] = g_opt;
	       if(m == nmod && s2 > scr){ score=scr=s2; max_r=r; }
	     }
	  } else for(scr=INT4_MIN,r=start; r<= end; r++) {	
	     s = score_smtrx(seq,r,M[m]);
	     end2 = r-start-overlap[m-1]; 
	     for(s2=INT4_MIN,r2=r-lenM, g=-overlap[m-1]; g <= end2; g++,r2--){  
		if(g >= 0) g0=g; else g0=0; // set gaps less than 0 to zero for now
		if((gs=gapscore[m-1][g0]) == SHRT_MIN) break; // within gap limit? 
		s1=D[m-1][r2]+ s + gs;
		if(s2 < s1){ s2 = s1; g_opt = g; }
	     }
	     D[m][r] = s2; T[m][r] = g_opt;
	     if(m == nmod && s2 > scr){ score=scr=s2; max_r=r; }
	  }
	  lenM = LenSMatrix(M[m])-overlap[m]; end += lenM;
	}
   if(pos != NULL){
	/*** 3. Trace back step. ***/
	for(m=nmod,r=max_r; m > 0 && r > begin; ){
	   	g = T[m][r];
		pos[m]=r;	// record optimum alignment positions
		m--; if(m <= 0) break;
		r = r - g - LenSMatrix(M[m]); // don't subtract overlap!!
	}
   }
	/*** 4. free allocated memory ***/
	for(m=0; m<=nmod; m++) {
		D[m]+=begin; free(D[m]);T[m]+=begin; free(T[m]);
	}
	free(D); free(T);
	return score;
}

////////// END CtermFix OVERLAP with GAPS /////////////////

Int4	LocalAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
	Int4 *pos)
{
// char overlap[]={3,3,3,3,3,3,3,3,3,3,3,3,3,3};
	if(pos==NULL) return fast_aln_seq_smatrix(len,seq,nmod,M,local_score_smatrix);
	else return aln_seq_smatrix(len, seq, nmod, M, pos,local_score_smatrix); 
	// else return aln_seq_overlap_smatrix(len,seq,nmod,M,pos,overlap,local_score_smatrix);
}

Int4	fast_aln_seq_smatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
	Int4 (*score_smtrx)(register unsigned char *, register Int4, 
	register Int4 **))
/*******************************************************************
 Perform a dynamic programming alignment of sequence E and m colinear 
 scoring matrices M and return the optimal score.  
 from calling environment: len = LenSeq(E); seq = XSeqPtr(E);
 Note: this uses modified Gotoh method -> O(m*n).
 Note: this does not allow models to be deleted.
********************************************************************/
{
	Int4	r,m,s,lenM,score,start,end,**D,**F,adjust,scr;
	Int4	**scrmtx,K;
	unsigned char *sq;
	double	p;

	// 0. Make sure that sequence is Int4 enough to accommodate motifs
	for(s=0,score=0,m=1; m<=nmod; m++){
		s+=LenSMatrix(M[m]);
	} if(s > len) return SHRT_MIN;
	
	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(F,nmod+3,Int4*);
	for(m=0; m<= nmod; m++) { 
		MEW(D[m],len+3,Int4); MEW(F[m],len+3,Int4);
	}
	/*** 2. Dynamic programming step. ***/
	for(r=0; r<=len; r++) { D[0][r] = 0; }
	for(start=0,m=1; m<=nmod; m++){ 
	   D[m][start]=INT4_MIN; F[m][start]=INT4_MIN; 
	   start+=LenSMatrix(M[m]);
	}
	end = len-start+1;
	sq = seq-1;
	for(lenM=0,start=1,m=1; m<=nmod; m++,start+=lenM) {
	   K = M[m]->K; scrmtx = M[m]->score; 
	   for(scr=INT4_MIN,r=start; r<= end; r++) {	
		s=D[m-1][r-lenM]+score_smtrx((sq+r),K,scrmtx); 
		/*** compare with site at previous seq. position ***/
                F[m][r] = MAXIMUM(Int4,D[m][r-1],F[m][r-1]);
                if(s < F[m][r]) s = F[m][r];
		D[m][r] = s; 
		if(s > scr) score=scr=s; 
	   }
	   lenM = LenSMatrix(M[m]);
	   end += LenSMatrix(M[m]);
	}
	/*** 4. free allocated memory ***/
	for(m=0; m<=nmod; m++) {free(D[m]); free(F[m]);}
	free(D); free(F); 
	return score;
}


Int4	aln_seq_smatrix(Int4 len, unsigned char *seq, Int4 nmod, smx_typ *M,
	Int4 *pos, Int4 (*score_smtrx)(register unsigned char *, register Int4,
        register Int4 **))
/*******************************************************************
 Perform a dynamic programming alignment of sequence E and m colinear 
 scoring matrices M and return the optimal score.  

 from calling environment: len = LenSeq(E); seq = XSeqPtr(E);

 Find maximum score and trace back to get alignment.

        |  0  |  1  |  2  |  3  |  
      --+-----+-----+-----+-----+
      0 |  0  |-inf |  ?  |  ?  |  
      --+-----+-----+-----+-----+
      1 |  0  |s1,1 |  ?  |  ?  |  
      --+-----+---- +-----+-----+
        :     :     :     :     :
      --+-----+-----+-----+-----+
      m1|  0  |s1,m1|-inf |  ?  |  
      --+-----+-----+-----+-----+
      m1|  0  |s1,  |s2,  |  ?  |  
      +1|     |m1+1 |m1+1 |     |  
      --+-----+-----+-----+-----+
        :     :     :     :     :
      --+-----+-----+-----+-----+
      m1|  0  |s1,  |s2,  |-inf |  
     +m2|     |m1+m2|m1+m2|     |  
      --+-----+-----+-----+-----+
     m1+|  0  |s1,m1|s2,m1|s3,m1|  
    m2+1|     |+m2+1|+m2+1|+m2+1|  
      --+-----+-----+-----+-----+

 Note: this uses modified Gotoh method -> O(m*n).
 Note: this does not allow models to be deleted.
********************************************************************/
{
	Int4	r,m,s,lenM,score,max_r,start,end,**D,**F,adjust,scr;
	BooLean	**T,t;
	Int4	**scrmtx,K;
	unsigned char *sq;
	double	p;

	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(T,nmod+3,BooLean*); MEW(F,nmod+3,Int4*);
	for(m=0; m<= nmod; m++) { 
		NEW(T[m],len+3,BooLean); 
		// NEW(D[m],len+3,Int4); NEW(F[m],len+3,Int4);// DEBUG
		MEW(D[m],len+3,Int4); 
		MEW(F[m],len+3,Int4);
	}
	/*** 2. Dynamic programming step. ***/
	for(r=0; r<=len; r++) { D[0][r] = 0; }
	for(start=0,m=1; m<=nmod; m++){ 
	   D[m][start]=INT4_MIN; 
	   F[m][start]=INT4_MIN; 
	   start+=LenSMatrix(M[m]);
	}
	end = len-start+1;
	sq = seq-1;
	for(lenM=0,start=1,m=1; m<=nmod; m++,start+=lenM) {
	   K = M[m]->K; scrmtx = M[m]->score; 
	   for(scr=INT4_MIN,r=start; r<= end; r++) {	
		s=D[m-1][r-lenM]+score_smtrx((sq+r),K,scrmtx); t=TRUE;
		/*** compare with site at previous seq. position ***/
                F[m][r] = MAXIMUM(Int4,D[m][r-1],F[m][r-1]);
                if(s < F[m][r]) { s = F[m][r]; t=FALSE; }
		T[m][r] = t; D[m][r] = s; 
		if(s > scr){ score=scr=s; max_r=r; }
	   }
	   lenM = LenSMatrix(M[m]);
	   end += LenSMatrix(M[m]);
	}
#if 0
	printf("score = %d; max_r = %d\n",score,max_r);
	put_dp_matrix(D, T, len, nmod, M);
	put_dp_matrix(F, T, len, nmod, M);
#endif
    if(pos!=NULL){
	/*** 3. Trace back step. ***/
	for(adjust=-1,m=1; m<=nmod; m++) adjust+=LenSMatrix(M[m]);
	for(m=nmod,r=max_r; m > 0 && r > 0; ){
	   if(T[m][r]){
		pos[m]=r;	/** record optimum alignment positions **/
#if 0
		s = score_smtrx((sq+r),M[m]->K,M[m]->score); 
		p = SMatrixProb(s, M[m]);
		p = -log10(p*(double)(len-adjust));
		fprintf(stderr,"model %d: site %d-%d; score %d; pval=%.1f\n",
		   m,r,r+LenSMatrix(M[m])-1,s,p);
#endif 	/*** DEBUG ****/
		m--;
		if(m <= 0) break;
		r-=LenSMatrix(M[m]);
	   } else r--;
	}
    }
	/*** 4. free allocated memory ***/
	for(m=0; m<=nmod; m++) {free(D[m]);free(T[m]);free(F[m]);}
	free(D); free(T); free(F); 
	return score;
}

/************************ Spouge modifications ***************************/

Int4	LocalGapFuncAlnSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos, Int4 **gapscore)
{ return gap_func_aln_seq_smatrix(len,seq,nmod,M,pos,gapscore,LocalScoreSMatrix); }

Int4	GapFuncAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos, Int4 **gapscore)
{ return gap_func_aln_seq_smatrix(len, seq, nmod, M, pos, gapscore, ScoreSMatrix); }

Int4	gap_func_aln_seq_smatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos, Int4 **gapscore, 
	Int4 (*score_smtrx)(unsigned char *, Int4 , smx_typ))
/*******************************************************************
 Perform a dynamic programming alignment of sequence E and m colinear 
 scoring matrices M using the gapscores provided; returns the optimal 
 score and motif positions in *pos.  

 from calling environment: len = LenSeq(E); seq = XSeqPtr(E);

 Find maximum score and trace back to get alignment.

        |  0  |  1  |  2  |  3  |  
      --+-----+-----+-----+-----+
      0 |  0  |-inf |  ?  |  ?  |  
      --+-----+-----+-----+-----+
      1 |  0  |s1,1 |  ?  |  ?  |  
      --+-----+---- +-----+-----+
        :     :     :     :     :
      --+-----+-----+-----+-----+
      m1|  0  |s1,m1|-inf |  ?  |  
      --+-----+-----+-----+-----+
        :     :     :     :     :
      --+-----+-----+-----+-----+
      m1|  0  |s1,  |s2,  |-inf |  
     +m2|     |m1+m2|m1+m2|     |  
      --+-----+-----+-----+-----+
     m1+|  0  |s1,m1|s2,m1|s3,m1|  
    m2+1|     |+m2+1|+m2+1|+m2+1|  
      --+-----+-----+-----+-----+

 Note: this uses modified Gotoh method -> O(m*n).
 Note: this does not allow models to be deleted.

  T[m][r] = traceback matrix giving optimum gap back to previous motif.   
  D[m][r] = best previous score up to starting point r for mth column.

********************************************************************/
{
	Int4	r,m,s,lenM,score,max_r,start,end,**D,adjust,scr;
	Int4	g,s2,r2,s1,g_opt,end2,gs;
	unsigned short	**T,t;
	double	p;

	BooLean	debug=FALSE;

	/*** 1. Allocate memory. ***/
	MEW(D,nmod+3,Int4*); MEW(T,nmod+3,unsigned short*);
	for(m=0; m<= nmod; m++) { 
		NEW(T[m],len+3,unsigned short); 
		if(debug) NEW(D[m],len+3,Int4);	/**DEBUG**/
		else MEW(D[m],len+3,Int4);
	}
	/*** 2. Initialization step. ***/
	for(r=0; r<=len; r++) { D[0][r] = 0; }
	for(start=0,m=1; m<=nmod; m++){ 	/*** set seq ends to -inf ***/
	   D[m][start]=INT4_MIN; 
	   start+=LenSMatrix(M[m]);
	}
	end = len-start+1; start=1; lenM=0;
	/*** 3. Dynamic programming step. ***/
	for(m=1; m<=nmod; m++,start+=lenM) {
	  for(scr=INT4_MIN,r=start; r<= end; r++) {	
	     /*** compare with site at all previous sites using gapscores ***/
	     s = score_smtrx(seq,r,M[m]);
	     /*** NEW: reverse loop so can exit early ***/
	     end2 = r-start;
	     for(s2=INT4_MIN,r2=start-lenM+end2, g=0; g <= end2; g++,r2--){  
		if((gs = gapscore[m-1][g]) == SHRT_MIN) break; /** w/in gap limits? **/
		s1=D[m-1][r2]+ s + gs;
		if(s2 < s1){ s2 = s1; g_opt = g; }
	     } /**** NEW ****/
#if 0
	     /** OLD ***
	     for(s2=INT4_MIN,r2=start-lenM, g=r-start; g >= 0; g--,r2++){
		s1=D[m-1][r2]+ s + gapscore[m-1][g]; 
		if(s2 < s1){ s2 = s1; g_opt = g; }
	     } /**** OLD ****/
#endif
	     D[m][r] = s2; T[m][r] = g_opt;
	     if(m == nmod && s2 > scr){ score=scr=s2; max_r=r; }
	  }
	  lenM = LenSMatrix(M[m]); end += lenM;
	}
	if(debug){
		printf("\n\nscore = %d; max_r = %d\n",score,max_r);
	} /**** DEBUG ****/
   if(pos != NULL){
	/*** 3. Trace back step. ***/
	for(adjust=-1,m=1; m<=nmod; m++) adjust+=LenSMatrix(M[m]);
	for(m=nmod,r=max_r; m > 0 && r > 0; ){
	   	g = T[m][r];
		pos[m]=r;	/** record optimum alignment positions **/
		/**********/
	   if(debug){
		s = score_smtrx(seq,r,M[m]);
		p = SMatrixProb(s, M[m]);
		p = -log10(p*(double)(len-adjust));
		printf("(gapscore[%d]=%d) model %d: site %d-%d; score %d; pval=%.1f\n",
		   g,gapscore[m-1][g],m,r,r+LenSMatrix(M[m])-1,s,p);
	   }	/*** DEBUG ****/
		m--;
		if(m <= 0) break;
		r = r - g - LenSMatrix(M[m]);
	}
   }
	/*** 4. free allocated memory ***/
	for(m=0; m<=nmod; m++) {free(D[m]);free(T[m]);}
	free(D); free(T);
	return score;
}

Int4	**BlastCircularSMatrix(Int4 maxlen, smx_typ M)
/*********************************************************************
 Create a matrix bounded by SHRT_MIN values on ends; for use by
     gapxdrop.c alignment functions.  
/*********************************************************************/
{
	Int4	lenM,i,j,k,K,r;
	a_type	A = SMatrixA(M);
	unsigned char	*qseq,*maxseq;

	lenM = LenSMatrix(M);
	if(M->posMatrix != NULL){
		i = M->posMtrxLen;
		free(M->posMatrix[0]); free(M->posMatrix[i+1]);
		free(M->posMatrix);
		if(M->qseq != NULL) free(M->qseq);
	}
	M->posMtrxLen = K = 2*maxlen + lenM;
	i = (maxlen + lenM) % lenM;
	M->posMtrxOffset = maxlen - i;
printf("offset = %d\n\n",M->posMtrxOffset);
	NEWP(M->posMatrix, K +3, Int4);
	NEW(M->posMatrix[0], M->nlet+3, Int4);
        for(r=0; r<=(M->nlet +1); r++) M->posMatrix[0][r]=SHRT_MIN;
	NEW(M->posMatrix[K+1], M->nlet+3, Int4);
        for(r=0; r<=(M->nlet +1); r++) M->posMatrix[K+1][r]=SHRT_MIN;
	
	NEW(maxseq, lenM +2, unsigned char);
	MaxSegSMatrix(maxseq, M);
	NEW(M->qseq, K +3, unsigned char);
	qseq =  M->qseq; qseq[0] = M->nlet +1;
        for(k=0,i=1; i<=K;i++) {
		if(k == lenM) k = 1; else k++;
		qseq[i] = maxseq[k];
		M->posMatrix[i] = M->score[k];
        }
	qseq[i] = M->nlet +1;
	free(maxseq);
	return M->posMatrix;
}

unsigned char	*CircularSMatrixQuery(smx_typ M)
{
	if((M->posMatrix == NULL) || M->qseq == NULL) 
		print_error("input error in CircularSMatrixQuery( )");
	return M->qseq;
}

Int4	SMatrixQueryLen(smx_typ M) { return M->posMtrxLen; }

Int4	CircularSMatrixPos(Int4 seqlen, Int4 mtxpos, smx_typ M)
/** returns the position in the "circular" matrix to begin the extension **/
{
	if((M->posMatrix == NULL) || M->qseq == NULL ||
		(2*seqlen + LenSMatrix(M)) > M->posMtrxLen) {
		print_error("input error in CircularSMatrixPos( )");
	}
	return (M->posMtrxOffset + mtxpos);
}

double	LocalHistSeqSMatrix(FILE *fp, Int4 len, unsigned char *seq, smx_typ M)
{
	Int4	i,s,max=-9999;
	double	ave,sd,score;
	h_type  H;

	H = Histogram("smatrix local scores",-100,200,2.0);
	for(i=-5; i<=len-M->K+5; i++){
		s = LocalScoreSMatrix(seq, i, M);
		if(s > max) max = s;
		IncdHist(s, H);
	}
	ave = MeanHist(H); sd = sqrt(VarianceHist(H));
	score = ((double)max-ave)/sd;
	if(fp != NULL) PutHist(fp,60,H); 
	NilHist(H);
	if(fp != NULL) fprintf(fp,"maximum (%d) is %.1f sd above mean\n\n",max,score);
	return score;
}

double	HistSeqSMatrix(FILE *fp, Int4 len, unsigned char *seq, smx_typ M)
{
#if 1
	Int4	i,s;
	double	ave,sd,score,prob,max=-999999;
	h_type  H;

	H = Histogram("smatrix probabilities",-100,200,0.25);
	for(i=-5; i<=len-M->K+5; i++){
		s = ScoreSMatrix(seq, i, M);
		prob = -log10((double)(len - M->K+10)*SMatrixProb(s, M));
		if(prob > 1.0) {
		   fprintf(stderr,"length = %d; i = %d; prob = %.3f\n",
			len,i,prob);
		}
		if(prob > max) max = prob;
		IncdHist(prob, H);
	}
	fprintf(stderr,"\n");
	ave = MeanHist(H); sd = sqrt(VarianceHist(H));
	score = ((double)max-ave)/sd;
	if(fp != NULL) PutHist(fp,60,H); NilHist(H);
	if(fp != NULL) fprintf(fp,"maximum (%.2f) is %.1f sd above mean\n\n",max,score);
	return score;
#endif
#if 0
	Int4	i,s,max=-9999;
	double	ave,sd,score;
	h_type  H;

	// H = Histogram("smatrix scores",-5,len+5,1.0);
	H = Histogram("smatrix scores",-100,200,2.0);
	for(i=-5; i<=len-M->K+5; i++){
		s = ScoreSMatrix(seq, i, M);
		if(s > max) max = s;
		// IncdMHist(i, s, H);
		IncdHist(s, H);
	}
	ave = MeanHist(H); sd = sqrt(VarianceHist(H));
	score = ((double)max-ave)/sd;
	if(fp != NULL) PutHist(fp,60,H); NilHist(H);
	if(fp != NULL) fprintf(fp,"maximum (%d) is %.1f sd above mean\n\n",max,score);
	return score;
#endif
}

Int4	BestAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos)
/*******************************************************************
 Finds the best match for all nmod smatrices...
********************************************************************/
{
	Int4		r,m,s,lenM,end,j,**scrmtx,K,score,s_max;
	unsigned char	*sq;

	for(score=0,m=1; m<=nmod; m++) {
	   lenM = LenSMatrix(M[m]);
	   end = len - lenM + 3;
	   scrmtx = M[m]->score; K = M[m]->K;
	   for(s_max=-9999,sq=seq-4,r=-3; r<=end; r++,sq++) {
		for(s=0,j=K; j; j--) s+= scrmtx[j][sq[j]]; 
		if(s > s_max) { s_max = s; pos[m]=r; }
	   }
	   score += s_max;
	}
	return score;
}

Int4	BestLocalAlnSeqSMatrix(Int4 len, unsigned char *seq, Int4 nmod, 
	smx_typ *M, Int4 *pos)
/*******************************************************************
 Finds the best match for all nmod smatrices...
********************************************************************/
{
	Int4		r,m,s,lenM,end,j,**scrmtx,K,score,s_max,max;
	unsigned char	*sq;

	for(score=0,m=1; m<=nmod; m++) {
	   lenM = LenSMatrix(M[m]);
	   end = len - lenM + 3;
	   scrmtx = M[m]->score; K = M[m]->K;
	   for(s_max=-9999,sq=seq-4,r=-3; r<=end; r++,sq++) {
                for(s=max=0,j=K; j; j--) {
                    s += scrmtx[j][sq[j]];
                    if(s < 0) s=0; else if(s > max) max = s;
                } s = max;
		if(s > s_max) { s_max = s; pos[m]=r; }
	   }
	   score += s_max;
	}
	return score;
}

Int4	GapXDropSMatrix(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	smx_typ M, Int4 *mpos, Int4 maxgap)
////////////////////////////////////////////////////////////////////////
//  Similar to AlnSeqSMatrixSW( )
//
//  Example for maxgap = 2  			 D = DEL; M = MATCH; I = INS;
//       D   R   F   R   G   L   V   L   I   S   (concensus)
//  ---+---+---+---+---+---+---+---+.                 +---+
//   , - ,,- ,,| ,,| ,,| ,,| ,,| ,,|             0    |DMI|
//  ---+---+---\---+---+---+---+---+..                +---+
// I,, | 0 | 1 |[2]|** |   |   |   |             1   
//  ---+---+---+---\---+---+---+---+...                ' ' = undefined.
// N,, | 1 | 0 | 1 |[2]|** |   |   |             2                  
//  ---+---+---+---+---\---+---+---+...                '*' = NEG_INF.
// S,, | 2 | 1 | 0 | 1 |[2]|** |   |   :         :   
//  ---+---+---+---+---+---\---+---+---+---+---+-      ',' = zero.
// A,, | **| 2 | 1 | 0 | 1 |[2]|** |   |   |   | j-2  
//  ---+---+---+---+---+---+---\---+---+---+---+-     S
// V,, |   | **| 2 | 1 | 0 | 1 |[2]|** |   |   | j-1  E
//  ---+---+---+---+---+---+---+---\---+---+---+-     Q
// L,, |   |   | **| 2 | 1 |_0_| 1 |[2]|** |   | j    U
//  ---+---+---+---+---+---+---+---+---\---+---+-     E
// I   :   :   :   | **| 2 | 1 | 0 | 1 |[2]|** | j+1  N
//              ...+---+---+---+---+---+---\---+-     C
// S               |   | **| 2 | 1 | 0 | 1 |[2]| j+2  E
//               ..+---+---+---+---+---+---+---+-    
// P               |   |   | **| 2 | 1 | 0 | 1 | j+3  
//                .+---+---+---+---+---+---+---+-    
// N               |   |   |   | **| 2 | 1 | 0 | j+4 (n2)
//                 +---+---+---+---+---+---+---+-    
//   0   1   2  ... i-2 i-1  i  i+1 i+2 i+3 i+4
//                     MODEL                (n1)
//
//        DRFRGLVLIS model concensus
//          . ..::::
//        --INSAVLIS sequence
//
/////////////////////////////////////////////////////////////////////////
{
	Int4	m,i,j,k,r1,s,t,v,n1,score,max_i,max_j;
	Int4	*pos[3],p,w,**MAT,**T,**DEL,**INS;
	unsigned char	*seq1;
	char	*out[3];
	a_type	A = SMatrixA(M);
	const Int4	NEG_INF = -999999;
	Int4	end2,igaps,jgaps;

	/** get concensus sequence for smatrix **/
	n1 = LenSMatrix(M);
	NEW(seq1, n1+3, unsigned char); 
	MaxSegSMatrix(seq1, M); s = LenSMatrix(M); 

	/*** 1. Allocate memory. ***/
	MEW(MAT,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		MEW(MAT[i],n2+3,Int4); MEW(T[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
	}
	/*** 2a. Dynamic programming step. ***/
	// make GLOBAL ALIGNMENT with respect to profile (len profile = n1)
	MAT[0][0]=0;
  	for(i=1; i<=n1; i++){ 
		MAT[i][0] = INS[i][0] = 0; T[i][0] = -1;
	}
  	for(j=1; j<=n2; j++) { // no gap opening penalty for deletions from ends.
		DEL[0][j] = MAT[0][j] = 0;
	}
	for(w=a+b, score=0, i=1; i<= n1; ) {
	   j = MAXIMUM(Int4,i-maxgap,1);
	   end2 = MINIMUM(Int4,i+maxgap,n2);
	   if(j > 1) INS[i][j-1] = MAT[i][j-1] = NEG_INF;   // sets INS[i][0] = -inf;

	   for(  ; j<= end2; j++) {
	 	s = MAT[i-1][j-1] + ValSMatrix(i,seq2[j],M); t=0;  // match

                DEL[i][j] = MAXIMUM(Int4,MAT[i-1][j]-w,DEL[i-1][j]-b); // deletion
                if(s < DEL[i][j]) { s = DEL[i][j]; t=-1; }

                INS[i][j]=MAXIMUM(Int4,MAT[i][j-1]-w,INS[i][j-1]-b); // insertion 
                if(s < INS[i][j]) { s = INS[i][j]; t=1; }

		T[i][j] = t; MAT[i][j] = s; 
	   }
	   if(j <= n2) DEL[i][j]=MAT[i][j]=NEG_INF;
	   i++;
	}
	// 2b. Find optimum global score with respect to profile 
	max_i = n1; /*** global with respect to profile ***/
	j = MAXIMUM(Int4,max_i-maxgap,1);
	end2 = MINIMUM(Int4,max_i+maxgap,n2);
	for(score=INT4_MIN; j<= end2; j++) { // requires global alignment to model
		s = MAT[max_i][j]; 
		if(s > score){ score = s; max_j = j; }
	}
	/*** 3. Trace back step. ***/
	MEW(out[0],n1+n2+3,char); MEW(out[1],n1+n2+3,char);
	MEW(out[2],n1+n2+3,char);
	NEW(pos[1],n1+n2+3,Int4); NEW(pos[2],n1+n2+3,Int4);
	igaps=jgaps=0;
	for(p=1,i=max_i,j=max_j; i > 0 && j >= 0; ){
		switch(T[i][j]){
		  case -1:		// deletion in seq.
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			pos[1][p] = i; p++; i--;
			jgaps++;
			break;
		  case 0:		// match state
			v = ValSMatrix(i,seq2[j],M);
			pos[1][p] = i; pos[2][p] = j;
			out[1][p] = AlphaChar(seq1[i],A);
			out[2][p] = AlphaChar(seq2[j],A);
			if(seq1[i]==seq2[j]) out[0][p] = ':';
			else if(v >  0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; i--; j--;
			break;
		  case 1:		// insertion in seq.
			out[2][p] = AlphaChar(seq2[j],A);
			out[0][p] = ' '; out[1][p] = '-';
			pos[2][p] = j; p++; j--;
			igaps++;
			break;
		  default: 
		  fprintf(stderr,"i=%d; j=%d; t=%d; max_i=%d; max_j=%d; n1=%d; n2=%d; score=%d\n",
			i,j,T[i][j],max_i,max_j,n1,n2,score);
		  print_error("this should not happen"); break;
		}
	}
	i++;j++;
	*mpos = (j-1) - (i-1);
#if 1
printf("max_i = %d; max_j= %d; lenM = %d; lenS=%d; i=%d; j=%d; igaps=%d; jgaps=%d\n",
 	max_i,max_j,LenSMatrix(M),n2,i,j,igaps,jgaps);
	/** 4. Print out alignment **/
	printf("\n\n");
	for(i=p-1; i > 0 ; ){
	   v = MINIMUM(Int4,50,i); 
	   printf("     ");
	   for(j=i; j > i-v ; j--) printf("%c",out[1][j]); 
	   printf(" concensus\n     ");
	   for(j=i; j > i-v ; j--) printf("%c",out[0][j]); printf("\n     ");
	   for(j=i; j > i-v ; j--) printf("%c",out[2][j]);
	   printf(" sequence\n  ");
	   for(j=i; j > i-v ; j--) {
		if(pos[2][j] && pos[2][j] % 10 == 0) {
			printf("%4d",pos[2][j]); j-=3;
		} else printf(" ");
	   }
	   printf("\n\n"); i-=v;
	}
printf("score = %d\n",score);
#endif
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(T); free(DEL); free(INS); free(seq1); 
	free(out[0]); free(out[1]); free(out[2]);
	free(pos[1]); free(pos[2]);
	return score;
}

///////////////// SW-ALIGNEMNT WITH GAP FUNCTION //////////////////
//
// Perform a Smith-Waterman alignment on multiple smatrix .M sequence seq2
// returns optimum subalignment score.  W(k) = a + bk.
// D[0][j]=D[i][0]=0;
//              MAX{ D[i-1][j-1] + S(r1,r2),
// D[i][j] =            D[i-k][j] - W(k) for k=1..i-1,
//                      D[i][j-k]- W(k) for k=1..j-1}.
//
// a = gap open; b = gap extend;
// Find maximum score and trace back to get alignment.
// see T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197.
//
//   Uses O(m*n) Gotoh method for penalties where W(k)=u*k+v where u,v >= 0.
//             .....
//             :DMI:                D = DEL; M = MATCH; I = INS;
//             .....    gap = TB[m][j]
//      .......m....... | ..........m+1.......... 
//     : D   R   F   R :v: G   L   V   L   I   S : (concensus)
//  ....................................                
//   , : ,,: ,,: ,,: ,,: : ,,: ,,: ,,:             0          
//  ...:...:...:...:...:.:...:...:...:...                'a' = affine penalty
// I,, :   :   :   :   : :   :   : I : I :         1
//  ...:...:...:...:...:.:...:...:.|.\.|.:               ' ' = undefined.
// N,, :   :   :   :   : :   :   : H :\H :         2
//  ...:...:...:...0...:.:...:...:.|.:.|.:               '*' = NEG_INF.
// S,, :   :   : ->1`?--+->? :   : a : a :         :
//  ...:...:...:...:...:|:...\...:.|.:.|.:.........      ',' = zero.
// A,, :   :   :   :   :|:   : ? : d_: d :   :   : j-2
//  ...:...:...:...:...:|:...:...:...\.|.:...:...:.     S
// V,, :   :   :   :   :|:   :   :   :\H :   :   : j-1  E
//  ...:...:...:...:...:|:...:...:...:.|.:...:...:.     Q
// L,, :   :   :   :   :|:   :   :   : e :   :   : j    U
//  ...:...:...:...:...:|....:...:...:.|.:...:...:.     E
// I   :   :   :   :   :+->? .   :   : k_:   :   : j+1  N
//             :...:...:.:...\...:...:...0...:...:.     C
// S               :   : :   : ? :   :   :\G_:   : j+2  E
//               ..:...:.:...:...:...:...:...\...:.
// P               :   : :   :   :   :   :   :   : j+3
//                .:...:.:...:...:...:...:...:...:.
// N               :   : :   :   :   :   :   :   : j+4 (n2)
//                 :...:.:...:...:...:...:...:...:.
//   0   1   2  ... i-2  i-1  i  i+1 i+2 i+3 i+4
//   k = 1   2  ...lenM   1   2   ...       lenM
//                         MODEL            (n1)
//
//        AAAA   BBBBBB
//        DRFR...GLVLIS model concensus
//          .    ..::::
//       XXINLAESAVLISKHI sequence
//
//                  :
//    INS: [////]---|   Don't use this for k == lenM!
//         [\\\\\\\]|
//
//                  :
//    DEL: [///////]|   Don't use this for k == 1
//         [\\\\]---|
//
//                  :
//    MAT: [///////]|
//         [\\\\\\\]|
//
///////////////////////////////////////////////////////////////////
char	*gapped_aln_trace_smatrixSW(Int4 a,Int4 b,Int4 n2,unsigned char *seq2, 
	Int4 nmod,smx_typ *M,Int4 **gapscore,Int4 *J,Int4 *alnscore,Int4 *start)
{
	Int4	m,i,j,k,r1,s,t,v,n1,score,max_i,max_j,mscore;
	Int4	**MAT,**T,**DEL,**INS,**TB;
	Int4	j2,s1,s0,g,gs,g_opt,w=a+b,J0;
	unsigned char	*seq1;
	short	*mtf;
	a_type	A = SMatrixA(M[1]);
	Int4	gINS,*gDEL,t0;
	char	*operation;

	if(gapscore != NULL) gapscore--; // want gapscore[m-1] for block m.
	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	NEW(mtf, n1+3, short); NEWP(TB, nmod+3, Int4); 
	for(s=0, m=1; m <= nmod; m++){
	    NEW(TB[m], n2+3, Int4); 
	    // MaxSegSMatrix(seq1+s, M[m]); // NEED ONLY FOR OUTPUT.
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; mtf[s] = m; }
	}
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
#if 1 
		MEW(MAT[i],n2+3,Int4); MEW(T[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
#endif
#if 0 // debug using matrix output.
		NEW(MAT[i],n2+3,Int4); NEW(T[i],n2+3,Int4); 
		NEW(DEL[i],n2+3,Int4); NEW(INS[i],n2+3,Int4); 
#endif
	}
	// Make GLOBAL ALIGNMENT with respect to profile,
	MAT[0][0] = 0; DEL[0][0] = 0; INS[0][0] = 0;
	for(i=1; i<= n1; i++) {  
		DEL[i][0] = DEL[i-1][0] - b;
		INS[i][0] = INS[i-1][0] - b;
		MAT[i][0] = INS[i-1][0];  // Local on ends of profile 
		T[i][0] = 1; // for full alignment
	}
	for(j=1; j<= n2; j++) { // Make LOCAL with respect to sequence.
		MAT[0][j] = 0; T[0][j]=-1; // traceback full alignment 
		DEL[0][j] = INS[0][j] = SHRT_MIN;
	}
	// 2. Dynamic programming step. 
	NEW(gDEL,n2+3,Int4);
	for(k=m=i=1; i<= n1; i++) {
	   gINS=0;
	   if(k==LenSMatrix(M[m])) {
	     for(j=1; j<= n2; j++) {  // Eliminate insert state...
                t=0;
	 	s = MAT[i-1][j-1] + ValSMatrix(k,seq2[j],M[m]);
		if(gapscore == NULL){ // no affine penalty for gapscore.
                   if((s0=MAT[i][j-1]) > (s1=INS[i][j-1])){
                        gINS=-1;  INS[i][j] = s0;
                   } else { gINS--; INS[i][j] = s1; }
                   if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }
		}

                if((s0=MAT[i-1][j]-w) > (s1=DEL[i-1][j]-b)){
                        gDEL[j]=1;  DEL[i][j] = s0;
                } else { gDEL[j]++;  DEL[i][j] = s1; }
                if(s < DEL[i][j]){ s=DEL[i][j]; t=gDEL[j]; }

		T[i][j] = t; MAT[i][j] = s;
	     }
	   } else if(k==1 && m > 1 && gapscore != NULL) {
	     for(j=1; j<= n2; j++) {  // Use gap function here.
		// Don't jump over gap function regions! Can't use t=gDEL[j].
                DEL[i][j] = MAXIMUM(Int4,MAT[i-1][j]-w,DEL[i-1][j]-b);
                s = DEL[i][j]; t=1;

                if((s0=MAT[i][j-1]-w) > (s1=INS[i][j-1]-b)){
                        gINS=-1;  INS[i][j] = s0;
                } else { gINS--; INS[i][j] = s1; }
                if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }

		v = ValSMatrix(k,seq2[j],M[m]);
		for(s1=SHRT_MIN,j2=j-1, g=0; j2 > 0; g++,j2--){  
		    if((gs = gapscore[m][g]) == SHRT_MIN) break; // over max gap?
		    s0 = MAT[i-1][j2] + v + gs;
		    if(s1 < s0){ s1 = s0; g_opt = g; }
		}
                if(s < s1) { s = s1; t=0; TB[m][j] = -g_opt; }
		T[i][j]=t; MAT[i][j] = s;
	     }
	   } else { 		      // Use affine gap penalty
	      for(j=1; j<= n2; j++) {
                t=0; s = MAT[i-1][j-1] + ValSMatrix(k,seq2[j],M[m]);

                if((s0=MAT[i-1][j]-w) > (s1=DEL[i-1][j]-b)){
                        gDEL[j]=1;  DEL[i][j] = s0;
                } else { gDEL[j]++;  DEL[i][j] = s1; }
                if(s < DEL[i][j]){ s=DEL[i][j]; t=gDEL[j]; }

                if((s0=MAT[i][j-1]-w) > (s1=INS[i][j-1]-b)){
                        gINS=-1;  INS[i][j] = s0;
                } else { gINS--; INS[i][j] = s1; }
                if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }

                T[i][j] = t; MAT[i][j] = s;
	      }
	   }
	   k++;
	   if(k > LenSMatrix(M[m])) { k = 1; m++; }
	}
#if 0
Int4 startS=135,endS=155,motif=4;
	for(i=1; i<=n1; i++) if(mtf[i]==motif) break; 
	fprintf(stderr,"DEL[i][j]:\n");
	put_dp_matrixSW(DEL, T, n2, M[motif],seq1,seq2,A,i,startS,endS);
	fprintf(stderr,"INS[i][j]:\n");
	put_dp_matrixSW(INS, T, n2, M[motif],seq1,seq2,A,i,startS,endS);
	fprintf(stderr,"MAT[i][j]:\n");
	put_dp_matrixSW(MAT, T, n2, M[motif],seq1,seq2,A,i,startS,endS);
	fprintf(stderr,"SMX[i][j]:\n");
	put_dp_matrixSMX(n2, M[motif],seq1,seq2,A,i,startS,endS);
	fprintf(stderr,"T[i][j]:\n");
	put_dp_matrixSW(T, T, n2, M[motif],seq1,seq2,A,i,startS,endS);
#endif
	// 2b. Find optimum global score with respect to profile.
	max_i=n1;		// global with respect to profile.
	for(score=INT4_MIN, j=1; j<= n2; j++)
		if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
	/*** 3. Trace back step. ***/
	// operations: 
	// 'i' = insertion in sequence outside of motifs
	// 'I' = insert in sequence within motif (need to delete this)
	//  'M' = match to start of a motif block
	//  'm' = match to other sites in motif
	//  'D' = deletion of sequence within motif
	//  'd' = deletion of sequence outside of motif (not used)
	NEW(operation,n1+n2+3,char); operation[0]='E';
	m=nmod; k=LenSMatrix(M[m]);
	for(J0=0,j=n2; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }
	// for(i=n1; i > 0 || j > 0; ){  // full alignment...
	for(i=n1; i > 0; ){  // full alignment...Stop at end of profile
	  t0 = T[i][j];
	  do {
	    if(t0 > 0){ t=1; t0--; }
            else if(t0 < 0){ t=-1; t0++; } else { t=0; }
	    switch(t){
		case 0:
		   i--; J0++;
		   if(k==1) operation[J0] = 'M'; else operation[J0] = 'm';
		   if(gapscore != NULL && k==1 && m > 1){
			for(g=TB[m][j], j--; g < 0; g++){
			    J0++; operation[J0] = 'i'; j--;
			}
		   } else j--; k--; break;
		case 1:  // Gap ('-') is in sequence; add X's for gaps
			J0++; 
			if(k==1) operation[J0] = 'D'; else operation[J0] = 'd';
			i--;  k--; break;
		case -1:	// Gap ('-') is in profile
			if(m < 1 || k==LenSMatrix(M[m]))
				{ J0++; operation[J0] = 'i'; }
			else { J0++; operation[J0] = 'I'; }
			j--; break;
		default: print_error("this should not happen"); break;
	     }
	     if(k==0){ m--; if(m > 0) k=LenSMatrix(M[m]); }
	  } while(t0 != 0);
	}
#if 1	// Temporary FIX for "mIDE" at end (beginning once reversed).
	if(operation[J0]== 'D' && operation[J0-1] == 'I'){
	   operation[J0-1]='M'; operation[J0]='E';
	   *start=j+1; 
	} else {
	   *start=j+1; J0++; operation[J0]='E';
	}
#else	// original code that allows "mIDE" at ends...
	*start=j+1; J0++; operation[J0]='E';
#endif
#if 0	// reverse in calling function.
	while(j > 0) { J0++; operation[J0] = 'i'; j--; }
	/*** 4b. reverse operations ***/
        char *operation2;
        NEW(operation2,J0+3,char); operation2[0] = 'E';
        for(j=J0,i=1; i <= J0; j--,i++){
                operation2[i] = operation[j];
        } free(operation); operation2[i] = 'E'; i++; operation2[i] = 0;
#endif
	/*** 5. free allocated memory ***/
	for(m=1; m <= nmod; m++) free(TB[m]); free(TB);
	for(i=0; i<=n1; i++) {free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(T); free(DEL); free(INS); 
	free(seq1); free(mtf); free(gDEL);
	*alnscore = score; *J = J0;
// printf("%s\n",operation2);
// printf("gapscore = %d\n",(Int4)gapscore);
        // return operation2;
        return operation;
}

char	*gapped_aln_seq_smatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *J, Int4 *alnscore)
{ 
  Int4	i,j,start,tracelen;
  char	*operation,*tmp;

  tmp=gapped_aln_trace_smatrixSW(a,b,n2,seq2,nmod,M,gapscore,J,alnscore,&start);
// std::cerr << "tmp traceback"; std::cerr << std::endl; std::cerr << tmp; std::cerr << std::endl;
  // Traceback returned is reversed and ends at start; needs to be reversed.
  tracelen=*J; 
  NEW(operation,tracelen+start+3,char); operation[0] = 'E';
  for(j=1;j < start; ) { operation[j] = 'i'; j++; }
  for(i=tracelen-1; i > 0; i--){ operation[j] = tmp[i]; j++; }
  free(tmp); operation[j] = 'E'; j++; operation[j] = 0;
  *J = j; // length of new trace;
// std::cerr << "real traceback"; std::cerr << std::endl; std::cerr << operation; std::cerr << std::endl;
  return operation;
}

char	*GapAlnTraceSMatrix(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *start)
{ Int4	s; return GapAlnTraceSMatrix2(a,b,len,seq2,nmod,M,gapscore,start,&s); }

char	*GapAlnTraceSMatrix2(Int4 a, Int4 b, Int4 len, unsigned char *seq2,
	Int4 nmod, smx_typ *M, Int4 **gapscore, Int4 *start,Int4 *score)
// return the trace, the start position of the trace, and the trace length.
{
	Int4    i,j,newlen,lenTrace;
	char	*operation,*tmp;

	if(gapscore)
	  tmp=gapped_aln_trace_smatrixSW(a,b,len,seq2,nmod,M,gapscore,
				&lenTrace,score,start); 
	else tmp=gap_aln_trace_smatrixSW(a,b,len,seq2,nmod,M,&lenTrace,score,start);
#if 0
std::cerr << "tmp traceback"; std::cerr << std::endl;
std::cerr << tmp; std::cerr << std::endl;
#endif
	for(i=1,newlen=lenTrace; tmp[i]=='i'; i++) newlen--;
	NEW(operation,newlen+3,char); operation[0]='E';
	for(i=1,j=lenTrace-1; i < newlen; i++,j--) operation[i]=tmp[j];
	operation[i]='E'; i++;
	operation[i]=0; free(tmp);
#if 0
std::cerr << "real traceback"; std::cerr << std::endl;
std::cerr << operation; std::cerr << std::endl;
std::cerr << *start; std::cerr << " = start\n";
#endif
	return operation;
}

Int4	AlnSeqSMatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 *mpos, Int4 cutoff)
{ return GappedAlnSeqSMatrixSW(a, b, n2, seq2, nmod, M, NULL); }

Int4	GappedAlnSeqSMatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore)
/// SW dynamic programming with gap function 
{
	Int4  alnscore,J;
	char  *operation=gapped_aln_seq_smatrixSW(a,b,n2,seq2,nmod,M,NULL,&J,&alnscore);
	free(operation); return alnscore;
}

char	*GapOperationsSMatrix(Int4 a, Int4 b, Int4 len, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore)
{ Int4    s,N; return gapped_aln_seq_smatrixSW(a,b,len,seq2,nmod,M,gapscore,&N,&s); }

void    PutGappedSeqAlnSMatrix(FILE *fp, char *operation, Int4 offset, Int4 lenseq,
        unsigned char *seq, Int4 nmod, smx_typ *M)
{ put_seqaln_smatrixSW(fp,operation,lenseq,seq,offset,strlen(operation),nmod, M); }

Int4	PutSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore)
{ return PutSeqAlnSMatrixSW(fp, a, b, n2, seq2, nmod, M, gapscore,0); }

Int4	PutSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore,Int4 offset)
{ return PutSeqAlnSMatrixSW(fp,a,b,n2,seq2,nmod,M,gapscore,offset,' '); }

Int4	PutSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore,Int4 offset,char mode)
{
	char	*operation;
	Int4	J,alnscore;

        operation=gapped_aln_seq_smatrixSW(a,b,n2,seq2,nmod,M,gapscore,&J,&alnscore);
	// fprintf(fp,"Put: %s\n\n",operation);
	put_seqaln_smatrixSW(fp, operation, n2, seq2, offset, J, nmod, M, mode);
	// fprintf(stderr,"operations = %s\n",operation);
	free(operation);
	return alnscore;
}

#if 1
Int4	PutSeqAlnSMatrixSW2(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
	UInt4 offset, smx_typ M, Int4 *start)
{
        Int4    o,i,j,r1,s,t,v,n1,score,*pos,p,p0;
	Int4	j2,s1,s0,g,gs,g_opt,J;
	unsigned char	*seq1;
	char	*out[3];
	a_type	A = SMatrixA(M);

	J = strlen(operation);
	/** get total length of profile **/
	n1 = LenSMatrix(M); 
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	MaxSegSMatrix(seq1, M);
	/*** 3. Trace back step. ***/
	MEW(out[0],J+3,char); MEW(out[1],J+3,char);
	MEW(out[2],J+3,char); NEW(pos,J+3,Int4); 

	// operations = i,I,d,D,m,M,E;
	fprintf(fp,"\n\n");
	for(p=p0=1,o=1,j=start[1],i=start[2]; operation[o] != 'E'; o++){ 
		switch(operation[o]){
		    case 'M': case 'm':
			v = ValSMatrix(i,seq2[j],M);
			pos[p] = j;
			out[1][p] = AlphaChar(seq1[i],A);
			out[2][p] = AlphaChar(seq2[j],A); 
			p0++;
			if(seq1[i]==seq2[j]) out[0][p] = ':';
			else if(v > 0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; i++; j++; 
			break;
		    case 'D': 
		    case 'd': // Gap ('-') is in sequence
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			p0++; pos[p]=j;
			p++; i++;  
			break;
		    case 'I':   // Gap ('-') is in profile
			out[1][p] = '-'; out[0][p] = ' ';
			out[2][p] = AlphaChar(seq2[j],A);
			pos[p] = j; p++; j++; break;
		    case 'i': 
#if 1
			out[1][p] = '-'; out[0][p] = ' ';
			out[2][p] = AlphaChar(seq2[j],A);
			pos[p] = j; p++; 
#endif
			j++; break;
			break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			print_error("this should not happen"); 
			break;
		}
	}
	fprintf(fp,"\nQuery: ");
	for(s=1; s < p ; s++) fprintf(fp,"%c",out[1][s]); 
	fprintf(fp,"\n       ");
	for(s=1; s < p ; s++) fprintf(fp,"%c",out[0][s]); 
	fprintf(fp,"\nSbjct: ");
	for(s=1; s < p ; s++) fprintf(fp,"%c",out[2][s]);
	fprintf(fp," %d-%d ",pos[1]+offset,j-1+offset); pos[0]=j-1;
	fprintf(fp,"(%d)\n",n2-pos[0]);
	/*** 5. free allocated memory ***/
	free(seq1); free(pos);
	free(out[0]); free(out[1]); free(out[2]);
	return o;
}

#endif

Int4	put_seqaln_smatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
	UInt4 offset, Int4 J, Int4 nmod, smx_typ *M)
{ return put_seqaln_smatrixSW(fp, operation, n2, seq2, offset, J, nmod, M,' '); }

Int4	put_seqaln_smatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
	UInt4 offset, Int4 J, Int4 nmod, smx_typ *M, char mode)
// mode == 'C' output cma style fake sequence; else regular output...
{
        Int4    o,m,i,j,k,r1,s,t,v,n1,score,*pos,p,p0;
	Int4	j2,s1,s0,g,gs,g_opt;
	unsigned char	*seq1,*seq;
	char	*out[3];
	short	*mtf;
	a_type	A = SMatrixA(M[1]);
	Int4	fake_len=0,true_len=0;
	char	*BigStr=0;
	Int4	BigLen=0;

	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	NEW(mtf, n1+3, short); 
	for(s=0, m=1; m <= nmod; m++){
	    MaxSegSMatrix(seq1+s, M[m]);
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; mtf[s] = m; }
	}
	fake_len=s;
	if(mode == 'C'){
		NEW(BigStr,strlen(operation)+10,char); BigLen=0;
	}
	/*** 3. Trace back step. ***/
	NEW(seq, J+3, unsigned char); 
	MEW(out[0],J+3,char); MEW(out[1],J+3,char);
	MEW(out[2],J+3,char); NEW(pos,J+3,Int4); 

	// operations = i,I,d,D,m,M,E;
	// fprintf(fp,"\n\n");
	m=1;  k=1; score=0;
	for(p=p0=1,o=j=i=1; operation[o] != 'E'; o++){ 
		switch(operation[o]){
		    case 'M': case 'm':
			v = ValSMatrix(k,seq2[j],M[m]); score+=v;
			pos[p] = j;
			out[1][p] = AlphaChar(seq1[i],A);
			out[2][p] = AlphaChar(seq2[j],A); 
			seq[p0]=seq2[j]; p0++;
			if(seq1[i]==seq2[j]) out[0][p] = ':';
			else if(v > 0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; i++; j++; k++;
			break;
		    case 'D': 
		    case 'd': // Gap ('-') is in sequence
			// score+=ValSMatrix(k,UndefAlpha(A),M[m]);
			score+=(Int4)floor(ExpScoreSMatrix(k,M[m])+0.5);
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			seq[p0]=0; p0++; pos[p]=j;
			p++; i++;  k++;
			break;
		    case 'I':   // Gap ('-') is in profile
			out[1][p] = '-'; out[0][p] = ' ';
			out[2][p] = tolower(AlphaChar(seq2[j],A));
			pos[p] = j; p++; j++; break;
		    case 'i': j++; break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			print_error("this should not happen"); 
			break;
		}
		if(k > LenSMatrix(M[m])){ // 4. Print out alignment 
if(mode=='C'){		// output a fake cma file entry
		   // if(m == 1)  fprintf(fp,"{()");
	   	   for(s=1; s < p ; s++){ 
		      if(isalpha(out[2][s])) true_len++;
		      // fprintf(fp,"%c",out[2][s]);
		      BigStr[BigLen]=out[2][s]; BigLen++;
		   }
} else {
		   if(m > 1)  fprintf(fp,"(%d)\n",pos[1]-pos[0]-1);
	   	   fprintf(fp,"\nmtf %c: ",(char)(m + 'A' - 1));
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[1][s]); 
	   	   fprintf(fp,"\n       ");
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[0][s]); 
		   fprintf(fp,"\n%5.1f: ", -log10(SMatrixProb(score,M[m])));
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[2][s]);
	   	   fprintf(fp," %d-%d [%d]",pos[1]+offset,j-1+offset,score);
	   	   // fprintf(fp," %d-%d ",pos[1]+offset,j-1+offset);
}
		   pos[0]=j-1;
		   if(m < nmod) { m++; k=1; p=p0=1; score=0; }
		   else break;
		}
	}
	if(mode=='C'){
	   fprintf(fp,"%d(%d):\n",true_len,fake_len); 
	   fprintf(fp,">fake goscan seq\n{()"); 
	   fprintf(fp,"%s()}*\n\n",BigStr); 
	   free(BigStr);
	} else fprintf(fp,"(%d)\n",n2-pos[0]);
	/*** 5. free allocated memory ***/
	free(seq1); free(mtf); free(pos);  free(seq);
	free(out[0]); free(out[1]); free(out[2]);
	return o;
}

Int4	PutFullSeqAlnSMatrixSW(FILE *fp, Int4 a, Int4 b, Int4 n2, 
	unsigned char *seq2, Int4 nmod, smx_typ *M, Int4 **gapscore)
{
	char	*operation;
        Int4    alnscore,J,o,m,i,j,k,r1,s,t,v,n1,score,*pos[3],p;
	Int4	j2,s1,s0,g,gs,g_opt,w=a+b;
	unsigned char	*seq1;
	char	*out[3];
	short	*mtf;
	a_type	A = SMatrixA(M[1]);

        operation=gapped_aln_seq_smatrixSW(a,b,n2,seq2,nmod,M,gapscore,&J,&alnscore);
	// fprintf(fp,"Put: %s\n\n",operation);
	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	NEW(mtf, n1+3, short); 
	for(s=0, m=1; m <= nmod; m++){
	    MaxSegSMatrix(seq1+s, M[m]);
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; mtf[s] = m; }
	}
	/*** 3. Trace back step. ***/
	MEW(out[0],J+3,char); MEW(out[1],J+3,char);
	MEW(out[2],J+3,char);
	NEW(pos[1],J+3,Int4); NEW(pos[2],J+3,Int4);

	m=nmod; k=LenSMatrix(M[m]);
	// operations = i,I,d,D,m,M,E;
	for(p=1,o=J,j=n2,i=n1; operation[o] != 'E'; o--){ 
		switch(operation[o]){
		    case 'M': case 'm':
			v = ValSMatrix(k,seq2[j],M[m]);
			pos[1][p] = i; pos[2][p] = j;
			out[1][p] = AlphaChar(seq1[i],A);
			out[2][p] = AlphaChar(seq2[j],A);
			if(seq1[i]==seq2[j]) out[0][p] = ':';
			else if(v > 0) out[0][p] = '.';
			else out[0][p] = ' ';
			p++; i--; j--; k--;
			break;
		    case 'D': case 'd': // Gap ('-') is in sequence
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			pos[1][p] = i; p++; i--;  k--;
			break;
		    case 'I': case 'i': // Gap ('-') is in profile
			out[1][p] = '-'; out[0][p] = ' ';
			out[2][p] = AlphaChar(seq2[j],A);
			pos[2][p] = j; p++; j--;
			break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			print_error("this should not happen"); 
			break;
		}
		if(k==0){ if(m > 1){ m--; k=LenSMatrix(M[m]); } else m--; }
	}
	/** 4. Print out alignment **/
	fprintf(fp,"\n\n");
	for(i=p-1; i > 0 ; ){
	   v = MINIMUM(Int4,50,i); fprintf(fp,"     ");
	   for(j=i; j > i-v ; j--) {
		if(pos[1][j]) {
		   fprintf(fp,"%c",(char)(mtf[pos[1][j]] + 'A' - 1));
		} else fprintf(fp," ");
	   }
	   fprintf(fp,"\n     ");
	   for(j=i; j > i-v ; j--) fprintf(fp,"%c",out[1][j]); 
	   fprintf(fp," concensus\n     ");
	   for(j=i; j > i-v ; j--) fprintf(fp,"%c",out[0][j]); fprintf(fp,"\n     ");
	   for(j=i; j > i-v ; j--) fprintf(fp,"%c",out[2][j]);
	   fprintf(fp," sequence\n  ");
	   for(j=i; j > i-v ; j--) {
		if(pos[2][j] && pos[2][j] % 10 == 0) {
			fprintf(fp,"%4d",pos[2][j]); j-=3;
		} else fprintf(fp," ");
	   }
	   fprintf(fp,"\n\n"); i-=v;
	}
	/*** 5. free allocated memory ***/
	free(seq1); free(mtf);
	free(out[0]); free(out[1]); free(out[2]);
	free(pos[1]); free(pos[2]);
	free(operation);
	return alnscore;
}

Int4	ScoreGappedAlnSeqSmatrixSW(Int4 a, Int4 b, Int4 n2, unsigned char *seq2, 
	Int4 nmod, smx_typ *M, Int4 **gapscore, BooLean local)
{
	Int4	m,i,j,k,r1,s,t,v,n1,score,max_i,max_j,mscore,maxscore=0;
	Int4	**MAT,**DEL,**INS;
	Int4	j2,s1,s0,g,gs,g_opt,w=a+b,J0;
	unsigned char	*seq1;
	short	*mtf;
	a_type	A = SMatrixA(M[1]);
	Int4	gINS,*gDEL,t0;

	if(gapscore != NULL) gapscore--; // want gapscore[m-1] for block m.
	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	NEW(mtf, n1+3, short); 
	for(s=0, m=1; m <= nmod; m++){
	    MaxSegSMatrix(seq1+s, M[m]);
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; mtf[s] = m; }
	}
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,n1+3,Int4*); 
	MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		MEW(MAT[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
	}
	// Make GLOBAL ALIGNMENT with respect to profile,
	MAT[0][0] = 0; DEL[0][0] = 0; INS[0][0] = 0;
	for(i=1; i<= n1; i++) {  
		DEL[i][0] = DEL[i-1][0] - b;
		INS[i][0] = INS[i-1][0] - b;
		if(local) MAT[i][0] = 0;  // New: Local on ends of profile 
		else MAT[i][0] = INS[i-1][0];  // Global for profile 
	}
	for(j=1; j<= n2; j++) { // Make LOCAL with respect to sequence.
		MAT[0][j] = 0; // traceback full alignment 
		DEL[0][j] = INS[0][j] = SHRT_MIN;
	}
	// 2. Dynamic programming step. 
	NEW(gDEL,n2+3,Int4);
	for(k=m=i=1; i<= n1; i++) {
	   gINS=0;
	   if(k==LenSMatrix(M[m])) {
	     for(j=1; j<= n2; j++) {  // Eliminate insert state...
                t=0;
	 	s = MAT[i-1][j-1] + ValSMatrix(k,seq2[j],M[m]);
		if(gapscore == NULL){ // no affine penalty for gapscore.
                   if((s0=MAT[i][j-1]) > (s1=INS[i][j-1])){
                        gINS=-1;  INS[i][j] = s0;
                   } else { gINS--; INS[i][j] = s1; }
                   if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }
		}

                if((s0=MAT[i-1][j]-w) > (s1=DEL[i-1][j]-b)){
                        gDEL[j]=1;  DEL[i][j] = s0;
                } else { gDEL[j]++;  DEL[i][j] = s1; }
                if(s < DEL[i][j]){ s=DEL[i][j]; t=gDEL[j]; }

		MAT[i][j] = s;
	        if(local && maxscore < s) maxscore=s;
	     }
	   } else if(k==1 && m > 1 && gapscore != NULL) {
	     for(j=1; j<= n2; j++) {  // Use gap function here.
		// Don't jump over gap function regions! Can't use t=gDEL[j].
                DEL[i][j] = MAXIMUM(Int4,MAT[i-1][j]-w,DEL[i-1][j]-b);
                s = DEL[i][j]; t=1;

                if((s0=MAT[i][j-1]-w) > (s1=INS[i][j-1]-b)){
                        gINS=-1;  INS[i][j] = s0;
                } else { gINS--; INS[i][j] = s1; }
                if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }

		v = ValSMatrix(k,seq2[j],M[m]);
		for(s1=SHRT_MIN,j2=j-1, g=0; j2 > 0; g++,j2--){  
		    if((gs = gapscore[m][g]) == SHRT_MIN) break; // over max gap?
		    s0 = MAT[i-1][j2] + v + gs;
		    if(s1 < s0){ s1 = s0; g_opt = g; }
		}
                if(s < s1) { s = s1; t=0; }
		MAT[i][j] = s;
	        if(local && maxscore < s) maxscore=s;
	     }
	   } else { 		      // Use affine gap penalty
	      for(j=1; j<= n2; j++) {
                t=0; s = MAT[i-1][j-1] + ValSMatrix(k,seq2[j],M[m]);

                if((s0=MAT[i-1][j]-w) > (s1=DEL[i-1][j]-b)){
                        gDEL[j]=1;  DEL[i][j] = s0;
                } else { gDEL[j]++;  DEL[i][j] = s1; }
                if(s < DEL[i][j]){ s=DEL[i][j]; t=gDEL[j]; }

                if((s0=MAT[i][j-1]-w) > (s1=INS[i][j-1]-b)){
                        gINS=-1;  INS[i][j] = s0;
                } else { gINS--; INS[i][j] = s1; }
                if(s < INS[i][j]){ s=INS[i][j]; t=gINS; }

                MAT[i][j] = s;
	        if(local && maxscore < s) maxscore=s;
	      }
	   }
	   k++;
	   if(k > LenSMatrix(M[m])) { k = 1; m++; }
	}
	// 2b. Find optimum global score with respect to profile.
    if(local){ 
	score = maxscore; 
    }
    else {
	max_i=n1;		// global with respect to profile.
	for(score=INT4_MIN, j=1; j<= n2; j++) {
		if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
	}
    }
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(MAT[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(DEL); free(INS); 
	free(seq1); free(mtf); free(gDEL);
        return score;
}

//======================= FAST gapped align for sampling ================
char	*gap_aln_trace_smatrixSW(Int4 a,Int4 b,Int4 n2,unsigned char *seq2, 
	Int4 nmod,smx_typ *M,Int4 *J,Int4 *alnscore,Int4 *start)
{
	Int4	m,i,j,k,r1,s,t,v,n1,score,max_i,max_j,mscore;
	Int4	**MAT,**T,**DEL,**INS;
	Int4	j2,s1,s0,g,gs,g_opt,w=a+b,J0;
	a_type	A = SMatrixA(M[1]);
	Int4	gINS,t0;
	char	*operation;

	/** get total length of profile **/
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/*** 1. Allocate and initialize memory. ***/
	MEW(MAT,n1+3,Int4*); MEW(T,n1+3,Int4*);
	MEW(DEL,n1+3,Int4*); MEW(INS,n1+3,Int4*);
	for(i=0; i<= n1; i++) { 
		MEW(MAT[i],n2+3,Int4); MEW(T[i],n2+3,Int4); 
		MEW(DEL[i],n2+3,Int4); MEW(INS[i],n2+3,Int4); 
	}
	// Make GLOBAL ALIGNMENT with respect to profile,
	MAT[0][0] = 0; DEL[0][0] = 0; INS[0][0] = 0;
	for(i=1; i<= n1; i++) {  
		DEL[i][0] = DEL[i-1][0] - b;
		INS[i][0] = INS[i-1][0] - b;
		MAT[i][0] = INS[i-1][0];  // Local on ends of profile 
		T[i][0] = 1; // for full alignment
	}
	for(j=1; j<= n2; j++) { // Make LOCAL with respect to sequence.
		MAT[0][j] = 0; T[0][j]=-1; // traceback full alignment 
		DEL[0][j] = INS[0][j] = SHRT_MIN;
	}
	Int4 *gdel; NEW(gdel,n2+3,Int4);
	// 2. Dynamic programming step. 
	for(k=m=i=1; i<= n1; i++) {
	   gINS=0;
	   Int4 *mat_im1=MAT[i-1]; 
	   Int4 *mat_i=MAT[i];
	   Int4 *ins_i=INS[i];
	   Int4 *del_i=DEL[i];
	   Int4 *del_im1=DEL[i-1];
	   register Int4 jj,jm1;
	   smx_typ  smx=M[m];
	   if(k==LenSMatrix(smx)) {
	     for(jm1=0,jj=1; jj<= n2; jj++,jm1++) {  // Eliminate insert state...
                t=0; s = mat_im1[jm1] + ValSMatrix(k,seq2[jj],smx);
                if((s0=mat_i[jm1]) > (s1=ins_i[jm1])){
                        gINS=-1;  ins_i[jj] = s0;
                } else { gINS--; ins_i[jj] = s1; }
                if(s < ins_i[jj]){ s=ins_i[jj]; t=gINS; }

                if((s0=mat_im1[jj]-w) > (s1=del_im1[jj]-b)){
                        gdel[jj]=1;  del_i[jj] = s0;
                } else { gdel[jj]++;  del_i[jj] = s1; }
                if(s < del_i[jj]){ s=del_i[jj]; t=gdel[jj]; }

		T[i][jj] = t; mat_i[jj] = s;
	     }
	   } else { 		      // Use affine gap penalty
	      for(jm1=0,jj=1; jj<= n2; jj++,jm1++) {
                t=0; s=mat_im1[jm1] + ValSMatrix(k,seq2[jj],smx);
                if((s0=mat_im1[jj]-w) > (s1=del_im1[jj]-b)){
                        gdel[jj]=1;  del_i[jj] = s0;
                } else { gdel[jj]++;  del_i[jj] = s1; }
                if(s < del_i[jj]){ s=del_i[jj]; t=gdel[jj]; }

                if((s0=mat_i[jm1]-w) > (s1=ins_i[jm1]-b)){
                        gINS=-1;  ins_i[jj] = s0;
                } else { gINS--; ins_i[jj] = s1; }
                if(s < ins_i[jj]){ s=ins_i[jj]; t=gINS; }
                T[i][jj] = t; mat_i[jj] = s;
	      }
	   } k++;
	   if(k > LenSMatrix(M[m])) { k = 1; m++; }
	} free(gdel);
	// 2b. Find optimum global score with respect to profile.
	max_i=n1;		// global with respect to profile.
	for(score=INT4_MIN, j=1; j<= n2; j++)
		if((s=MAT[max_i][j]) > score){ score=s; max_j=j; }
	/*** 3. Trace back step. ***/
	// operations: 
	// 'i' = insertion in sequence outside of motifs
	// 'I' = insert in sequence within motif (need to delete this)
	//  'M' = match to start of a motif block
	//  'm' = match to other sites in motif
	//  'D' = deletion of sequence within motif
	//  'd' = deletion of sequence outside of motif (not used)
	NEW(operation,n1+n2+3,char); operation[0]='E';
	m=nmod; k=LenSMatrix(M[m]);
	for(J0=0,j=n2; j > max_j; ){ J0++; operation[J0] = 'i'; j--; }
	// for(i=n1; i > 0 || j > 0; ){  // full alignment...
	for(i=n1; i > 0; ){  // full alignment...Stop at end of profile
	  t0 = T[i][j];
	  do {
	    if(t0 > 0){ t=1; t0--; }
            else if(t0 < 0){ t=-1; t0++; } else { t=0; }
	    switch(t){
		case 0:
		   i--; J0++;
		   if(k==1) operation[J0] = 'M'; else operation[J0] = 'm';
		   j--; k--; break;
		case 1:  // Gap ('-') is in sequence; add X's for gaps
			J0++; 
			if(k==1) operation[J0] = 'D'; else operation[J0] = 'd';
			i--;  k--; break;
		case -1:	// Gap ('-') is in profile
			if(m < 1 || k==LenSMatrix(M[m]))
				{ J0++; operation[J0] = 'i'; }
			else { J0++; operation[J0] = 'I'; }
			j--; break;
		default: print_error("this should not happen"); break;
	     }
	     if(k==0){ m--; if(m > 0) k=LenSMatrix(M[m]); }
	  } while(t0 != 0);
	}
	*start=j+1; J0++; operation[J0]='E';
	/*** 5. free allocated memory ***/
	for(i=0; i<=n1; i++) {free(MAT[i]);free(T[i]);free(DEL[i]);free(INS[i]);}
	free(MAT); free(T); free(DEL); free(INS); 
	*alnscore = score; *J = J0;
        return operation;
}

//======================= END FAST gapped align for sampling ================

Int4	put_cmaseq_smatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
	UInt4 offset, Int4 J, Int4 nmod, smx_typ *M)
{
        Int4    o,m,i,j,k,r1,s,t,n1,*pos,p,p0;
	Int4	j2,s1,s0,g,gs,g_opt;
	Int4	real,fake;
	unsigned char	*seq;
	char	c,*out;
	a_type	A = SMatrixA(M[1]);

	/** get concensus sequence for smatrix **/
	for(s=0, m=1; m <= nmod; m++){
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; }
	}
	/*** 3. Trace back step. ***/
	NEW(seq, J+3, unsigned char); 
	MEW(out,J+3,char); NEW(pos,J+3,Int4); 

	// operations = i,I,d,D,m,M,E;
	// fprintf(fp,"\n\n");
	m=1;  k=1; real=fake=0;
	for(p=p0=1,o=j=i=1; operation[o] != 'E'; o++){ 
		switch(operation[o]){
		    case 'M': case 'm':
			pos[p] = j;
			out[p] = AlphaChar(seq2[j],A); 
			real++; fake++;
			seq[p0]=seq2[j]; p0++;
			p++; i++; j++; k++;
			break;
		    case 'D': 
		    case 'd': // Gap ('-') is in sequence
			out[p] = '-';
			fake++;
			seq[p0]=0; p0++; pos[p]=j;
			p++; i++;  k++;
			break;
		    case 'I':   // Gap ('-') is in profile
			c = AlphaChar(seq2[j],A);
			out[p] = tolower(c);
			real++;
			pos[p] = j; p++; j++; break;
		    case 'i': j++; break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			print_error("this should not happen"); 
			break;
		}
		if(k > LenSMatrix(M[m])){ // 4. Print out alignment 
		   // if(m > 1)  fprintf(fp,"(%d)\n",pos[1]-pos[0]-1);
	   	   fprintf(fp,"\n$%d=%d(%d):\n",m,real,fake); real=fake=0;
	   	   fprintf(fp,">xxx {|%d(0)|}\n",pos[1]+offset-1); 
		   pos[0]=j-1;
		   fprintf(fp,"{()");
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[s]);
		   fprintf(fp,"()}*\n");
		   if(m < nmod) { m++; k=1; p=p0=1; }
		   else break;
		}
	}
	fprintf(fp,"(%d)\n",n2-pos[0]);
	/*** 5. free allocated memory ***/
	free(pos);  free(seq); free(out);
	return o;
}

//============== Aleksandar's Complex alignment routines ================

static void	complex_affine_aln_smx(Int4 *matj, Int4 *matjm1, 
	   Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
           Int4  delo, Int4  dele, Int4  inso, Int4 ie, 
	   Int4 seq_len, unsigned char *seq, Int4 *scorej)
{
   register Int4	i,im1,s,Mjm1_im1,Djm1_im1,Ijm1_im1;

   for(im1=0,i=1;i<=seq_len;im1++,i++){
	Ijm1_im1 = insjm1[im1]; Djm1_im1 = deljm1[im1]; Mjm1_im1 = matjm1[im1];

        if(Mjm1_im1>=Djm1_im1 && Mjm1_im1>=Ijm1_im1) matj[i] = scorej[seq[i]]+Mjm1_im1;
        else if(Djm1_im1 >= Ijm1_im1) matj[i]=scorej[seq[i]]+Djm1_im1;
        else matj[i] = scorej[seq[i]]+Ijm1_im1;

        if((s=delo+matjm1[i]) > deljm1[i]) delj[i]=s+dele;
	else delj[i]=deljm1[i]+dele;

        if((s=inso+matj[im1]) > insj[im1]) insj[i]=s+ie;
	else insj[i]=insj[im1]+ie;
   }
}

static void	complex_affine_aln_smx(Int4 *matj, Int4 *matjm1, 
	   Int4 *insj, Int4 *insjm1, Int4 *delj, Int4 *deljm1,
           Int4  delo, Int4  dele, Int4  inso, Int4 ie, 
	   Int4 seq_len, unsigned char *seq, Int4 *scorej, Int4 *grease_pen)
// about 80% of compute time is spent here...
{
   register Int4	i,im1,s,Mjm1_im1,Djm1_im1,Ijm1_im1;

   for(im1=0,i=1;i<=seq_len;im1++,i++){
	Ijm1_im1 = insjm1[im1]; Mjm1_im1 = matjm1[im1]; Djm1_im1 = deljm1[im1]; 

        if(Mjm1_im1>=Djm1_im1 && Mjm1_im1>=Ijm1_im1) matj[i] = scorej[seq[i]]+Mjm1_im1;
        else if(Djm1_im1 >= Ijm1_im1) matj[i]=scorej[seq[i]]+Djm1_im1;
        else matj[i] = scorej[seq[i]]+Ijm1_im1;

        if((s=delo+matjm1[i]) > deljm1[i]) delj[i]=s+dele;
	else delj[i]=deljm1[i]+dele;

        if((s=inso+matj[im1]) > insj[im1]) insj[i]=s+grease_pen[i]+ie;
	else insj[i]=insj[im1]+grease_pen[i]+ie;
   }
}

static void	complex_gap_align_smx(Int4 *matj, register Int4 *matjm1, 
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

static char	*get_traceback_operation_smx2(Int4 seq_len, Int4 prof_len,
	Int4 nblks, Int4 *block_lengths, Int4 *start_prof, Int4 **MAT, 
	Int4 **DEL, Int4 **INS, Int4 *oper_len, Int4 *grease_pen, 
	Int4 *alignscore, Int4 *start, idp_typ *idp)
// Traceback without gpen...
{
        Int4    l,s,i,i0,j,k,bl_nmbr,best_i,inse,mxm, g;
        char    state,*operation, *back_operation,flag;
        Int4	*io=idp->InsOpen();
        Int4	*ie=idp->InsExtend();
        Int4	*od=idp->DelOpen();
        Int4	*de=idp->DelExtend();

        NEW(operation,seq_len+prof_len+3,char);
        NEW(back_operation,seq_len+prof_len+3,char);
        j=prof_len; i=seq_len;

        for(state='E',k=0;state != 'X'; ){
           switch(state){
                case 'B': // begin state (exit traceback).
                        back_operation[k++]='E'; *start=i; state='X'; break;
                case 'E': // end state.  
                        back_operation[k++]='E'; 
                        bl_nmbr=nblks; s=start_prof[bl_nmbr]; flag='T';
                        mxm=MAT[j][i];state='M'; i0=i;
                        for(i=1;i<=seq_len;i++){
                           if(MAT[j][i]>mxm){ state='M'; mxm=MAT[j][i];i0=i;}
                           if(DEL[j][i]>mxm){ state='D'; mxm=DEL[j][i];i0=i;}
// fprintf(stderr,"%d: max score = %d\n",i,mxm);
                        } i=i0; *alignscore=mxm;
                        break;
                case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        if(j==s){
                           bl_nmbr--; s=start_prof[bl_nmbr]; flag='T';
                           j--; back_operation[k++]='M'; 
                        } else { // else we are within a block & use affine gaps.
                           j--; i--; back_operation[k++]='m';
			}
                        if(MAT[j][i]>INS[j][i] && MAT[j][i]>DEL[j][i]) state='M';
                        else if(INS[j][i]>DEL[j][i]) state='I';
			else state='D'; 
                        break;
                case 'D': // previously sampled deletion
                        if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                        if(j==s){
                            j--; back_operation[k++]='D'; flag='T';
                            bl_nmbr--; s=start_prof[bl_nmbr];
                        } else { j--; back_operation[k++]='d'; }
#if 1	// OLD... od and de here but not up there...
                        if(DEL[j][i]+de[j] >= MAT[j][i]+od[j]+de[j]) state='D';
                        else state='M';
#endif
#if 0	// AFN...  allow insertions after deletions and vise versa (below)
                        if(DEL[j][i]+de[j] <= MAT[j][i]+od[j]+de[j]) state='M';
                        else if(INS[j][i] +io[j]+ie[j] > DEL[j][i] +de[j]) state='I';
                        else state='D';
#endif
                        break;
                case 'I': // previously sampled insertion.
                        if(j==start_prof[bl_nmbr]+block_lengths[bl_nmbr]-1){
			    if(flag=='F') back_operation[k++]='i';  flag='F';
			} else back_operation[k++]='I'; 
			i--;
                        if(grease_pen){
                            if(j==start_prof[bl_nmbr]+block_lengths[bl_nmbr]-1) {inse=0;}
                            else { inse = ie[j] + grease_pen[i]; }
                        } else inse = ie[j];
                        if (INS[j][i]+inse > MAT[j][i]+io[j]+inse) state='I';
                        else state='M';
                        break;                   
                default:  print_error("this should not happen");
           }
        } *oper_len = k;
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
	return operation;
}

static char	*get_traceback_operation_smx0(Int4 seq_len, Int4 prof_len,
	Int4 nblks, Int4 *block_lengths, Int4 *start_prof, Int4 **MAT, 
	Int4 **DEL, Int4 **INS, Int4 *oper_len, Int4 *grease_pen, 
	Int4 *alignscore, Int4 *start, idp_typ *idp)
{
        Int4    l,s,i,i0,j,k,bl_nmbr,best_i,inse;
        Int4    *gapfunct,gapend,i_gap,gapstart;
        Int4    mxm, gap,g;
        char    state,*operation, *back_operation;
        Int4 *io=idp->InsOpen();
        Int4 *ie=idp->InsExtend();
        Int4 *od=idp->DelOpen();
        Int4 *de=idp->DelExtend();
        Int4    **gpen=idp->GapFunct();

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
// fprintf(stderr,"%d: max score = %d\n",i,mxm);
                        } i=i0;*alignscore=mxm;
                        break;
                case 'M': // previously sampled matched state.
                        if(j <= 1) { back_operation[k++]='M'; state='B'; break; }
                        if(j==s){
                           bl_nmbr--; s=start_prof[bl_nmbr];
                           j--; back_operation[k++]='M'; 
                           state='G'; gapstart=i-1;  
                        } else { // else we are within a block & use affine gaps.
                           j--; i--; back_operation[k++]='m';
                           if(MAT[j][i]>INS[j][i] && MAT[j][i]>DEL[j][i]) state='M';
                           else if(INS[j][i]>DEL[j][i]) state='I';
			   else state='D'; 
                        } break;
                case 'D': // previously sampled deletion
                        if(j<=1){ back_operation[k++]='D'; state='B'; i++; break; }
                        if(j==s){
                            j--; back_operation[k++]='D'; 
                            bl_nmbr--; s=start_prof[bl_nmbr];
                            state='G'; gapstart=i; 
                        } else { j--; back_operation[k++]='d';
                             if(DEL[j][i]+de[j]>=MAT[j][i]+od[j]+de[j]) state='D';
                             else state='M';
                        } break;
                case 'I': // previously sampled insertion.
                          back_operation[k++]='I'; i--;
                          if(grease_pen){
                                if(j==start_prof[bl_nmbr]+block_lengths[bl_nmbr]-1) {inse=0;}
                                else { inse = ie[j] + grease_pen[i]; }
                          } else inse = ie[j];
                          if (INS[j][i]+inse > MAT[j][i]+io[j]+inse) state='I';
                          else state='M'; 
                          break;                   
                case 'G':  // use gap function (between blocks).
                  if(gpen[bl_nmbr]){
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
                  } else {
                      gapend = i_gap = gapstart;
                      mxm=MAT[j][i_gap];
                      for(g=0; g <= gapend; g++,i_gap--){
                         gap=MAT[j][i_gap];
                         if(MAT[j][i_gap]>=mxm) {best_i=i_gap; state='M';mxm=gap;}
                         gap=DEL[j][i_gap];
                         if(DEL[j][i_gap]>=mxm) {best_i=i_gap; state='D';mxm=gap;}
                      } g=(gapstart-best_i);
                  }
                  for(l=0;l < g; l++) { back_operation[k++]='i';} i =best_i; 
                  break;
                default:  print_error("this should not happen");
           }
        } *oper_len = k;		// for alex's code: *oper_len = k-2;
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
	return operation;
}

char    *ComplexAlnSMatrix0(idp_typ *idp, Int4 seq_len, unsigned char *seq,
        Int4 Rpts, smx_typ *M, Int4 *Score, Int4 *Start, Int4 *Oper_len)
// From Alex...
{
        Int4    pen,i,j,k,l,im1,jm1,i0,n,r,s,t,prof_len;
        Int4    maxMAT,maxDEL,alph_len,total,g,s1,s0;
        Int4    *block_lengths, *start_prof;
        Int4    **MAT,**DEL,**INS,**score;
        Int4    nblks = Rpts*idp->nBlks();
         
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

        Int4 *io=idp->InsOpen();
        Int4 *ie=idp->InsExtend();
        Int4 *od=idp->DelOpen();
        Int4 *de=idp->DelExtend();
        Int4    **gpen=idp->GapFunct();

        DEL[0][0]=INS[0][0]=MAT[0][0]=0;
        DEL[1][0]=INS[1][0]=MAT[1][0]=od[1]+de[1];
        for(s=2,jm1=1,j=2;j<=prof_len;jm1++,j++) {
                if(j!=start_prof[s]){ MAT[j][0]=INS[j][0]=DEL[j][0]=DEL[jm1][0]+de[j];
                } else {
                   if(gpen[s-1]){
                        MAT[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0]; 
                        DEL[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
                        INS[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
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
           if(J!=start_prof[s+1] || gpen[s]==0){
// if(J==End) fprintf(stderr,"%d: io = %d; ie = %d\n",J,io[J],ie[J]);
              if((grease_pen && J!=End)){
                complex_affine_aln_smx(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J],grease_pen);
              } else { // don't use grease penalities between blocks...
                complex_affine_aln_smx(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J]);
              }
           } else { // impose gap penalties between blocks when at first position...
              complex_gap_align_smx(MAT[J],MAT[Jm1],INS[J],DEL[J],DEL[Jm1],
                od[J]+de[J],io[J],ie[J],seq_len,seq,score[J],gpen[s],grease_pen);
           }
           if(J==start_prof[s+1]) s++;
        }
	char *operation=get_traceback_operation_smx0(seq_len,prof_len,nblks,
		block_lengths,start_prof,MAT,DEL,INS,Oper_len,grease_pen, 
		Score, Start, idp);
	free(block_lengths);free(start_prof);
        for(i=0;i<=prof_len;i++){free(MAT[i]);free(DEL[i]);free(INS[i]);}
        free(MAT);free(DEL);free(INS);
        for(i=0;i<=prof_len;i++) free(score[i]);
        free(score);
        if(grease_pen) free(grease_pen);
// std::cerr << operation; std::cerr << std::endl;
// fprintf(stderr,"Score = %d\n",*Score);
        return operation;
}

static char *get_traceback_operation_smx(Int4 seq_len, Int4 prof_len,
        Int4 nblks, Int4 *block_lengths, Int4 *start_prof, Int4 **MAT, 
        Int4 **DEL, Int4 **INS, Int4 *oper_len, Int4 *grease_pen, 
        Int4 *alignscore, Int4 *start, idp_typ *idp)
{
        Int4    l,s,i,i0,j,k,bl_nmbr,best_i,inse,end;
        Int4    *gapfunct,gapend,i_gap,gapstart;
        Int4    mxm, gap,g;
        char    state,*operation, *back_operation;
        Int4 	*io=idp->InsOpen();
        Int4 	*ie=idp->InsExtend();
        Int4 	*od=idp->DelOpen();
        Int4 	*de=idp->DelExtend();
        Int4 	**gpen=idp->GapFunct();

        NEW(back_operation,seq_len+prof_len+3,char);

        for(j=prof_len,i=seq_len,state='E',k=0;state != 'X'; ){
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
        NEW(operation,k+3,char);
        for(i=0;i<k;i++) operation[i]=back_operation[k-i-1];
        free(back_operation);
        return operation;
}

char *ComplexAlnSMatrix(idp_typ *idp, Int4 seq_len, unsigned char *seq,
        Int4 Rpts, smx_typ *M, Int4 *Score, Int4 *Start, Int4 *Oper_len)

{
        Int4    pen,i,j,k,l,im1,jm1,i0,n,r,s,t,prof_len;
        Int4    maxMAT,maxDEL,alph_len,total,g,s1,s0;
        Int4    *block_lengths, *start_prof;
        Int4    **MAT,**DEL,**INS,**score;
        Int4    nblks = Rpts*idp->nBlks();
         
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

        Int4 *io=idp->InsOpen();
        Int4 *ie=idp->InsExtend();
        Int4 *od=idp->DelOpen();
        Int4 *de=idp->DelExtend();
        Int4 **gpen=idp->GapFunct();

        DEL[0][0]=INS[0][0]=MAT[0][0]=0;
        DEL[1][0]=INS[1][0]=MAT[1][0]=od[1]+de[1];
        for(s=2,jm1=1,j=2;j<=prof_len;jm1++,j++) {
                if(j!=start_prof[s]){ MAT[j][0]=INS[j][0]=DEL[j][0]=DEL[jm1][0]+de[j];
                } else {
                   if(gpen[s-1]){
                        MAT[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0]; 
                        DEL[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
                        INS[j][0]=DEL[jm1][0]+od[j]+de[j]+gpen[s-1][0];
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
           if(J!=start_prof[s+1] || gpen[s]==0){
              if((grease_pen && J!=End)){
                complex_affine_aln_smx(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J],grease_pen);
              } else { // don't use grease penalities between blocks...
                complex_affine_aln_smx(MAT[J],MAT[Jm1],INS[J],INS[Jm1],DEL[J],DEL[Jm1],
                  od[J],de[J],io[J],ie[J],seq_len,seq,score[J]);
              }
           } else { // impose gap penalties between blocks when at first position...
              complex_gap_align_smx(MAT[J],MAT[Jm1],INS[J],DEL[J],DEL[Jm1],
                od[J]+de[J],io[J],ie[J],seq_len,seq,score[J],gpen[s],grease_pen);
           }
           if(J==start_prof[s+1]) s++;
        }
        char *operation=get_traceback_operation_smx(seq_len,prof_len,nblks,
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

