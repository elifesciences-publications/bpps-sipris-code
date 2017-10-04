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

smx_typ	MkSMatrixN(Int4 N, Int4 K, double *freq, a_type A)
{	
	double  p,nsd;

        p = 10.0/(double)N;  p = MINIMUM(double,p,0.05);
        for(nsd=2.0; nsd <= 5.0; nsd += 0.1){
             if(p >= (0.5*gammq(0.5,(nsd*nsd*0.5)))) { nsd -= 0.1; break; }
        }
	/*** fprintf(stderr,"\nN = %d; nsd = %g; p = %g\n",N,nsd,p); /****/
	return MkSMatrix(nsd,K,freq,A);
}

Int4    **SMatrix2GPSI(smx_typ smx, Int4 neg_inf)
// Call as: SMatrix2GPSI(smx, GBLAST_SCORE_MIN);
{
	a_type	A = smx->A;
        Int4 **posMatrix,i,j,r,length=smx->K;
        NEWP(posMatrix, length+3, Int4);
        for(j=0,i=1; j < length; j++,i++){
	   NEW(posMatrix[j],nAlpha(A)+3,Int4);
           for(r=0; r <= nAlpha(A); r++){
                posMatrix[j][r] = ValSMatrix(i,r,smx);
           }
        } NEW(posMatrix[j],nAlpha(A)+3,Int4);
	for(r = 0; r <= nAlpha(A); r++) { posMatrix[j][r] = neg_inf; }
	return posMatrix;
}

smx_typ	GPSI2SMatrix(Int4 length, int **posMatrix, a_type A)
{
	Int4	i,j,k,r,s;

	assert(posMatrix);
	smx_typ M=MkSMatrix(2.0,length,blosum62freq,A);
	for(j=0,i=1; j < length; j++,i++){
	   for(r=0; r <= nAlpha(A); r++){ s=(Int4) posMatrix[j][r]; SetSMatrix(r,i,s,M); }
	} return M;
}

smx_typ ReverseSMatrix(smx_typ smx)
// create and return the reverse of matrix smx.
{
	Int4	i,j,r,s;
	smx_typ smx2 = MkSMatrix(smx->nsd, smx->K, smx->freq, smx->A);
	for(j=1,i=smx->K; j <= smx->K; j++,i--){
	   for(r=0; r <= nAlpha(smx->A); r++){
		s=ValSMatrix(j,r,smx); SetSMatrix(r,i,s,smx2);
	   }
	} return smx2;
}

smx_typ	MkSMatrix(double nsd, Int4 K, double *freq, a_type A)
/* create a scoring matrix with all zero scores. */
{
	smx_typ	M;
	Int4	i,c;
	
	NEW(M,1,smatrix_type);
	M->A = A;
	M->posMatrix = NULL;
	M->qseq = NULL;
	M->nlet = nAlpha(A);
	M->K = K;
	M->nsd = nsd;
	M->neginf = (INT4_MIN/K) + 1;	/*** NEW ***/
	M->calc_prob = TRUE;
	M->calc_stats = TRUE;
	M->changed = TRUE;
	M->f0 = NULL;
	NEWP(M->score,K+2,Int4);
	for(i=1; i<=K;i++) {
		NEW(M->score[i],M->nlet+3,Int4);
		for(c=0; c<=M->nlet; c++) M->score[i][c]=0;
		M->score[i][c]=SHRT_MIN;  /** for blast matrix **/
	}
	NEW(M->max,K+2,Int4); NEW(M->min,K+2,Int4); 
	NEW(M->cmax,K+2,Int4); NEW(M->cmin,K+2,Int4);
	NEW(M->freq,M->nlet +2,double);
	for(c=0; c<=M->nlet; c++) M->freq[c]=freq[c];
	return M;
}

void	NilSMatrix(smx_typ M)
{
	Int4	i;

	for(i=1; i <= M->K;i++) { free(M->score[i]); }
	free(M->score); 
	if(M->f0 != NULL) free(M->f0);
	if(M->posMatrix != NULL){
		i = M->posMtrxLen;
		free(M->posMatrix[0]); free(M->posMatrix[i+1]);
		free(M->posMatrix);
		if(M->qseq != NULL) free(M->qseq);
	}
	free(M->cmin); free(M->cmax);
	free(M->min); free(M->max);
	free(M->freq); 
	free(M);
}

void	SetSMatrix(Int4 r, Int4 row , Int4 score, smx_typ M)
/* set the matrix for residue r in row equal to score */
{
	if(row < 1 || row > M->K) smatrix_error("input error");
#if 0
/*****
if(row<4) fprintf(stderr,"score[%d][%d] = %d; neginf=%d\n", row,r,score,M->neginf);
/*****/
#endif
	M->score[row][r]= MAXIMUM(Int4,score,M->neginf); 
	M->changed = TRUE;
}

double	ExpScoreSMatrix(double **freq, smx_typ M)
/** return the expected score for smatrix M **/
{
	Int4	j,d;
	double	score;
	a_type	A = M->A;

	for(score=0.0, j=1; j <= M->K; j++){
		for(d=0; d<=nAlpha(A); d++){
			score += freq[j][d]*M->score[j][d];
		}
	}
	return score;
}

double	ExpScoreSMatrix(Int4 j,smx_typ M)
/** return the expected score for smatrix M **/
{
	Int4	d;
	double	score;
	a_type	A = M->A;

	for(score=0.0,d=1; d<=nAlpha(A); d++){
		score += M->freq[d]*M->score[j][d];
	} return score;
}

double	ExpScoreSMatrix(smx_typ M)
/** return the expected score for smatrix M **/
{
	Int4	j,d;
	double	score;
	a_type	A = M->A;

	for(score=0.0, j=1; j <= M->K; j++){
		for(d=0; d<=nAlpha(A); d++){
			score += M->freq[d]*M->score[j][d];
		}
	}
	return score;
}

void	PutCircularSMatrix(FILE *fptr, smx_typ M)
{
	Int4	j,d,v;
	a_type	A = M->A;

  if(M->posMatrix != NULL){
	fprintf(fptr,"\n       ");
	for(d=0; d<=nAlpha(A); d++){
		fprintf(fptr,"  %c ", AlphaChar(d,A));
	}
	fprintf(fptr,"  $ \n");
	for(j=0; j <= M->posMtrxLen+1; j++){
		fprintf(fptr,"%c%4d: ",AlphaChar(M->qseq[j],A),j);
		for(d=0; d<=nAlpha(A); d++){
			v = M->posMatrix[j][d];
			if(v==SHRT_MIN) fprintf(fptr,"inf ");
			else fprintf(fptr,"%3d ", v);
		}
		v = M->posMatrix[j][d];
		if(v==SHRT_MIN) fprintf(fptr,"inf ");
		else fprintf(fptr,"%3d ", v);
		fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n");
   }
}

void	PutSMatrix(FILE *fptr, smx_typ M)
{
	Int4	j,d;
	a_type	A = M->A;

	fprintf(fptr,"\n    ");
	for(d=0; d<=nAlpha(A); d++){
		fprintf(fptr,"   %c ", AlphaChar(d,A));
	}
	fprintf(fptr,"\n");
	for(j=1; j <= M->K; j++){
		fprintf(fptr,"%2d: ",j);
		for(d=0; d<=nAlpha(A); d++){
			fprintf(fptr,"%4d ", M->score[j][d]);
		}
		fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n");
#if 0	// put pseudo informtion content of smatrix...
	double	p,q,tot_info,total,info,sum;
	h_type 	H;

	H = Histogram("smatrix information",0,M->K+1,1.0);
#if 0  // precompute total info
        for(tot_info=0.0, j=1; j<= M->K; j++){
            for(total=0.0, d=1; d <= nAlpha(M->A); d++){
              total += exp((double)M->score[j][d]);
	    }
            for(info=0.0,d=1; d <= nAlpha(M->A); d++){
                p = exp((double)M->score[j][d])/total;
                if(p > 0.0){ q = M->freq[d]; info += p*log(p/q); }
		else print_error(" PutSMatrix( ): this should not happen!");
            }
            tot_info += info;
	}
#endif
// WARNING: NEED TO ADD PERNATS PARAMETER.
	fprintf(fptr,"Information (relative entropy) contribution in 1/100th nats\n");
	fprintf(fptr,"POS  ");
        for(sum=0.0,j=1; j<= M->K; j++){
            for(total=0.0, d=1; d <= nAlpha(M->A); d++){
              total += exp((double)M->score[j][d]);
              if(j==1) fprintf(fptr,"%3c", AlphaChar(d, M->A));
            }
            if(j==1) fprintf(fptr,"  Info (%% tail)\n");
            fprintf(fptr,"%4d ",j);
            for(info=0.0,d=1; d <= nAlpha(M->A); d++){
                p = exp((double)M->score[j][d])/total;
                if(p > 0.0){ q = M->freq[d]; info += p*log(p/q); }
		else print_error(" PutSMatrix( ): this should not happen!");
                fprintf(fptr,"%3d", (Int4)(100*p+0.5));
            }
	    sum+=info; 
	    if(tot_info >= 2.0*sum) p = (100.*sum/tot_info);
	    else p = (100.*((tot_info - sum) + info)/tot_info);
            fprintf(fptr,"   %3d (%d)\n",(Int4)(100*info+0.5),(Int4)(p+0.5));
	    IncdMHist(j, (Int4)(100*info+0.5), H);
        }
	fprintf(fptr,"     %61s %4d (%d ave)\n", "total:",
		(Int4)(100*tot_info+0.5), (Int4)(100*(tot_info/(double)M->K)+0.5)); 
	PutHist(fptr,60,H); NilHist(H);
#endif
}

Int4	SplitScoreSMatrix(unsigned char *seq, Int4 n, Int4 *start, Int4 *leng, 
	smx_typ M)
/* determine score by split regions. WARNING: assumes lengths sum to M->K */
{
	register Int4	m,j,c,s,score;

	for(score=0,m=1,c=1; m <= n; m++){
	    for(s=start[m],j=0; j < leng[m]; s++,j++,c++){
		score += M->score[c][seq[s]];
	    }
	}
	return score;
}

Int4    LocalScoreSMatrix2(register unsigned char *seq, register Int4 start, 
	register smx_typ M)
/*** Get the optimum local score of seq at start for matrix M ****/
{
        register Int4   j,score,max;

        for(seq += start-1, score=max=0,j = 1; j <= M->K; j++){
                score += M->score[j][seq[j]];
		if(score < 0) score = 0;
		else { max = MAXIMUM(Int4,score, max); }
        }
        return max;
}

Int4    LocalScoreSMatrix(register unsigned char *seq, register Int4 start, 
	register smx_typ M)
{ return local_score_smatrix(((seq)+(start)-1),M->K,M->score); }

Int4	local_score_smatrix(register unsigned char *seq, register Int4 j, 
	register Int4 **scr) /*** inner loop function ***/
{
	register Int4	score=0,max=0;

	while(j){ 
		score += scr[j][seq[j]]; 
		if(score < 0) score = 0;
		else if(score > max) max = score;
		j--; 
	} 
	return max;
}

Int4	OptInsertScoreSMatrix(unsigned char *seq, Int4 start, Int4 overlap, 
	smx_typ M1, smx_typ M2)
// find the optimum imsertion point 
/************************************************************************
   Example (overlap = 3):

   --[_M1__|//]----------          --[_M1__|//]---------------  return:
                                            abc---                 0
            abc             ----->          ab---c                 1
                                            a---bc                 2
                                            ---abc                 3
   --------[//|__M2__]---          ------------[//|__M2__]-----

   abc---,  ab---c, a---bc, ---abc      (find the best configuration)
   i=M1->K...             ...M1->K - overlap
   j=overlap...           ...0

 ************************************************************************/
{
	Int4		s,s1,s2,i,j,k,max,s_max,offset,end;
	unsigned char	*sq;
	Int4		**scr,**scr0;

	if(overlap > M1->K || overlap > M1->K || overlap < 1)
		print_error("OptInsertScoreSMatrix( ) input error!");
	seq = seq + (start + M1->K - overlap - 1); 
	end = 2*overlap;
	sq = new unsigned char[end+2];
	scr = new Int4*[end+2];
	scr0 = M1->score + (M1->K - overlap - 1);
	for(j=1; j<=overlap; j++){ sq[j]=seq[j]; scr[j]=scr0[j]; }
	for(i=1; j <= end; i++,j++){ sq[j]=0; scr[j]=M2->score[i]; }

	Int4	t;
	for(s_max=-9999,offset=overlap,t=end; offset >= 0; offset--,t--){
	   for(s=0,j=1; j <= end; j++){ s += scr[j][sq[j]]; }
	   if(s > s_max) { s_max=s; max=offset; }
	   for(j=1; j <= end; j++){
		i=AlphaChar(sq[j],M1->A);
		fprintf(stderr,"%c",i);
	   }
	   if(offset > 0){ sq[t]=sq[offset]; sq[offset]=0; }
	   fprintf(stderr,"\noffset = %d; s = %d; s_max = %d\n",offset,s,s_max);
	}
	delete []sq;
	delete []scr;
	return overlap - max;
}

Int4	SubScoreSMatrix(unsigned char *seq, Int4 start, Int4 smx_start,
	Int4 smx_end, smx_typ M)
// compute a subscore for SMatrix from smx_start (1..LenSMX) to 
// smx_end (1..LenSMX) && smx_start <= smx_end.
{ 
	return score_smatrix((seq+start-1)+smx_start,
		(smx_end - smx_start), M->score + smx_start); 
}


Int4	ScoreSMatrix(register unsigned char *seq, Int4 start, smx_typ M)
{ return score_smatrix(((seq)+(start)-1),M->K,M->score); }

Int4	score_smatrix(register unsigned char *seq, register Int4 j, 
	register Int4 **scr) /*** inner loop function ***/
{
	register Int4	score=0;

	while(j){ score += scr[j][seq[j]]; j--; } 
	return score;
}

double	SMatrixProbFast(Int4 score, smx_typ M)
{
	double	f0[3];

	if(M->changed){M->changed=FALSE;M->calc_prob=M->calc_stats=TRUE;} 
	if(M->calc_stats) stats_smatrix(M);
	if(score <= M->cmin[M->K]) return 1.0;
	if(M->calc_prob){
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
#if 0
		fprintf(stderr,
			"score (%d) = %g s.d.;",
			  score, n=((double)score - M->mean)/M->sd);
		if(score >= M->cmin[M->K]) {
			fprintf(stderr," p=%g",M->f[score]);
			fprintf(stderr," pnorm=%g\n",0.5*gammq(0.5,(n*n*0.5)));
		} else  fprintf(stderr," p=??\n");
#endif
		M->calc_prob = FALSE; 
	}
	if(score < M->cmin[M->K]) return 1.0;
	else if(score > M->cmax[M->K]) return 0.0;
	else return M->f[score];
}

double	SMatrixProb(Int4 score, smx_typ M)
{
	double	f0[3];

	if(M->changed){
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = M->calc_prob = FALSE; 
	}
	if(score < M->cmin[M->K]) return 1.0;
	else if(score > M->cmax[M->K]) return 0.0;
	else return M->f[score];
}

Int4	MinSegSMatrix(unsigned char *minseg, smx_typ M)
{
	Int4	s,i,r,K=M->K;
	double	f0[3];

	if(M->changed){	
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = FALSE; M->calc_prob = FALSE; 
	}
	for(i=1; i<=K; i++){
	   for(r = 0; r <= M->nlet; r++){
	    	if(M->score[i][r] == M->min[i]){
			minseg[i]=r; break;
		}
	   }
	} return i;
}

Int4	MaxSegSMatrix(unsigned char *maxseg, smx_typ M)
{
	Int4	s,i,r,K=M->K;
	double	f0[3];

	if(M->changed){	
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = FALSE; M->calc_prob = FALSE; 
	}
	for(i=1; i<=K; i++){
	   for(r = 0; r <= M->nlet; r++){
		s = M->score[i][r];
	    	if(s == M->max[i]){
			maxseg[i]=r; break;
			// maxseg[i-1]=AlphaChar(r,M->A); break; 
		}
	   }
	}
	// maxseg[i-1]='\0'; 
	return i;
}

Int4	MaxScoreSMatrix(smx_typ M)
{
	double	f0[3];

	if(M->changed){
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = M->calc_prob = FALSE; 
	}
	return M->cmax[M->K];
}

Int4	MinScoreSMatrix(smx_typ M)
{
	double	f0[3];

	if(M->changed){
		stats_smatrix(M);
		if(M->f0!= NULL){ free(M->f0); M->f0=NULL; }
		f0[0] = 1.0; min_max_smatrix(M);
		M->fx = NULL;
		smatrix_prob(0, f0, M);
		M->changed = M->calc_prob = FALSE; 
	}
	return M->cmin[M->K];
}

/****************************** private **********************************/
void	stats_smatrix(smx_typ M)
{
	Int4	i,r;
	double	E,mean,variance,V;
	
	for(mean=variance=0.0,i=1; i <= M->K; i++){
	    for(E=0.0,r=0; r <= M->nlet; r++){
		E += (double) M->score[i][r] * M->freq[r];
	    }
	    mean += E;
	    for(V=0.0,r=0; r <= M->nlet; r++){
		V += pow((E - (double) M->score[i][r]),2.0);
	    }
	    variance += V/(M->nlet - 1);
	}
	M->mean = mean; M->var = variance; M->sd = sqrt(variance);
	M->calc_stats = FALSE;
}

Int4	min_max_smatrix(smx_typ M)
/********************************************************************
      cmin0[K]                       n sd              cmax[K] 
	|-----------------------------|-----------------|

	cmin0[K-1]                cmin[K-1]        cmax[K-1]
	   |------------------------|-----------------|
		:		:		:
	        cmin0[1] cmin[1]     cmax[1]
	             |----|------------|

	                cmin[0]=cmax[0]=0
			       |

	cmax[0]   = 0
	cmax[i]   = cmax[i-1] + max[i];

	cmin[K]   = mean + n*sd;  where n = nsd
	cmin[i-1] = cmin[i] - max[i];
**********************************************************************/
{
	Int4	s,i,r,min,K=M->K;
	double	minscore;

	if(M->calc_stats) stats_smatrix(M);
	M->cmax[0] = 0;
	for(M->max[0]=M->min[0]=0,i=1; i<=K; i++){
	   M->max[i]=INT4_MIN; M->min[i]=INT4_MAX;
	   for(r = 0; r <= M->nlet; r++){
		s = M->score[i][r];
	   	M->max[i]=MAXIMUM(Int4,M->max[i],s);
	   	M->min[i]=MINIMUM(Int4,M->min[i],s);
	   }
	   M->cmax[i] = M->cmax[i-1] + M->max[i];
	}
	minscore = (M->mean + M->nsd*M->sd);
	if(minscore  >  -99.0){
		M->cmin[K] = (Int4)minscore - 1;
	} else M->cmin[K] = 0;	/** beware of negative infinity socres **/
	// for short matrices cmin[K] = mean + 2 std_dev may be > cmax[K]!!!
	if(M->cmin[K] > M->cmax[K]) M->cmin[K] = MINIMUM(Int4,0,M->cmax[K]);
	for(i=K-1; i>=1; i--){ M->cmin[i] = M->cmin[i+1] - M->max[i+1]; }
	M->cmin[0] = 0;
	for(min = 0,i=1; i<=K; i++){
	    min += M->min[i]; M->cmin[i] = MAXIMUM(Int4,M->cmin[i],min);
	    // assert(M->cmin[i] <= M->cmax[i]);
	}
#if 0
	// fprintf(stderr,"cmax=%d; cmin=%d; mean = %f\n",M->cmax[K],M->cmin[K],M->mean);
	// M->cmin[i] = min; // DEBUG.
#endif
	return M->cmax[K];
}

void	smatrix_prob(register Int4 k, register double *f, smx_typ M)
/*************************************************************************

  Call: smatrix_prob(0, f[0]=1.0, M)

 *************************************************************************/
{
	register double	*f2,*freq = M->freq;
	register Int4	s,s2,r,min,**S=M->score;
	double 	*f0;

	if(k == M->K){ 
		s=M->cmax[k];
		f[s+1]=0.0;
		for( ; s>=M->cmin[k]; s--){
			f[s] = f[s] + f[s+1];
#if 0	// WARNING: Insure is complaining here about Read and Write bad indices.
	// This is occuring when cmin is > 0; in which case the beginning of the
	// array points to unallocated memory (see below). This is not a problem, 
	// however, as the indices s and s+1 are always set to point to an
	// allocated position in the array.  
			fprintf(stderr,
			  "s=%d; k = %d; cmin = %d; cmax = %d; f[s]=%g; f[s+1]=%g\n",
			   s,k,M->cmin[k],M->cmax[k],f[s],f[s+1]);
#endif
		} // PutSMatrix(stderr,M);
		// assert(k!=7);
	} else {
	    k++;
#if 0
if(M->cmax[k] < M->cmin[k]) {
	fprintf(stderr,"k = %d; cmin = %d; cmax = %d\n",k,M->cmin[k],M->cmax[k]);
	PutSMatrix(stderr,M);
}
#endif
	    NEW(f0,(M->cmax[k] - M->cmin[k]) + 4,double); 
	    f2 = f0 - M->cmin[k]; min = M->cmin[k];
	    for(s=M->cmin[k-1]; s <= M->cmax[k-1]; s++){
		if(f[s] > 0.0){
		    for(r = M->nlet; r >= 0; r--){
			s2 = s + S[k][r];
			if(s2 >= min) f2[s2] += freq[r]*f[s];
#if 0
if(s2 > M->cmax[k])
	fprintf(stderr,"s2 = %d > max: s = %d; S[k][r] =%d; r=%d\n",
		s2,s,S[k][r],r);
#endif
		    }
		}
	    }
	    if(M->fx != NULL) free(M->fx);
	    M->fx = f0;
	    smatrix_prob(k, f2, M);
	    if(k==M->K){M->f=f2; M->f0=f0; }
	} 
}

#if 0
e_type ConsensusSMatrix(smx_typ M)
// From libalex...
{
        e_type          E;
        unsigned        char *cons;
        Int4            i,k,s,p=1,max_k,len;
        double          max_val;
        a_type 		A = SMatrixA(M);

        len = LenSMatrix(M);
        NEW(cons,len+3,unsigned char);
	for(s=1;s<=LenSMatrix(M);s++){
	  max_k=1; max_val = ValSMatrix(s,1,M);
	  for(k=2; k<=nAlpha(A); k++){
		if(ValSMatrix(s,k,M) > max_val) {max_val=ValSMatrix(s,k,M); max_k=k;}
	  } cons[p++] = max_k;
	}   
        E = MkSeq("cons", len, cons); free(cons);
        return E;
}
#endif

e_type ConsensusSMatrix(smx_typ *M, Int4 nbr, Int4 rpts)
// From libalex...
{
        e_type          E;
        unsigned        char *cons;
        Int4            i,k,s,j,n,p=1,max_k,totlen=0;
        double          max_val;
        a_type 		A = SMatrixA(M[1]);

        for(i=1;i<=nbr;i++) totlen += LenSMatrix(M[i]);
        NEW(cons,rpts*totlen+3,unsigned char);
        for(j=1;j<=nbr;j++){
             for(s=1;s<=LenSMatrix(M[j]);s++){
                 max_k=1; max_val = ValSMatrix(s,1,M[j]);
                 for(k=2; k<=nAlpha(A); k++){
                    if(ValSMatrix(s,k,M[j]) > max_val){
			max_val=ValSMatrix(s,k,M[j]); max_k=k;
		    }
                 }
                 for(n=0;n<rpts;n++){ cons[n*totlen+p] = max_k; }
                 p++;
             }   
        } E = MkSeq("cons", rpts*totlen, cons); free(cons);
        return E;
}

void    smatrix_error(char *s){fprintf(stderr,"SMatrix: %s\n",s);exit(1);}

