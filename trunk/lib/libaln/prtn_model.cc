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

#include "prtn_model.h"

ptm_typ *MakePrtnModels(char *snfile, a_type A, Int4 maxrpts, char method,
	char mode, Int4 maxLength, BooLean weight,float minmap,
	double pseudo, double minprob, double *freq)
{
	ptm_typ	*PM;
	char	str[1000],*s2;
	FILE	*fp,*fp2;
	Int4	i,n;

	print_error("MakePrtnModels() needs further testing: Comment this print_error out");
	fp=open_file(snfile,"","r");
	for(n=0; (s2=fgets(str,999,fp)) != NULL; ){
	   if(strncmp(s2,"//",2)==0) n++;
	} fclose(fp);	
	if(n==0) return NULL;

	NEW(PM,n+3,ptm_typ);
	fp=open_file(snfile,"","r");

	fp2=tmpfile();
	s2=fgets(str,999,fp); fprintf(fp2,"%s",s2);
	for(i=1; (s2=fgets(str,999,fp)) != NULL; ){
	   if(strncmp(s2,"//",2)==0){
		rewind(fp2);
		PM[i]=MakePrtnModel(fp2,A,maxrpts,method,mode,maxLength,
			weight,minmap,pseudo,minprob,freq);
		fclose(fp2); fp2=tmpfile(); i++;
	   }
	   fprintf(fp2,"%s",s2);
	} 
	rewind(fp2); fclose(fp);
	PM[i]=MakePrtnModel(fp2,A,maxrpts,method,mode,maxLength,
			weight,minmap,pseudo,minprob,freq);
	fclose(fp2); 
	return PM;
}

ptm_typ MakePrtnModel(char *snfile, a_type A, Int4 maxrpts, char method,
	char mode, Int4 maxLength, BooLean weight,float minmap,
	double pseudo, double minprob, double *freq)
{
        FILE *fptr=open_file(snfile,"","r");
	ptm_typ PM=MakePrtnModel(fptr,A,maxrpts,method,mode,maxLength,weight,minmap,
			pseudo,minprob,freq);
	fclose(fptr);
	return PM;
}

ptm_typ MakePrtnModel(FILE *fptr, a_type A, Int4 maxrpts, char method,
	char mode, Int4 maxLength, BooLean weight,float minmap,
	double pseudo, double minprob, double *freq)
{
	ptm_typ F;
	Int4	end,m,c,len,totcol,i,j,*map[2],total;
	double	info;
	BooLean **use;
	Int4	*numcol;

dh_type H;

	NEW(F,1,prtn_model_type);

	if(islower(mode)) F->gapfunct= TRUE;
	else F->gapfunct= FALSE;
	F->segmask=TRUE;
	F->mode=mode;
	F->minmap=minmap;
	F->max_gap_score=0;
	F->maxLength=maxLength;
	F->weights=weight;
	F->observedGap=NULL;
	F->totGaps=0;
	F->A = A; 
	F->N = 0; 
	NEW(F->temp,nAlpha(A)+2, double); 
	F->method = method; 
	F->maxrpts= MAXIMUM(Int4,1,maxrpts);
	BooLean okay=ReadPrtnModelSMA(fptr,pseudo,freq,minprob,F);
// assert(F != 0);
// fprintf(stderr,"DEBUG 2: okay = %d\n",okay);
	if(!okay){ NilPrtnModel(F); return NULL; }
	// if(!okay) print_error("MakePrtnModel( ) input error");
/**************************************************************************
  try using 1/4rd most informative columns using a dheap  
 **************************************************************************/
	for(totcol=0,m=1; m<=F->N;m++) totcol+=LenWModel(F->M[m]);
	H = dheap(totcol+1,4);
	NEW(map[0],totcol+2 ,Int4);
	NEW(map[1],totcol+2 ,Int4);
	NEWP(use,F->N+2 ,BooLean);
	NEW(numcol,F->N+2 ,Int4);
	for(m=1; m<=F->N;m++){
	    len = LenWModel(F->M[m]);
	    info = ExpectedScoreWModel(NULL, F->M[m]);
	    insrtHeap(m,(keytyp)-info,H);
	}
	i = F->N/2;
	while(ItemsInHeap(H) > i){
	    m = delminHeap(H);
	    len = LenWModel(F->M[m]);
	    NEW(use[m],len+2,BooLean);
	}
	while(ItemsInHeap(H) > 0) delminHeap(H);
	for(i=1,m=1; m<=F->N;m++){
	   if(use[m] != NULL){
	    len = LenWModel(F->M[m]);
	    for(c=1; c<=len;c++){
		info = InfoColWModel(c, F->M[m]);
		insrtHeap(i,(keytyp)-info,H);
		map[0][i] = m; map[1][i] = c;
		i++;
	    }
	  }
	}
	total = i-1; end = 1+(total/3);
	for(j=1; j <= end && (i=delminHeap(H)) !=NULL; j++){
		m=map[0][i]; c=map[1][i];
		use[m][c]=TRUE; numcol[m]++;
	}
	for(m=1; m<=F->N;m++){
	  if(use[m]!=NULL){
	    len = LenWModel(F->M[m]);
#if  0
	    for(j=0,c=1; c<=len;c++){
   		fprintf(stderr,"info[%d][%d] = %.2f",m,c,InfoColWModel(c, F->M[m]));
   		if(use[m][c]) fprintf(stderr,"*\n"); else fprintf(stderr,"\n");
	    }
#endif
	  } 
#if  0
	  else fprintf(stderr,"\n"); 
#endif
	}
	for(i=1,m=1; m<=F->N;m++){
	   if(use[m]!= NULL){
#if  0
	     fprintf(stderr,
	    "Model %d (%d columns): expected full score = %.2f (heuristic = %.2f)\n\n",
		m,numcol[m],ExpectedScoreWModel(NULL, F->M[m]),
		ExpectedScoreWModel(use[m], F->M[m]));
#endif
	     free(use[m]);
	   } else {
#if  0
	      fprintf(stderr,
		"Model %d (%d columns): expected full score = %.2f (heuristic = 0.0)\n\n",
		  m,numcol[m],ExpectedScoreWModel(NULL, F->M[m]));
#endif
	   }
	}
	free(use); free(map[0]); free(map[1]);
	Nildheap(H);
	if(F->overlap != NULL){
		MEW(F->X[0],maxLength+F->overlap[0]+3,Int4);
		MEW(F->X[1],maxLength+F->overlap[0]+3,Int4);
	} else {
		MEW(F->X[0],maxLength+3,Int4);
		MEW(F->X[1],maxLength+3,Int4);
	}
// fprintf(stderr,"maxLength = %d\n",maxLength);
	free(numcol);
	return F;
}

double  *ObservedPrtnModel(Int4 m,Int4 i,ptm_typ PM)
{
	if(m <= PM->N && m > 0){
		 return ObservedWModel((i), (PM)->M[(m)]);
	} else return NULL;
}

double  *ObservedFilterPrtnModel(Int4 m,Int4 i,ptm_typ PM)
{
	double	*observed;
	BooLean *conserved;
	Int4	j;

	if(m <= PM->N && m > 0){
		observed=ObservedWModel((i), (PM)->M[(m)]);
		if((conserved=ConservedSMA(m,i,PM->MA)) != NULL){
		   for(j=1; j<=nAlpha(PM->A); j++){
			if(conserved[j]) PM->temp[j] = observed[j];
			else PM->temp[j] = 0;
		   }
		}
		free(conserved);
		return PM->temp;
	} else return NULL;
}

void	NilPrtnModel(ptm_typ F)
{
	Int4	i;

        // if(F->Spg != NULL) { NilSpouge(F->Spg); free(F->scoreGap); }
	if(F->w) free(F->w); 
	if(F->sM) free(F->sM); 
	if(F->fsM) free(F->fsM); 
	if(F->smx) NilSMatrix(F->smx);
	if(F->M) for(i=1;i<=F->N; i++) { 
		NilWModel(F->M[i]); 
		if(F->M2[i] != NULL) NilWModel(F->M2[i]); 
	}
	if(F->overlap != NULL) free(F->overlap);
	if(F->X[0]) free(F->X[0]); 
	if(F->X[1]) free(F->X[1]); 
	if(F->M2) free(F->M2); 
	if(F->M) free(F->M); 
	if(F->freq) free(F->freq); 
	if(F->temp) free(F->temp);
	if(F->observedGap != NULL){
		for(i=0;i<=F->N; i++) free(F->observedGap[i]);
		free(F->observedGap);
	} 
	if(F->MA) NilSMA(F->MA);
	free(F);
}

float	**InfoPrtnModel(Int4 purge_cutoff, ptm_typ PM)
{
	float   **info = InfoSMA(purge_cutoff, PM->MA);
	return edit_info_protein_model(info, PM);
}

float	**ExcessInfoPrtnModel(char *string, Int4 purge_cutoff, ptm_typ PM)
{
	float   **info = ExcessInfoSMA(string, purge_cutoff, PM->MA);
	return edit_info_protein_model(info, PM);
}

float	**edit_info_protein_model(float **info0, ptm_typ PM)
{
	float   **info;
	Int4	n,m;

	NEWP(info,NumModelsPrtnModel(PM) +2,float);
	for(n=0,m=1; m<= ntypSMA(PM->MA); m++){
		if(fieldmapSMA(m,PM->MA) >= PM->minmap){
		    n++; info[n] = info0[m];
		} else free(info0[m]);
	}
	free(info0);
	if( n != NumModelsPrtnModel(PM)) print_error("This should not happen!");
	return info;
}

BooLean	ReadPrtnModelSMA(FILE *fptr, double pseudo,double *freq, double minprob,
	ptm_typ F)
{ 
    Int4		m,n,s0,s,t,total,*t_used,i,length;
    char		*null0=NULL,*null=NULL;
    unsigned char	*seq=NULL;
    BooLean		ignore,use_null=TRUE;
    double		*wt;
    a_type		A = F->A;
    sma_typ		MA;
    
    MA = ReadSMA(fptr);
    if(MA == NULL) return FALSE;   
    if(isupper(F->method)){ use_null=FALSE; F->method=tolower(F->method); }
    NEW(F->freq,nAlpha(A)+2,double);
    for(i=0; i<=nAlpha(A); i++) F->freq[i]= freq[i];
    
    if(F->weights) wt = WeightsSMA(MA); 
    else {
	NEW(wt,nseqSMA(MA)+2,double);
	for(n=1; n<= nseqSMA(MA); n++) wt[n]=1.0;
    }
    NEW(F->M,ntypSMA(MA)+2,wm_type);
    NEW(F->M2,ntypSMA(MA)+2,wm_type);
    NEW(t_used,ntypSMA(MA)+3,Int4);
    for(F->N=m=0,t=1; t<=ntypSMA(MA); t++){
        length = lengthSMA(t,MA);
	if(null != NULL) free(null); NEW(null, length+3, char);
	null0 = nullSMA(t,MA);
	ignore = FALSE;
	for(s0=0,s=1; null0[s] != 0; s0++,s++){
	    switch(null0[s0]){
	      case '.': if(use_null) null[s]='.'; else null[s]='*'; break;
	      case '^': null[s]='^'; ignore=TRUE; break;
	      case '!': null[s]='!'; break;
	      case '*': null[s]='*'; break;
	      default: print_error("ReadPrtnModelCMSA( ) error 1"); break;
	    }
        }
	if(fieldmapSMA(t,MA) >= F->minmap){
	  F->N++; m++; t_used[m]=t; 
	  F->M[m]=MkWModel(null,length,pseudo,F->freq,A);
	  SetMethodWModel(F->method,F->M[m]);
	  if(ignore) {
	    F->M2[m]=MkWModel(NULL,length,pseudo,F->freq,A);
	    SetMethodWModel(F->method,F->M2[m]);
	  } else F->M2[m] = NULL;
          for(n=1; n<= nseqSMA(MA); n++){
	    seq=seqSMA(t,n,MA);
	    if(probSMA(n,t,MA) >= minprob){
		Add2WModel(seq,1,wt[n],F->M[m]);
	    	if(ignore) Add2WModel(seq,1,wt[n],F->M2[m]);
	    }
	  }
        } 
    }
    free(wt);
    if(null != NULL) free(null);
    AllocatePrtnModel(t_used,MA,F->maxLength, F);
    free(t_used); F->MA = MA;
    return TRUE;
}

BooLean	CheckPrtnModel(Int4 *Score, e_type E,UInt8 total, double maxEval, ptm_typ F)
/**** DO AN INITIAL FAST CHECK ****/ 
{
	Int4	len,n,m,score;
	double	factor,target,pval;

	len = LenSeq(E);
	n = len - F->tot_lenOL; m = F->N;
	factor = bico(n+m,m)*(double)total/(double) len; /* Altschul's way */
	target = maxEval/((double)total/(double) len); 
	switch (F->mode) {
	     case 'D': case 'd':
#if 0
	   	score=LocalAlnSeqSMatrix(len,SeqPtr(E),F->N,F->sM,NULL);
#else
	   	score=FastLocalAlnSeqSMatrix(len,SeqPtr(E),F->N,F->sM,F->X);
#endif
		break;
	     case 'G': case 'g':
#if 1
	     case 'a':
#endif
	   	score=FastAlnSeqSMatrix(len,SeqPtr(E),F->N,F->sM,F->X);
		break;
	     case 'O': case 'o': case 'e':
	   	score=FastAlnSeqOverlapSMatrix(len,SeqPtr(E),F->N,
			F->overlap,F->sM,F->X);
		break;
	     default: print_error("PrtnModelScan( ) mode input error");
	} 
	*Score = score;
	if(F->mode == 'G' || F->mode == 'O'){
	   pval = SMatrixProb(score,F->smx)*factor; 
	   if(pval > maxEval) return FALSE;
	   else return TRUE;
	} else { // NOTE: totalOL = 0 if overlap == 0;
	    print_error("Spouge gap functions deactivated");
#if 0
	    if((score + F->max_gap_score) < 
			SpougeThreshold(F->Spg,target,len+F->totalOL)){
				return FALSE; 
	    } else return TRUE; 
#endif
	} 
}

void	PutSmxPrtnModel(Int4 t, ptm_typ PM)
{
	smx_typ smx;

	if(t > 0 && t <= PM->N){
		smx = GetSMatrixWModel(PM->M[t]);
		PutSMatrix(stdout, smx);
		// if(PM->Spg != NULL) PutHistSpouge(stdout,t,1,PM->Spg);
	}
}

wm_type	MergePrtnModels(e_type E, ptm_typ PM)
// Create and return a model spanning sequence E (use Null sites between spaces).
{
        Int4    score,totN,*p,len;
        char    mode = 'G';
        float   *pv;
        double  pvalue;
	sni_typ sI;
	wm_type	M;

        totN = NmaxPrtnModel(PM);
        NEW(p,totN+3,Int4); NEW(pv,totN+3,float);
	len=LenSeq(E);
	if(len < TotLenPrtnModel(PM)) 
		print_error("MergePrtnModels( ): seqlen < model leng");
	switch (PM->mode) {
	   case 'G': AlnSeqSMatrix(len,SeqPtr(E),PM->N,PM->sM,p); 
		break;
	   case 'g': print_error("Spouge gap functions deactivated");
		// GapFuncAlnSeqSMatrix(len, SeqPtr(E), PM->N, PM->sM,p,GapScoreSpouge(PM->Spg));
		break;
	   default: print_error("MergePrtnModels( ) mode input error");
	}
#if 0
        double  maxEval=1.;
        if(ComparePrtnModel(E,LenSeq(E),maxEval,p,pv,&pvalue,PM) != 1){
		free(p); free(pv); return NULL;
	}
        sI=MakeScanInfo(E,pvalue,NumModelsPrtnModel(PM),1,pv,p,1,mode);
        PutScanInfo(stderr, sI, PrtnModelLengths(PM), FALSE, PM->A);
        fprintf(stderr,"pvalue = %.2f\n",pvalue);
        NilScanInfoRtnE(sI); 
#endif

	M=MergeWModels(NumModelsPrtnModel(PM),p,LenSeq(E),PM->freq, PM->M);
	free(pv); free(p);
	if(M == NULL) print_error("input error in MergePrtnModels( )");
	return M;
}

Int4	ComparePrtnModel(e_type E,UInt8 total,double maxEval,Int4 *p,float *pv,
	double *pvalue,ptm_typ F)
{ return ComparePrtnModelRpts(E,total,maxEval,p,pv,pvalue,maxEval,F); }

double	EvaluePrtnModel(Int4 score, UInt4 len, ptm_typ PM)
// return evalue for a single sequence.
{
	double factor;
	Int4	n,m;
	if(len < PM->tot_lenOL) return 0;
	n = len - PM->tot_lenOL; 
	if(n < 0) n = 0; m = PM->N;
	factor = bico(n+m,m);
	return SMatrixProb(score,PM->smx)*factor;
}

Int4	ComparePrtnModelRpts(e_type E,UInt8 total,double maxEval,Int4 *p,float *pv,
	double *pvalue,double rptsEval,ptm_typ F)
/** return the number of significant repeats found; zero if none are found. **/
{
	Int4	m,m0,repeats=0,rmax,score,score2,totN,last_sum_score,r,sum_score;
	Int4	n,field,len = LenSeq(E);
	char	mode;
	a_type	A=F->A;
	double	factor,pval,pval2,mask_drop=10.0,prob;
	BooLean	masked;

	if(len < F->tot_lenOL) return 0;
	if(!CheckPrtnModel(&score,E,total,maxEval,F)) return 0;
	if(F->nonglobular) {
           masked= ProcessSeqPSeg(45,3.4,3.75,100,E,A);
	   if(F->segmask) { /* mask compositionally biased regions. */
		masked= (masked || ProcessSeqPSeg(14,2.2,2.5,100,E,A));
	   }
	} else if(F->segmask) { /* mask compositionally biased regions. */
	        masked=ProcessSeqPSeg(14,2.2,2.5,100,E,A); 
	} else masked=FALSE;
	n = len - F->tot_lenOL; m = F->N;
	factor = bico(n+m,m)*(double)total/(double) len; /* Altschul's way */
	switch (F->mode) {
	   case 'D':   
		print_error("Spouge gap functions deactivated");
#if 0
		score=LocalAlnSeqSMatrix(len,SeqPtr(E),F->N,F->sM,NULL);
		pval = SpougePvalue(F->Spg,total,len,score); 
		if(masked){
			score2=LocalAlnSeqSMatrix(len,XSeqPtr(E),F->N,F->sM,NULL);
		        pval2 = SpougePvalue(F->Spg,total,len,score2); 
		} else pval2=pval; 
#endif
		break;
	   case 'd': 
		print_error("Spouge gap functions deactivated");
#if 0
	        score=LocalGapFuncAlnSMatrix(len,SeqPtr(E),
				F->N,F->sM,NULL,GapScoreSpouge(F->Spg));
		pval = SpougePvalue(F->Spg,total,len+F->totalOL,score); 
		if(masked){
	                score2=LocalGapFuncAlnSMatrix(len,XSeqPtr(E),
				F->N, F->sM,NULL,GapScoreSpouge(F->Spg));
		        pval2 = SpougePvalue(F->Spg,total,len+F->totalOL,score2); 
		} else pval2=pval; 
#endif
		break;
	   case 'G':   /** score and positions already computed **/
#if 0	// for heuristic need to recompute 
fprintf(stderr,"heuristic score = %d --> ",score);
	   	score=AlnSeqSMatrix(len,SeqPtr(E),F->N,F->sM,NULL);
fprintf(stderr,"; true score = %d\n",score);
#endif
		pval = SMatrixProb(score, F->smx)*factor; 
		if(masked){
		       score2=AlnSeqSMatrix(len,XSeqPtr(E),F->N,F->sM,NULL); 
		       pval2 = SMatrixProb(score2, F->smx)*factor; 
		} else pval2=pval; 
		break;
	   case 'a':
	   case 'g':
		print_error("Spouge gap functions deactivated");
#if 0
		score=GapFuncAlnSeqSMatrix(len, SeqPtr(E),
				F->N, F->sM,NULL,GapScoreSpouge(F->Spg));
		pval = SpougePvalue(F->Spg,total,len,score); 
		if(masked){
	               score2=GapFuncAlnSeqSMatrix(len, XSeqPtr(E),
				F->N, F->sM,NULL,GapScoreSpouge(F->Spg));
		       pval2 = SpougePvalue(F->Spg,total,len,score2); 
		} else pval2=pval; 
#endif
		break;
	   case 'O':
		pval = SMatrixProb(score, F->smx)*factor; 
		if(masked){
		       score2=AlnSeqOverlapSMatrix(len,XSeqPtr(E),
				F->N,F->overlap,F->sM,NULL); 
		       pval2 = SMatrixProb(score2, F->smx)*factor; 
		} else pval2=pval; 
		break;
	   case 'o':
		print_error("Spouge gap functions deactivated");
#if 0
	        score = GapFuncAlnSeqOverlapSMatrix(len,SeqPtr(E),
			F->N,F->sM,NULL,GapScoreSpouge(F->Spg),F->overlap);
		pval = SpougePvalue(F->Spg,total,len+F->totalOL,score); 
		if(masked){
	               score2 = GapFuncAlnSeqOverlapSMatrix(len,XSeqPtr(E),
			 F->N,F->sM,NULL,GapScoreSpouge(F->Spg),F->overlap);
		       pval2 = SpougePvalue(F->Spg,total,len+F->totalOL,score2); 
		} else pval2=pval; 
#endif
		break;
	   case 'e':
		print_error("Spouge gap functions deactivated");
#if 0
	        score = GapFuncAlnSeqOverlapSMatrixCterm(len,SeqPtr(E),
			F->N,F->sM,NULL,GapScoreSpouge(F->Spg),F->overlap);
		pval = SpougePvalue(F->Spg,total,len+F->totalOL,score); 
		if(masked){
	               score2 = GapFuncAlnSeqOverlapSMatrixCterm(len,XSeqPtr(E),
			 F->N,F->sM,NULL,GapScoreSpouge(F->Spg),F->overlap);
		       pval2 = SpougePvalue(F->Spg,total,len+F->totalOL,score2); 
		} else pval2=pval; 
#endif
		break;
	   default: print_error("PrtnModelScan( ) mode input error");
	}
	// routine to capture 'G' & 'O' hits within 'g' & 'o' mode
	mode = F->mode;
	if(pval <= maxEval && (pval2 <= (maxEval*mask_drop))) repeats = 1;
	else if(F->mode == 'g' || F->mode == 'a'){
	   mode = 'G';
   	   score=FastAlnSeqSMatrix(len,SeqPtr(E),F->N,F->sM,F->X);
	   pval = SMatrixProb(score, F->smx)*factor; 
	   if(masked){
	       score2=AlnSeqSMatrix(len,XSeqPtr(E),F->N,F->sM,NULL); 
	       pval2 = SMatrixProb(score2, F->smx)*factor; 
	   } else pval2=pval; 
	   pval*=10.; pval2*=10.; // Phil Green type of adjustment
	   // should really multiple 'g' pval times 10/9 = 1.1111
	   // but this doesn't matter much (cutoff p <= 0.0100 -> p<=0.0111)
	   if(pval <= maxEval && (pval2 <= (maxEval*mask_drop))) repeats=1;
	   else repeats=0;
	} else if(F->mode == 'o'){
	   mode = 'O';
   	   score=FastAlnSeqOverlapSMatrix(len,SeqPtr(E),F->N,
		F->overlap,F->sM,F->X);
	   if(masked){
	       score2=FastAlnSeqOverlapSMatrix(len,XSeqPtr(E),
			F->N,F->overlap,F->sM,F->X); 
	       pval2 = SMatrixProb(score2, F->smx)*factor; 
	   } else pval2=pval; 
	   pval*=10.; pval2*=10.; // Phil Green type of adjustment
	   if(pval <= maxEval && (pval2 <= (maxEval*mask_drop))) repeats=1; 
	   else repeats=0;
	} else if(F->mode == 'e'){
	   mode = 'O';
   	   score=FastAlnSeqOverlapSMatrix(len,SeqPtr(E),F->N,
		F->overlap,F->sM,F->X);
	   if(masked){
	       score2=FastAlnSeqOverlapSMatrix(len,XSeqPtr(E),
			F->N,F->overlap,F->sM,F->X); 
	       pval2 = SMatrixProb(score2, F->smx)*factor; 
	   } else pval2=pval; 
	   pval*=10.; pval2*=10.; // Phil Green type of adjustment
	   if(pval <= maxEval && (pval2 <= (maxEval*mask_drop))) repeats=1; 
	   else repeats=0;
	} else repeats = 0;
	if(repeats != 0){	// replaces next if(pval...)
           /*** If there is a match then check for repeats ***/
	   *pvalue = log10(pval); repeats=1;
           if(F->maxrpts > 1){
                 rmax = MINIMUM(Int4,len/F->tot_lenOL,F->maxrpts);
                 factor = bico((len + F->N - F->tot_lenOL),F->N);
		 // ^adjusted for single sequence, not entire database. 
                 for(last_sum_score = score, r=2; r <= rmax; r++){
                    totN = r*F->N;
	            switch (mode) {
	        	case 'D': 
			  print_error("Spouge gap functions deactivated");
#if 0
	   	    	  sum_score=LocalAlnSeqSMatrix(len,SeqPtr(E),totN,F->sM,NULL);
		          prob=SpougePvalue(F->Spg,total,len,sum_score-last_sum_score); 
#endif
			  break;
	        	case 'd': 
			  print_error("Spouge gap functions deactivated");
#if 0
			  sum_score=LocalGapFuncAlnSMatrix(len,SeqPtr(E),totN,F->sM,
						NULL,F->scoreGap);
		          prob=SpougePvalue(F->Spg,total,len,sum_score-last_sum_score); 
#endif
			  break;
	        	case 'G': 
                    	  sum_score=AlnSeqSMatrix(len,SeqPtr(E),totN,F->sM,NULL);
                          prob=factor*SMatrixProb(sum_score-last_sum_score,F->smx);
			  break;
	        	case 'a': 
	        	case 'g': 
			  print_error("Spouge gap functions deactivated");
#if 0
	            	  sum_score=GapFuncAlnSeqSMatrix(len,SeqPtr(E),totN,F->sM,
						NULL,F->scoreGap);
		          prob=SpougePvalue(F->Spg,total,len,sum_score-last_sum_score); 
#endif
			  break;
	   		case 'O': 
                    	  sum_score=AlnSeqOverlapSMatrix(len,SeqPtr(E),totN,
				F->overlap,F->sM,NULL);
                          prob=factor*SMatrixProb(sum_score-last_sum_score,F->smx);
			  break;
			case 'o':
			  print_error("Spouge gap functions deactivated");
#if 0
	            	  sum_score=GapFuncAlnSeqOverlapSMatrix(len,SeqPtr(E),
				totN,F->sM,NULL,F->scoreGap,F->overlap);
		          prob=SpougePvalue(F->Spg,total,len+F->totalOL,
				  sum_score-last_sum_score); 
#endif
			  break;
			case 'e':
			  print_error("Spouge gap functions deactivated");
#if 0
	            	  sum_score=GapFuncAlnSeqOverlapSMatrixCterm(len,SeqPtr(E),
				totN,F->sM,NULL,F->scoreGap,F->overlap);
		          prob=SpougePvalue(F->Spg,total,len+F->totalOL,
				  sum_score-last_sum_score); 
#endif
			  break;
	        	default: print_error("PrtnModelScan( ) mode input error");
		    }
                    last_sum_score = sum_score;
                    if(prob <= rptsEval) repeats = r;
                    else break;
                 }
	   }
	   totN = repeats*F->N;
	   UnXSeq(E); // remove X'ed out regions 
// THIS IS STILL NOT WORKING CORRECTLY; I'D LIKE TO COMPUTE Spouge BOTH WITH
// AND WITHOUT IGNORED COLUMNS BUT DON"T WANT TO USE MORE COMPUTE TIME.  ASK 
// SPOUGE IF CAN BE DONE AS A SIDE JOB (I WOULD THINK IT CERTAINLY COULD).
	   switch (mode) {
		case 'D': LocalAlnSeqSMatrix(len,SeqPtr(E),totN,F->fsM,p); break;
		case 'd': LocalGapFuncAlnSMatrix(len,SeqPtr(E),totN,F->fsM,p,F->scoreGap);
			    break;
	        case 'G': AlnSeqSMatrix(len,SeqPtr(E),totN,F->fsM,p); break;
	        case 'a': 
	        case 'g': GapFuncAlnSeqSMatrix(len,SeqPtr(E),totN,F->fsM,p,F->scoreGap);
			    break;
	        case 'O': AlnSeqOverlapSMatrix(len,SeqPtr(E),totN,
				F->overlap,F->fsM,p); break;
	        case 'o': GapFuncAlnSeqOverlapSMatrix(len,SeqPtr(E),
				totN,F->fsM,p,F->scoreGap,F->overlap); break;
		case 'e': GapFuncAlnSeqOverlapSMatrixCterm(len, SeqPtr(E), totN,
       				F->fsM, p, F->scoreGap,F->overlap); break;
	        default: print_error("PrtnModelScan( ) mode input error");
	   }
	   p[0]=1; p[totN+1]=len+1; /** set ends **/
           for(m = 0; m < totN; ){
		   m0 = (m % F->N) + 1; m++;
                   field = p[m+1] - (p[m-1] + F->w[m0-1]) - F->w[m0] + 1;
#if 1  // temperary fix for overlap method...need to add overlap value!!
		if(field < 1) field =1;
#endif 
                   prob = PvalWModel(SeqPtr(E),p[m],F->M[m0]);
                   pv[m] = (float)-log10(prob*(double)field);
           }
        }
	if(repeats == 0) return repeats;
	else if(mode != F->mode) return -repeats;
	else return repeats;
}

Int4	ComparePrtnModelSW(e_type E,UInt8 total,double maxEval, double *pvalue, ptm_typ F)
{ return ComparePrtnModelSW(E,total,maxEval, pvalue, 10,5,F); }

Int4	ComparePrtnModelSW(e_type E,UInt8 total,double maxEval, double *pvalue, 
	Int4 open, Int4 extend,ptm_typ F)
{ return ComparePrtnModelSW(E,total,maxEval,pvalue,open,extend,F,' '); }

Int4	ComparePrtnModelSW(e_type E,UInt8 total,double maxEval, double *pvalue, 
	Int4 open, Int4 extend,ptm_typ F,char mode)
// return 1 if a significant hit is found, otherwise zero
{
	Int4	m,m0,score,n,len = LenSeq(E),cutoff=0;
	a_type	A=F->A;
	double	factor,pval,pval2,mask_drop=10.0,prob,target;
	BooLean	masked;

#if 0
	if(F->Spg == NULL) {
	   print_error("Spouge gap functions deactivated");
	   F->Spg = MkSpouge(nAlpha(A)+1,F->freq,F->N,F->sM,
			F->maxLength+F->tot_len+1, F->observedGap,F->mode);
	}
#endif
	if(F->nonglobular) {
           masked= ProcessSeqPSeg(45,3.4,3.75,100,E,A);
	   if(F->segmask) { /* mask compositionally biased regions. */
		masked= (masked || ProcessSeqPSeg(14,2.2,2.5,100,E,A));
	   }
	} else if(F->segmask) { /* mask compositionally biased regions. */
	        masked=ProcessSeqPSeg(14,2.2,2.5,100,E,A); 
	} else masked=FALSE;
	n = len - F->tot_lenOL; m = F->N;
	if(n < 0) n=0;
	factor = bico(n+m,m)*(double)total/(double) len; /* Altschul's way */

	target = maxEval/((double)total/(double)len); 
	if(mode != 'C'){
	  fprintf(stdout,"\n===================================================\n");
          PutSeqInfo(stdout,E);
	  fprintf(stdout,"\n---------------------------------------------------\n");
	}
	score = PutSeqAlnSMatrixSW(stdout,open,extend,LenSeq(E),SeqPtr(E),
				m,F->sM,F->scoreGap,OffSetSeq(E),mode);
	if(mode == 'C'){ 
	  // fprintf(stdout,">"); PutSeqInfo(stdout,E);
	} else { 
	   fprintf(stdout,"\n"); PutSeqInfo(stdout,E); 
	   fprintf(stdout,"\n===================================================\n");
	}
#if 0
	// cutoff=SpougeThreshold(F->Spg,target,MAXIMUM(Int4,len,F->tot_len));
        // score=AlnSeqSMatrixSW(10,5,LenSeq(E),XSeqPtr(E),m,F->sM,NULL,cutoff);
	// score=GappedAlnSeqSMatrixSW(10,5,LenSeq(E),SeqPtr(E),m,F->sM,F->scoreGap);
        // pval = factor*SMatrixProb(score, F->smx);
	// *pvalue = pval;
        if(score >= cutoff){
             PutSeqInfo(stdout,E);
             printf("\n !gapped score = %d (p = %g)\n",score,pval);
             printf("********************************\n\n");
	     return 1;
        } 
#endif

	return score;
}

// IMPORTANT: check to make sure that this is consistent with 
// default_overlap in sequence.h for static unsigned char *alloc_seq(Int4 length)!!
const char default_core=4;
const char default_overlap=9;

void	AllocatePrtnModel(Int4 *t_used,sma_typ MA, Int4 maxlen, ptm_typ F)
{
	Int4	i,j,c,r,m,totN,score;
	Int4	end,start;
	a_type	A=F->A;

	totN = NmaxPrtnModel(F);
	NEW(F->sM,totN+3,smx_typ); NEW(F->fsM,totN+3,smx_typ); 
	NEW(F->w,F->N+2,Int4); 
	for(F->w[0]=0, F->tot_len=0,m=1; m<=F->N; m++) { 
		F->w[m] = LenWModel(F->M[m]); F->tot_len += F->w[m];
	}
	if(F->mode == 'o' || F->mode == 'O' || F->mode == 'e'){
	   char	OL;
	   Int4	core;
	   NEW(F->overlap,totN+3 ,char);
	   F->tot_lenOL=0;
	   F->overlap[0] = default_overlap; 
	   for(m=1; m<=F->N; m++) { 
                core = F->w[m] - F->overlap[m-1] - default_overlap;
                core = MAXIMUM(char,core,default_core);
		OL = (F->w[m] - core)/2;
		if(OL < 0) { OL = 0; core = F->w[m]; }
	   	F->overlap[m-1] = MINIMUM(char, OL,F->overlap[m-1]);
		OL = MINIMUM(char,default_overlap,OL);
	   	F->overlap[m] = OL;
	   }
#if 1
	 // fprintf(stderr,"overlap[0] = %d\n",F->overlap[0]);
	 F->totalOL=F->overlap[0];
	 for(m=1; m<=F->N; m++) { 
	   core = F->w[m] - F->overlap[m];
	   if(m==1) core -= F->overlap[0];
	   F->tot_lenOL += core;
	   F->totalOL+=F->overlap[m];
	   // fprintf(stderr,"core[%d] = %d;len = %d\n",m,core,F->w[m]);
	   // fprintf(stderr,"overlap[%d] = %d\n",m,F->overlap[m]);
	   if(F->maxrpts > 1){
	     for(r=1; r < F->maxrpts; r++){ 
		j=m+r*F->N;
		if(m==1 || m==F->N)
		  F->overlap[j]=MINIMUM(char,F->overlap[1],F->overlap[F->N]);
		else F->overlap[j]=F->overlap[m];
	   	// fprintf(stderr,"overlap[%d] = %d\n",j,F->overlap[j]);
	     }
	   }
	 }
	 // fprintf(stderr,"total overlap = %d\n",F->totalOL);
	 F->maxLength=maxlen + F->totalOL;
#endif
	} else {
	   F->overlap = NULL;
	   F->tot_lenOL=F->tot_len; F->totalOL=0;
	   F->maxLength=maxlen;
	}
	/**** CONSTRUCT FULL SCORING MATRIX. ****/
	F->smx = MkSMatrix(2.0,F->tot_len,F->freq,A);
	for(i=1,m=1; m<=F->N; m++) { 
	    for(j=1; j<=F->w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
			score = CellScoreWModel(c, j, F->M[m]); 
			SetSMatrix(c,i,score,F->smx);
		}
	    }
	    F->sM[m] = GetSMatrixWModel(F->M[m]);
	    if(F->M2[m] == NULL) F->fsM[m] = F->sM[m];
	    else F->fsM[m] = GetSMatrixWModel(F->M2[m]);
	    for(r=1; r < F->maxrpts; r++){ 
		j=m+r*F->N; F->sM[j]=F->sM[m];
		F->fsM[j]=F->fsM[m]; 
	    }
	}
	/********************** Compute Gaps ****************************/
	if(F->gapfunct){
          F->totGaps = nseqSMA(MA);
	  NEWP(F->observedGap, F->N+2, Int4);
	  NEW(F->observedGap[0], F->maxLength+3, Int4);;
	  for(m=1; m< F->N; m++) {
	     F->observedGap[m]=GetGapsSMA(F->maxLength,1.3,t_used[m],
				t_used[m+1],MA);
	  }
	  NEW(F->observedGap[F->N], F->maxLength+3, Int4);;
#if 0	// make Cterm statistics
	  F->observedGap[F->N][0] = F->totGaps;
#endif
	}
	/********************* End Compute Gaps *************************/
	if(F->mode != 'G' && F->mode != 'O'){
	   print_error("Spouge gap functions deactivated");
#if 0
	   F->Spg = MkSpouge(nAlpha(A)+1,F->freq,F->N,F->sM,F->maxLength+1,
			F->observedGap,F->mode);
	   if(islower(F->mode)){
	   	NEWP(F->scoreGap,totN+3,Int4); 
	   	F->scoreGap[0] = GapScoreSpouge1(0,F->Spg);
	   	for(m=1; m<=F->N; m++) { 
	      	   F->scoreGap[m] = GapScoreSpouge1(m,F->Spg);
	      	   for(r=1; r < F->maxrpts; r++){ 
			j = m + r*F->N; 
			F->scoreGap[j] = GapScoreSpouge1(m,F->Spg); 
		   }
		}
	   } else F->scoreGap=NULL;
	   F->max_gap_score = MaxGapScoreSpouge(F->Spg);
	   fprintf(stderr,"max_gap_score = %d\n",F->max_gap_score);
#endif
	}
}

double	PutFullSWAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	s,m=NumModelsPrtnModel(PM);

	s=PutFullSeqAlnSMatrixSW(fp,a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
	fprintf(fp,"\n\tscore = %d\n\n",s);
	return 0;
}

Int4	LocalScoreSWAlnPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m=NumModelsPrtnModel(PM);
 return ScoreGappedAlnSeqSmatrixSW(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap,TRUE);
}

Int4	ScoreSWAlnPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m=NumModelsPrtnModel(PM);
 return ScoreGappedAlnSeqSmatrixSW(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap,FALSE);
}

double	PutSWAlnPrtnModelRpts(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM)
// same as PutSWAlnPrtnModel( ) but with N repeats 
{
	Int4	s,m=NmaxPrtnModel(PM);

	s=PutSeqAlnSMatrixSW(fp,a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
	fprintf(fp,"\n\tscore = %d\n\n",s);
	return 0;
}

BooLean	PutSeqAlnPrtnModel(FILE *fp, char *operation, Int4 start, e_type E, ptm_typ PM)
{
	// Int4	s,m=NumModelsPrtnModel(PM);  make the same as with repeats.
	Int4	s,m=NmaxPrtnModel(PM);

	fprintf(fp,"\n");
	PutSeqInfo2(fp,E);
	PutGappedSeqAlnSMatrix(fp, operation, OffSetSeq(E)+start-1, LenSeq(E)-start+1,
        	SeqPtr(E)+start-1, m, PM->sM);
	fprintf(fp,"\n\n\n");
	return TRUE;
}

Int4	PutSWAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	// Int4	s,m=NumModelsPrtnModel(PM);  make the same as with repeats.
	Int4	s,m=NmaxPrtnModel(PM);

	s=PutSeqAlnSMatrixSW(fp,a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
	// fprintf(fp,"\n\tscore = %d\n\n",s);
	return s;
}

double	PutSampledSWAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	s,m=NmaxPrtnModel(PM);

	s=PutSampledSeqAlnSMatrixSW(fp,a,b,LenSeq(E),SeqPtr(E),
		OffSetSeq(E),m,PM->sM,PM->scoreGap);
	fprintf(fp,"\n\tscore = %d\n\n",s);
	return s;
}

double	ShortGapAlnPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m,nmod=NumModelsPrtnModel(PM),score,score2,*posG,*posL,*posSW;
	Int4	end,lenM,s,len,max_s,ngaps=3,*pos,p;
	double	ave=0.;
	unsigned char	*seq;
	smx_typ	M[3];

	if(fp != NULL){ PutSeqID(fp,E); fprintf(fp,"\n"); }
	len = LenSeq(E); seq = SeqPtr(E);
	NEW(pos,nmod+3,Int4); NEW(posG,nmod+3,Int4); 
	NEW(posL,nmod+3,Int4); NEW(posSW,nmod+3,Int4);
	score2 = BestAlnSeqSMatrix(len,seq,nmod,PM->sM,posG);
	s = BestLocalAlnSeqSMatrix(len, seq, nmod, PM->sM, posL);
	score=AlnSeqOverlapSMatrix(LenSeq(E),SeqPtr(E),nmod,PM->overlap,PM->sM,pos);
	if(fp != NULL) fprintf(fp,"total ungapped score = %d\n",score);
        for(score = 0, m=1; m <= nmod; m++){
	     M[1] = PM->sM[m]; lenM = LenSMatrix(PM->sM[m]);
	     len = lenM + ngaps; 
	     seq = SeqPtr(E) + pos[m] - 1; 
	     // seq = SeqPtr(E) + pos[m] - 4;  // TEST  
	     // seq = SeqPtr(E) + pos[m] + 1;  // TEST  
             // s=AlnSeqSMatrixSW(a,b,len,seq,1,M,NULL,-999);
             s=GapXDropSMatrix(a,b,len,seq,PM->sM[m],&p,ngaps);
	     if(fp != NULL) fprintf(fp,"gapped score: %d\n",s); 
	     score += s;
	     s = ScoreSMatrix(SeqPtr(E),pos[m],M[1]);
	     if(fp != NULL) fprintf(fp,"ungapped score: %d\n",s); 
	     if(fp != NULL) 
		fprintf(fp,"pos: %d; offset = %d; posG = %d; posL = %d; posSW = %d\n",
					pos[m],p,posG[m],posL[m],posSW[m]); 
	     if(fp != NULL) fprintf(fp,"............................................\n");
	}
	if(fp != NULL) fprintf(fp,"total gapped score = %d\n",score);
	if(fp != NULL) fprintf(fp,"\n");
	free(pos); free(posG); free(posL); free(posSW);
	return score;
}

Int4    GapScorePrtnModel(Int4 block, Int4 gaplength, ptm_typ PM)
{
        if(PM->gapfunct == FALSE) print_error("gap scores undefined");
        if(block < 1 || block > NumModelsPrtnModel(PM)) 
                print_error("block undefined");
        if(gaplength > PM->maxLength) print_error("gap length too long");
        return PM->scoreGap[block][gaplength];
}

void	PutGapScoresPrtnModel(FILE *fp, ptm_typ PM)
{
	Int4	m,j,**gapscore=PM->scoreGap;

	if(gapscore != NULL){
	  for(m=1; m <= NumModelsPrtnModel(PM); m++){
            for(j=1; j<= PM->maxLength; j++) {
                if(gapscore[m][j] == SHRT_MIN) break; // in gap limit?
                fprintf(stdout,"gapscore[%d][%d] = %d\n", m,j,gapscore[m][j]);
            }
	  }
	}
}

char    *RptsSampleGapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m=NmaxPrtnModel(PM);
	return SampleGapOperationsSMatrix(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
}

char    *RptsGapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m=NmaxPrtnModel(PM);
	return GapOperationsSMatrix(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
}

char    *SampleGapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m=NumModelsPrtnModel(PM);   // one repeat...
	return SampleGapOperationsSMatrix(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
}

char    *GapOperationsPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	m=NumModelsPrtnModel(PM);   // one repeat...
	return GapOperationsSMatrix(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
}

void	CompareScoresPrtnModel(Int4 a, Int4 b, e_type E, ptm_typ PM)
// compare the p-values (or log-odds scores) for sequence hits w/ and w/o gaps
{
        Int4    m;
	char	*operation;

#if 0	// work in progress...
        m = NumModelsPrtnModel(PM);
	operation=GapOperationsSMatrix(a,b,LenSeq(E),SeqPtr(E),m,PM->sM,PM->scoreGap);
	free(operation);
#endif
}

double	SWHistPrtnModel(FILE *fp,Int4 a, Int4 b, e_type E, ptm_typ PM)
{
	Int4	i,m,n,nmod=NumModelsPrtnModel(PM),score;
	Int4	end,lenM,s,len,max_s,ngaps=0;
	double	ave=0.;
	unsigned char	*seq,r;
	smx_typ	M[3];
	h_type  H;
	BooLean	found;
	unsigned char seg[200];

	fprintf(fp,">simulated seq\n");
	for(i=1,n=1; n<=10; n++){
		r=GetBackGroundWModel(PM->M[1]);
		fprintf(fp,"%c",AlphaChar(r,PM->A));
		if(i%50==0) fprintf(fp,"\n");
	}
        for(m=1; m <= nmod; m++,i++){
		// n=MaxSegWModel(seg, PM->M[m]);
		n=GetSegWModel(seg, PM->M[m]);
		for(s=1; s<=n; s++){
			fprintf(fp,"%c",AlphaChar(seg[s],PM->A));
			if(i%50==0) fprintf(fp,"\n");
		}
		for(n=1; n<=10; n++){
			r=GetBackGroundWModel(PM->M[m]);
			fprintf(fp,"%c",AlphaChar(r,PM->A));
			if(i%50==0) fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n");
#if 0
	if(fp != NULL){ PutSeqID(fp,E); fprintf(fp,"\n"); }
        for(m=1; m <= nmod; m++){
	      M[1] = PM->sM[m]; lenM = LenSMatrix(PM->sM[m]);
	      H = Histogram("smatrix local scores",-100,200,2.0);
	      len = lenM + ngaps - 1; end = LenSeq(E) - lenM; 
	      seq = SeqPtr(E) - ngaps; 
	      found = FALSE;
	      max_s=AlnSeqSMatrixSW(a,b,(LenSeq(E)+2*ngaps),seq,1,M,NULL,999999);
	      for(Int4 i=-(ngaps-1); i < end; i++){
		seq = SeqPtr(E); seq += i-1; 
        	s=AlnSeqSMatrixSW(a,b,len,seq,1,M,NULL,max_s);
		IncdHist(s, H);
		if(s == max_s) { 
	      		found = TRUE;
			if(fp != NULL) fprintf(fp,"pos: %d\n",i); 
			s = ScoreSMatrix(SeqPtr(E), i, M[1]);
			if(fp != NULL) fprintf(fp,"ungapped score: %d\n",s); 
		}
	      }
	      if(fp != NULL && found) PutHist(fp,60,H); NilHist(H);
	      if(fp != NULL && found) fprintf(fp,"============================================\n");
	}
	if(fp != NULL) fprintf(fp,"\n");
	return ave/(double) nmod;
#endif
	return 0.0;
}

double	LocalHistPrtnModel(FILE *fp,e_type E, ptm_typ PM)
{
	Int4	m,nmod;
	double	ave=0.;

        nmod = NumModelsPrtnModel(PM);
	if(fp != NULL){ PutSeqID(fp,E); fprintf(fp,"\n"); }
        for(m=1; m <= nmod; m++) ave+=LocalHistSeqSMatrix(fp,LenSeq(E),SeqPtr(E),PM->sM[m]);
	if(fp != NULL) fprintf(fp,"\n");
	return ave/(double) nmod;
}

double	HistPrtnModel(FILE *fp,e_type E, ptm_typ PM)
{
	Int4	m,nmod;
	double	ave=0.;

        nmod = NumModelsPrtnModel(PM);
	if(fp != NULL){ PutSeqID(fp,E); fprintf(fp,"\n"); }
        for(m=1; m <= nmod; m++) {
		ave+=HistSeqSMatrix(fp,LenSeq(E),SeqPtr(E),PM->sM[m]);
	}
	if(fp != NULL) fprintf(fp,"\n");
	return ave/(double) nmod;
}

Int4	PutCRSSeqAlnSMatrixSW(FILE *fp, char *operation, Int4 n2, unsigned char *seq2, 
		UInt4 offset, Int4 J, Int4 nmod, smx_typ *M,char color, 
		Int4 m_use, Int4 rpt_use, fm_type *FM)
// Output a *crs file 
{
        Int4    o,m,i,j,k,r1,s,t,v,n1,score,*pos,p,p0;
	Int4	j2,s1,s0,g,gs,g_opt,rpt;
	unsigned char	*seq1,*seq;
	char	*out[3];
	short	*mtf;
	a_type	A = SMatrixA(M[1]);

	/** get total length of profile **/
	assert(m_use > 0 && m_use <= nmod);
	assert(nmod <= 26);
	for(n1=0, m = 1; m <= nmod; m++){ n1 += LenSMatrix(M[m]); }
	/** get concensus sequence for smatrix **/
	NEW(seq1, n1+3, unsigned char); 
	NEW(mtf, n1+3, short); 
	for(s=0, m=1; m <= nmod; m++){
	    MaxSegSMatrix(seq1+s, M[m]);
	    for(i=1; i<= LenSMatrix(M[m]); i++){ s++; mtf[s] = m; }
	}
	/*** 3. Trace back step. ***/
	NEW(seq, J+3, unsigned char); 
	MEW(out[0],J+3,char); MEW(out[1],J+3,char);
	MEW(out[2],J+3,char); NEW(pos,J+3,Int4); 

	// operations = i,I,d,D,m,M,E;
	fprintf(fp,"\n");
fprintf(fp,"#off //************** Motif %c(%d) *****************************************\n",
				   (char)(m_use + 'A' - 1),rpt_use);
	m=1;  k=1; score=0;
	for(p=p0=1,rpt=o=j=i=1; operation[o] != 'E'; o++){ 
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
			score+=ValSMatrix(k,UndefAlpha(A),M[m]);
			out[1][p] = AlphaChar(seq1[i],A);
			out[0][p] = ' '; out[2][p] = '-';
			seq[p0]=0; p0++; pos[p]=j;
			p++; i++;  k++;
			break;
		    case 'I':   // Gap ('-') is in profile
			out[1][p] = '-'; out[0][p] = ' ';
			out[2][p] = AlphaChar(seq2[j],A);
			pos[p] = j; p++; j++; break;
		    case 'i': j++; break;
		    default: 
			fprintf(stderr,"operation[%d] = %c\n",o,operation[o]);
			print_error("this should not happen"); 
			break;
		}
		if(k > LenSMatrix(M[m])){ // 4. Print out alignment 
		  if(m == m_use && rpt == rpt_use){
		   fprintf(fp,"%5.1f: ", -log10(SMatrixProb(score,M[m])));
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[1][s]); 

	   	   fprintf(fp,"\n       ");
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[0][s]); 
	   	   fprintf(fp,"\n");

		   if(m > 1){ 
	           	   char tmpstr[100];
			   sprintf(tmpstr,"(%d)",pos[1]-pos[0]-1);
			   fprintf(fp,"%6s ",tmpstr);
		   } else fprintf(fp,"       ");
	   	   for(s=1; s < p ; s++) fprintf(fp,"%c",out[2][s]);
	   	   fprintf(fp," %d-%d\n",pos[1]+offset,j-1+offset); 
		   if(FM && FM[m_use]) { fprintf(fp,"\n"); PutFModelShort(fp, FM[m_use]); }
fprintf(fp,"#on //*****************************************************************\n");
	   	   fprintf(fp,"%d-%d.%c\n",pos[1]+offset,j-1+offset,color); 
		  }
		  // fprintf(stderr,"m=%d;rpt=%d;k=%d;p=%d\n",m,rpt,k,p);
		  pos[0]=j-1;
		  if(m < nmod) { m++; k=1; p=p0=1; score=0; } 
		  else { m=k=1; rpt++; p=p0=1; score=0; i=1; }
#if 0
		  else break;
#endif
		}
	}
	// fprintf(fp,"(%d)\n",n2-pos[0]);
	fprintf(fp,"\n");
	/*** 5. free allocated memory ***/
	free(seq1); free(mtf); free(pos);  free(seq);
	free(out[0]); free(out[1]); free(out[2]);
	return o;
}


Int4	PlotCRSPrtnModelSW(FILE *fptr, Int4 open,Int4 extend, double info_cutoff, 
			e_type E, fm_type *FM, ptm_typ PM)
{ return PlotCRSPrtnModelSW(fptr,open,extend,info_cutoff,E,FM,PM,0); }

Int4	PlotCRSPrtnModelSW(FILE *fptr, Int4 open,Int4 extend, double info_cutoff, 
			e_type E, fm_type *FM, ptm_typ PM,char *InputColors)
// 'M' or 'D' = Start of Motif model 
{
        Int4    N,i,j,k,s,nr,n,exp,nc,m,o;
        FILE    *sfptr;
	float	**info;
	float	*zscore;  	// score at model position s
        unsigned char *seq;	// amino acid sequence sequence code at s
        char    r,c;
	char	*operation;
	Int4	color;
	char    *resname,res_char;
        char    colors[] = "BMROYGCBMROYGCB",*Colors;
	Int4	numBlks,rpt;
	a_type	A = PrtnModelA(PM);
	Int4 os = OffSetSeq(E);
	Int4 Score,Start,number_colors=7;

	if(InputColors){
	   Colors=InputColors;
	   number_colors=strlen(InputColors)-1;
	} else Colors=colors;
	// operation=RptsSampleGapOperationsPrtnModel(open, extend, E, PM);
	operation=RptsGapOperationsPrtnModel(open, extend, E, PM);
	// std::cerr << operation;
	Start = 1;
	numBlks = NumModelsPrtnModel(PM);
	if(FM){
		NEWP(info,numBlks+3,float);
		for(i=1; i <= numBlks; i++){
			assert(FM[i]); info[i]=InfoFModel(FM[i]);
		}
	} else info = InfoPrtnModel(250, PM);
#if 0
	fprintf(fptr,"File1=junk.pdb:A    // ( Angstroms)\n");
	fprintf(fptr,"#file1:R=x;T=x.\n");
	fprintf(fptr,"#file1:R=x;T=x.\n\n");
	fprintf(fptr,"1-1000.W\necho=trace 50\n\n");
#else
	fprintf(fptr,"1-1000.W\n\n");
#endif
#if 0
	Int4	score=PutSeqAlnSMatrixSW(fptr,open,extend,LenSeq(E),SeqPtr(E),
			NmaxPrtnModel(PM),PM->sM,PM->scoreGap);
	// fprintf(fptr,"\n\tscore = %d\n\n",s);
	// PutSWAlnPrtnModel(fptr,open, extend,E,PM);
#else
	// NumModelsPrtnModel(PM) == number blocks in one repeat.
#endif
	seq = SeqPtr(E)-os;
	Int4 score;
        for(m=0,rpt=1,o=k=1,s=os+Start; (c=operation[o]) != 'E'; o++){
	   switch(c){
	      case 'M': 	// Matching start position.
		k=1; m++; 
		if(m > numBlks) { m = 1; rpt++; }
		if(numBlks > 1) { color = m%number_colors; }
		else { color = rpt%number_colors; }
#if 0
		score=PutCRSSeqAlnSMatrixSW(fptr, operation, LenSeq(E),SeqPtr(E),
			OffSetSeq(E),strlen(operation),NmaxPrtnModel(PM),PM->sM,
			Colors[color],m,rpt,FM);
#else
		{
		 fm_type *dummy_FM=0;
		 score=PutCRSSeqAlnSMatrixSW(fptr, operation, LenSeq(E),SeqPtr(E),
			OffSetSeq(E),strlen(operation),numBlks,PM->sM,
			Colors[color],m,rpt,dummy_FM);
		}
#endif
		// fprintf(fptr,"// motif %d\n%d-%d.%c\n\n",m,s,s,Colors[color]);
	      case 'm':		// match within a block.
	        {
		 res_char = AlphaChar(seq[s],A);
#if 1
		 if(info[m][k] >= info_cutoff){ fprintf(fptr,"%c%d.M,",res_char,s); }
#elif 0
		 if(info[m][k] >= info_cutoff){
		   fprintf(fptr,"%c%d.M  // %.2f bits\n",res_char,s,info[m][k]);
		 }
#else
		 if(info[m][k] < info_cutoff) fprintf(fptr,"#");
		 fprintf(fptr,"%c%d.M  // %.2f bits\n",res_char,s,info[m][k]);
#endif
		} s++; k++; 
	        break;
	      case 'i': case 'I': s++; 
		break;
	      case 'D': k=1; m++; 
		if(m > numBlks) { m = 1; rpt++; }
		if(numBlks > 1) { color = m%number_colors; }
		else { color = rpt%number_colors; } 
	      case 'd': k++; 
		break;
	      default: print_error("input error");
	   }
	   // fprintf(fptr,"%c %d(%d): s=%d %c (%.2f)\n",c,m,k,s,AlphaChar(seq[s],A),info[m][k]);
        } 
	// fprintf(fptr,"\nclose:1.\n\n");
	fprintf(fptr,"\n\n");
	for(m=1; m <= NumModelsPrtnModel(PM); m++) free(info[m]); free(info);
	free(operation);
	return 0;
}

