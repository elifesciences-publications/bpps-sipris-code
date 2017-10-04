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

/* oscan.c - ordered scan program */
#include "oscan.h"

osn_typ MakeOScan(double pseudo, char *snfile, Int4 *counts, a_type A, 
	Int4 maxrpts, char method, char mode, double ecut,double Ecut,Int4 maxLength,
	float minmap, BooLean weight)
{ return make_oscan(pseudo,snfile,counts,A,maxrpts,method,ecut,
			Ecut,mode,maxLength,minmap,weight); }

osn_typ MkOScan(double pseudo, char *snfile, Int4 *counts, a_type A, 
	Int4 maxrpts, char method,double ecut,double Ecut)
{ return make_oscan(pseudo,snfile,counts,A,maxrpts,method,ecut,Ecut,'G',1000,0.0,FALSE); }

osn_typ make_oscan(double pseudo, char *snfile, Int4 *counts, a_type A, 
	Int4 maxrpts, char method,double ecut,double Ecut, char mode,
	Int4 maxLength, float minmap, BooLean weight)
{
	Int4	i;
	osn_typ F;

	NEW(F,1,oscan_type);
	NEW(F->M,MAX_NUM_MODELS+1,wm_type);
	NEW(F->M2,MAX_NUM_MODELS+1,wm_type);
	F->snfile = AllocString(snfile);
	if(islower(mode)) F->gapfunct= TRUE;
	else F->gapfunct= FALSE;
	F->shuffle=FALSE;
	F->segmask=TRUE;
	F->mode=mode;
	F->weights=weight;
	F->pdbfa=NULL;
	F->maxLength=maxLength;
	F->scoreGap=NULL;
	F->totGaps=0;
	F->dfp = NULL;
	F->gaps = FALSE;
	F->fp_gaps=NULL;
	F->recomb = NULL;
	F->permute=0;
	if(ecut <= 0.0 || Ecut <= 0.0)
		print_error("maximum E-value must be >= 0.0");
	F->maxEval = ecut; F->log10ecut = log10(ecut);
	F->singleEval = Ecut;
	F->A = A; 
	F->N = 0; 
	F->method = method; 
	if(!ReadOScanCMSA(F->snfile,minmap,pseudo,counts,F))
		ReadOScan(F,pseudo,F->snfile,counts);
	F->maxrpts= MAXIMUM(Int4,1,maxrpts);
	return F;
}

void	SetEvalOScan(double ecut, double Ecut, osn_typ F)
{
	if(ecut <= 0.0 || Ecut <= 0.0)
        	print_error("maximum E-value must be >= 0.0");
        F->maxEval = ecut; F->log10ecut = log10(ecut);
	F->singleEval = Ecut;
}

void	PDBhomOScan(char *infile, osn_typ F)
/*** Get a pdb structural file to create threaded sequences ***/
{ F->pdbfa=AllocString(infile); }

void    GapsOScan(osn_typ F)
{
	F->gaps=TRUE;
	F->fp_gaps=open_file(F->snfile,".gsq","w"); 
}

FILE	*OpenPermuteOScan(osn_typ F)
{ F->recomb=open_file(F->snfile,".pmt","w"); return F->recomb; }

FILE	*OpenDatabaseOScan(osn_typ F)
{ F->dfp=open_file(F->snfile,".dbs","w"); return F->dfp; }

void	NilOScan(osn_typ F)
{
	Int4	i;

	for(i=1;i<=F->N; i++) { 
		NilWModel(F->M[i]); 
		if(F->M2[i] != NULL) NilWModel(F->M2[i]); 
	}
	free(F->M2);
	free(F->M); free(F->freq); 
	free(F->snfile);
	if(F->dfp != NULL) fclose(F->dfp);
	if(F->recomb != NULL) fclose(F->recomb);
	if(F->scoreGap != NULL){
		for(i=0;i<=F->N; i++) free(F->scoreGap[i]);
		free(F->scoreGap);
	}
	if(F->gaps) fclose(F->fp_gaps);
	if(F->pdbfa != NULL) free(F->pdbfa);
	free(F);
}

BooLean	ReadOScanCMSA(char *msafile, float minmap, double pseudo,Int4 *counts, 
	osn_typ F)
{ 
    Int4		m,n,s,s0,t,N,*t_used,i,length;
    UInt8		total;
    char		*null0=NULL,*null=NULL;
    unsigned char	*seq=NULL;
    BooLean		ignore,use_null=TRUE;
    double		*wt;
    a_type		A = F->A;
    sma_typ		MA;
    
    if((MA = ReadSMA(msafile)) == NULL) return FALSE;   

    if(isupper(F->method)){ use_null=FALSE; F->method=tolower(F->method); }
    NEW(F->freq,nAlpha(A)+2,double);
    for(total=i=0; i<=nAlpha(A); i++) total += counts[i];
    for(i=0; i<=nAlpha(A); i++) {
        F->freq[i]= (double)counts[i]/(double)total;
    }
    F->total = (double) total;
    
    if(F->weights) wt = WeightsSMA(MA); 
    else {
	NEW(wt,nseqSMA(MA)+2,double);
	for(n=1; n<= nseqSMA(MA); n++) wt[n]=1.0;
    }
    if(ntypSMA(MA) >= MAX_NUM_MODELS)
	print_error("ReadOScanCMSA( ) input > MAX_NUM_MODELS!");
    NEW(t_used,ntypSMA(MA)+3,Int4);
    for(F->N=m=0,t=1; t<=ntypSMA(MA); t++){
        length = lengthSMA(t,MA);
	if(null != NULL) free(null); NEW(null, length+2, char);
	null0 = nullSMA(t,MA);
	ignore = FALSE;
	for(s0=0,s=1; null0[s] != 0; s0++,s++){
	    switch(null0[s0]){
	      case '.': if(use_null) null[s]='.'; else null[s]='*'; break;
	      case '^': null[s]='^'; ignore=TRUE; break;
	      case '*': null[s]='*'; break;
	      case '!': null[s]='!'; break;
	      default: print_error("ReadOScanCMSA( ) error 1"); break;
	    }
        }
	if(fieldmapSMA(t,MA) >= minmap){
	  F->N++; m++; t_used[m]=t; 
	  F->M[m]=MkWModel(null,length,pseudo,F->freq,A);
	  SetMethodWModel(F->method,F->M[m]);
	  if(ignore) {
	    F->M2[m]=MkWModel(NULL,length,pseudo,F->freq,A);
	    SetMethodWModel(F->method,F->M2[m]);
	  } else F->M2[m] = NULL;
          for(n=1; n<= nseqSMA(MA); n++){
	    seq=seqSMA(t,n,MA);
	    if(probSMA(n,t,MA) >= 1.0){
		Add2WModel(seq,1,wt[n],F->M[m]);
	    	if(ignore) Add2WModel(seq,1,wt[n],F->M2[m]);
	    }
	  }
        } 
    }
    free(wt);
    if(null != NULL) free(null);
    /*********************** Compute Gaps *******************************/
    if(F->gapfunct){
        F->totGaps = nseqSMA(MA);
	NEWP(F->scoreGap, F->N+2, Int4);
	NEW(F->scoreGap[0], F->maxLength+3, Int4);;
	for(m=1; m< F->N; m++) {
	   F->scoreGap[m]=GetGapsSMA(F->maxLength,1.3,t_used[m],t_used[m+1],MA);
	}
	NEW(F->scoreGap[F->N], F->maxLength+3, Int4);;
    }
    /*********************** End Compute Gaps ***************************/
    free(t_used);
    NilSMA(MA);
    return TRUE;
}

void	ReadOScan(osn_typ F, double pseudo,char *snfile, Int4 *counts)
{
	FILE	*fptr;
	Int4	m,m2,n,i,length,number,pos;
	UInt8	total;
	unsigned char	**seq;
	char	c = '\n',*null;
	BooLean	use_null=TRUE, verbose=FALSE,ignore;
	a_type	A = F->A;

   if(counts == NULL) print_error("ReadOScan( ) input error 1");
   if(isupper(F->method)){ use_null=FALSE; F->method=tolower(F->method); }
   NEW(null,MAX_BLOCK_LENGTH+2,char);
   NEWP(seq,MAXSCN_BLOCK_SIZE+2,unsigned char);
   for(i=0; i<=MAXSCN_BLOCK_SIZE; i++) NEW(seq[i],MAX_BLOCK_LENGTH+2,unsigned char);
   NEW(F->freq,nAlpha(A)+2,double);
   for(total=i=0; i<=nAlpha(A); i++) total += counts[i];
   for(i=0; i<=nAlpha(A); i++) {
	F->freq[i]= (double)counts[i]/(double)total;
	if(verbose) fprintf(stderr,"freq[%d] = %d/%d = %g\n",
		i,counts[i],total,F->freq[i]);
   }
   F->total = (double) total;
   if((fptr = fopen(snfile,"r")) == NULL) {
        fprintf(stderr,"Could not open file \"%s\"\n",snfile);
        print_error("File does not exist!");
   }
   for(c=' ',length = 0,F->N=0,m=1;c!=EOF && c!='@'; m++){
	ignore=FALSE;
        if(verbose) fprintf(stderr,"\n"); /****/
        for(; isspace(c); c=fgetc(fptr))	/** GO TO FIRST CHARACTER **/;
	if(c == '*' || c == '^'){			/** fragmentation row **/
		length = 0;
		do {
                   if(c=='^') {  /*** ignore these positions regardless ***/
			ignore = TRUE;
			length++; null[length] = c;
			if(verbose) fprintf(stderr,"^");
                   } else if(c=='!') {
			length++; null[length] = c;
			if(verbose) fprintf(stderr,"!");
                   } else if(c=='*') {
			length++; null[length] = c;
			if(verbose) fprintf(stderr,"*");
		   } else if(c=='.') {
			length++; 
			if(use_null) null[length] = '.';
			else null[length] = '*';
			if(verbose) fprintf(stderr,".");
		   } else if(!isspace(c)) {
                        fprintf(stderr,"illegal character -> %c",c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
		   if(length >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
                } while((c=fgetc(fptr))!='\n' && c!=EOF && c!='@');
	} else print_error("model file input error3");
	if(verbose) fprintf(stderr," length = %d\n",length);
	if(length == 0) print_error("input error");
        for(; isspace(c); c=fgetc(fptr))	/** GO TO FIRST CHARACTER **/;
	if(isalpha(c)){				/** = segment row **/
            for(i=1,pos=0;c!='*' && c!=EOF && c!='@' && c!='*';i++,pos=0) {
		do{
                   if(c=='*' || c == '^') { break; }
		   if(isalpha(c)) {
			if(islower(c)) c=toupper(c);
                        if(verbose) fprintf(stderr,"%c",c);/***/
		        if(++pos >= MAX_BLOCK_LENGTH) 
                   		print_error("reset MAX_BLOCK_LENGTH");
			seq[i][pos] = AlphaCode(c,A);
		   } else if(!isspace(c)) {
                        fprintf(stderr,"seq %d: illegal character -> %c",i,c);
                        fprintf(stderr,"\n");
			print_error("fatal error.");
                   }
                } while((c=fgetc(fptr))!='\n' && c!=EOF && c!='@');
		if(verbose) {
		    if(c == '\n' && pos > 0) fprintf(stderr," #%d\n",i);
		}
                if(i >= MAXSCN_BLOCK_SIZE)
                   print_error("aligned segments > MAXSCN_BLOCK_SIZE");
		else if(pos == 0)	i--;
		else if(pos != length) 
                   print_error("input segments have inconsistent lengths.");
	    }
	} else print_error("model file input error4");
        i--; number = i; 
	if(i < 1){
	    if(F->N==0) print_error("zero segments in snfile.");
	    else m--;	/** no segments in this pass **/
	} else if(i > 0){
		F->N++;
		if(F->N > MAX_NUM_MODELS) print_error("increase MAX_NUM_MODELS");
		F->M[m]=MkWModel(null,length,pseudo,F->freq,A);
		SetMethodWModel(F->method,F->M[m]);
		for(i=1; i<=number; i++){
			Add2WModel(seq[i],1,1.0,F->M[m]);
/*****
for(j=1; j<=length; j++){ printf("%c", AlphaChar(seq[i][j],A));
if(i%50 == 49) printf("\n"); } printf("\n\n");
/*****/
		}
		if(ignore){
		   F->M2[m]=MkWModel(NULL,length,pseudo,F->freq,A);
		   SetMethodWModel(F->method,F->M2[m]);
		   for(i=1; i<=number; i++) Add2WModel(seq[i],1,1.0,F->M2[m]);
		} else F->M2[m] = NULL;
		if(verbose){
		   fprintf(stderr,"\n\n");
		   PutWModel(stderr,F->M[m]); 
		   fprintf(stderr,"\n\n"); 
		}
	}
   }
   if(verbose) fprintf(stderr,"\n\t %d motif models\n",F->N);
   for(i=0; i<=MAXSCN_BLOCK_SIZE; i++) free(seq[i]);
   free(seq); free(null); 
   if(F->gapfunct){
	if(c != '@' || !ReadGapsOSCAN(fptr, F)){
		print_error("gap file read error in ReadGapsOSCAN( )");
	}
   }
   fclose(fptr);
}

BooLean	ReadGapsOSCAN(FILE *fptr, osn_typ F)
/*************************************************************************
 input a *.gap file for the gapped scan program.
 *************************************************************************/
{
        Int4    i,n,t,t2,N,len,mlen,end,begin,ntyp,g;
        unsigned char    *seq;
        float	p,prob;
	Int4	*gap;

	fscanf(fptr,"%d sequences; %d motifs; lengths:",&N,&ntyp);
	NEWP(F->scoreGap, ntyp+2, Int4);;
	NEW(F->scoreGap[0], F->maxLength+3, Int4);;
	F->totGaps = N;
	for(t=1; t<= ntyp; t++) {
	   NEW(F->scoreGap[t], F->maxLength+3, Int4);;
           fscanf(fptr," {%d}=%d",&t2,&mlen);
        }
	fscanf(fptr,".\n");
        for(n=1; n <= N; n++) {
           fscanf(fptr,"%d(%d aa): ",&n,&len);
	   for(i=1; i<= ntyp; i++){
               fscanf(fptr,"(%d)-{%d}[%f]-",&g,&t,&p);
	       /** if(g <=  F->maxLength) F->scoreGap[i-1][g]++;  /*** OLD ***/
	       if(g <=  F->maxLength && p > 1.3) F->scoreGap[i-1][g]++;  /*** NEW ***/
	   }
           fscanf(fptr,"(%d).\n",&g);
	   if(g <=  F->maxLength) F->scoreGap[ntyp][g]++; 
        }
	return TRUE;
}

void	PermuteOScan(Int4 permute, osn_typ F)
/*************** Circularly permute motifs ****************/
{
	Int4	m,m2,n,i;
	wm_type	*M;

   if(permute > 0){
   	NEW(M,F->N+2,wm_type);
	i = permute; n = F->N;
	if(i >= n) i = i%n;
	fprintf(stderr,"permute = %d\n",i);
	for(m=1; m<= n; m++){
		m2 = (m+i)%n;
		if(m2 == 0) m2 = n;
		M[m2] = F->M[m];
		fprintf(stderr,"motif %c -> %c\n",
			m + 'A' -1, m2 + 'A' -1);
	}
	for(m=1; m<= n; m++){ F->M[m] = M[m]; } free(M[m]);
   } else fprintf(stderr,"permute = %d\n",permute);
}

void	PutSMXOScan(FILE *fptr,osn_typ F)
{
	smx_typ	smx;
	Int4	m;

	for(m=1; m<=F->N; m++) { 
	    smx = GetSMatrixWModel(F->M[m]);
	    PutSMatrix(fptr, smx);
	}
}

snh_typ	OScanScan(FILE *fptr,Int4 number, unsigned short *nsize, osn_typ F)
/*** impose colinearity constraint and return heap with results. ***/
{
	e_type	E;
	wm_type	*M=F->M;
	a_type	A=F->A;
	unsigned char	*seq;
	char	c;
	Int4	best,s,i,j,end,*w,tot_len;
	Int4	score,gapscore,score2,sum_score,last_sum_score;
	Int4	field,hpsz,N,n,r,rmax,N0,m,m0,*p,*p0,rpts;
	double	pval,k,factor,prob,eval,target;
	float	*pv,singleEval;
	BooLean	okay;
	smx_typ	smx,*sM,sm;
	snh_typ	sH;
	sni_typ	I,I0,bestI,List=NULL;
	Int4	gap=0,ungap=0;
	Int4    **scoreGap;

/**************** SPOUGE I ************************/
	Int4	max_gap_score;
	// js_type	Spg;
	h_type	H[2];

	H[0] = Histogram("scores for sequences discarded", 0,800,10); 
	H[1] = Histogram("scores for sequences checked by Spouge functions",
			0,800,10); 
/**************** SPOUGE I ************************/
	singleEval = (float) -log10(F->singleEval);


	/**** ALLOCATE MEMORY & OBJECTS ****/
	rpts = F->maxrpts;
	N0 = F->N; N = rpts*N0; 
	NEW(p,N+3,Int4); NEW(pv,N+3,float); NEW(sM,N+3,smx_typ); 
	NEW(w,N0+2,Int4); 
	for(w[0]=0, tot_len=0,m=1; m<=N0; m++) { 
		w[m] = LenWModel(M[m]); tot_len += w[m];
	}
	smx = MkSMatrix(2.0,tot_len,F->freq,A);

	/**** CONSTRUCT FULL SCORING MATRIX. ****/
	for(i=1,m=1; m<=N0; m++) { 
	    for(j=1; j<=w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
			score = CellScoreWModel(c, j, M[m]); 
			SetSMatrix(c,i,score,smx);
		}
	    }
	    sM[m] = GetSMatrixWModel(M[m]);
	    for(r=1; r < rpts; r++){ j = m + r*N0; sM[j] = sM[m]; }
	}
	/*** PutSMatrix(stderr,smx); /******/
	/******/
	fprintf(stderr,"min = %d; max = %d --> p = %g; p = %g\n",
		MinScoreSMatrix(smx), MaxScoreSMatrix(smx),
	   	SMatrixProb(MinScoreSMatrix(smx), smx),
	   	SMatrixProb(MaxScoreSMatrix(smx), smx));
	/******/

/**************** SPOUGE II ************************/
  if(F->mode != 'G'){
	print_error("Spouge gap functions deactivated");
#if 0
        for(n=0,s=1; s<=number;s++) { n = MAXIMUM(Int4,n,nsize[s]); }
	Spg = MkSpouge(nAlpha(A)+1, F->freq, N0, sM, n+1, F->scoreGap,F->mode);
	NEWP(scoreGap,N+3,Int4); 
	scoreGap[0] = GapScoreSpouge1(0,Spg);
	for(m=1; m<=N0; m++) { 
	    scoreGap[m] = GapScoreSpouge1(m,Spg);
	    for(r=1; r < rpts; r++){ j = m + r*N0; scoreGap[j] = GapScoreSpouge1(m,Spg); }
	}
	max_gap_score = MaxGapScoreSpouge(Spg);
	fprintf(stderr,"n = %d; max_gap_score = %d\n",n,max_gap_score);
#endif
  }
/**************** SPOUGE II ************************/

	/**** SEARCH DATABASE FOR MATCHES. ****/
        for(s=1; s<=number;s++) {

	  /**** SCAN NEXT SEQUENCE FOR MOTIFS ****/ 
	  E = ReadSeq(fptr,s,nsize[s],A);
	  if(F->shuffle) ShuffleSeq(E);
	  if(tot_len <= (Int4) LenSeq(E)){
	   seq = SeqPtr(E);
	   target = F->maxEval/((double)F->total/(double) LenSeq(E));
	   n = LenSeq(E) - tot_len; m = N0;
	   /*** Altschul's method to adjust for multiple blocks ***/
	   factor = bico(n+m,m)*F->total/(double) LenSeq(E); /***/
	   /***  factor = bico(n+m,m)*(double) number; /** OLD **/

	   if(F->mode == 'D'){
		print_error("Spouge gap functions deactivated");
#if 0
	   	score=LocalAlnSeqSMatrix(LenSeq(E), SeqPtr(E), N0, sM,p);
		pval = log10(SpougePvalue(Spg,F->total,LenSeq(E),score)); 
#endif
	   } else {
	   	score=AlnSeqSMatrix(LenSeq(E), SeqPtr(E), N0, sM,p);
	   	pval = SMatrixProb(score, smx);
	   	pval = log10(pval*factor); 
	   }

	   if(pval > F->log10ecut) {
            if(F->gapfunct){
	     print_error("Spouge gap functions deactivated");
#if 0
	     if((score + max_gap_score) < SpougeThreshold(Spg,target,LenSeq(E))){
	        if(F->dfp != NULL) PutSeq(F->dfp,E,A); 
		NilSeq(E); E = NULL; IncdHist(score,H[0]);
	     } else { ProcessSeqPSeg(14,2.2,2.5,100,E,A); IncdHist(score,H[1]); }
#endif
	    } else { if(F->dfp != NULL) PutSeq(F->dfp,E,A); NilSeq(E); E = NULL; }
	   } else {	
	     if(F->recomb!= NULL) PutSeq(F->recomb,E,A); /** circular permuted **/
	     /**** IF P VALUE LOOKS PROMISING THEN DO RIGOROUS CHECK ****/ 
	     ProcessSeqPSeg(14,2.2,2.5,100,E,A);

	     /*** DO DYNAMIC PROGRAMMING ALIGNMENT WITH MOTIFS ***/
	     rmax = (Int4)LenSeq(E)/tot_len; rmax = MINIMUM(Int4,rmax,rpts);
	     I = NULL; 
	     if(F->mode == 'D'){
	   	score=LocalAlnSeqSMatrix(LenSeq(E), XSeqPtr(E), N0, sM,p);
	     } else {
	        score=AlnSeqSMatrix(LenSeq(E), XSeqPtr(E), N0, sM,p);
	     }

if(F->gaps){	/*** BEGIN TEST: FULLY GAPPED ALIGNMENT ****/
  prob = SMatrixProb(score, smx); 
  pval = log10(prob*factor); 
  if(TRUE || pval <= (3.0 + F->log10ecut)){
	  score2=AlnSeqSMatrixSW(15,10,LenSeq(E),XSeqPtr(E),N0,sM,p,score);
#if 0
	  score2=AlnSeqSMatrixSW(10,8,LenSeq(E),XSeqPtr(E),N0,sM,p,score);
#endif
	  if(score2 > score){
	     gap++;
	     PutSeqInfo(stdout,E);
	     eval = factor*SMatrixProb(score2, smx);
	     if(eval < 1.0) 
		printf("\n !gapped score = %d (p = %g)\n",score2,eval);
	     else printf("\n gapped score = %d (p = %g)\n",score2,eval);
	     printf("\n ungapped score = %d (p = %g)\n",score,
			factor*SMatrixProb(score, smx));
	     printf("********************************\n\n");
	  } else ungap++;
	  if(score2 > score){ PutSeq(F->fp_gaps,E,A); }
  }
} /*** END TEST: FULLY GAPPED ALIGNMENT ****/

	     p[0]=1; p[F->N+1]=LenSeq(E)+1;

	     /*** Calculate p value for ordered motifs ***/
	     if(F->mode == 'D'){
	        print_error("Spouge gap functions deactivated");
		// pval = log10(SpougePvalue(Spg,F->total,LenSeq(E),score)); 
	     } else {
		prob = SMatrixProb(score, smx); 
		pval = log10(prob*factor); 
	     }
	     if(pval <= F->log10ecut){
		for(m = 1; m<= F->N; m++){
		  field = p[m+1] - (p[m-1] + w[m-1]) - w[m] + 1;
                  prob = PvalWModel(SeqPtr(E),p[m],F->M[m]);
		  pv[m] = (float)-log10(prob*(double)field);
/**********
fprintf(stderr,"field = %d; pv[m]=%g;prob=%g\n",field,pv[m],prob);
/**********/
		}
		UnXSeq(E); /** remove X'ed out regions **/
		I=MakeScanInfo(E,pval,N0,1,pv,p,s,'G');

	       /****************** CHECK FOR REPEATS *********************
                For gapped Spouge routines need to make repeats gapscore structure, 
                use spouge threshold score instead of SMatrixProb( ), and use 
                GapFuncAlnSeqSMatrix( ) instead of AlnSeqSMatrix( ).  This should
                be all that is needed 
               ***********************************************************/ 
               if(rpts > 1){
                 if(F->mode != 'G' && F->mode != 'g') 
            		print_error("OScanScan( ) input error 'G'");
            	 N0 = F->N; n = LenSeq(E); 
            	 factor = bico((n + N0 - tot_len),N0);  /** single model & seq **/
            	 rmax = n/tot_len; rmax = MINIMUM(Int4,rmax,rpts);
            	 for(last_sum_score = score, r=2; r <= rmax; r++){
            	    N = r*N0;
            	    p[0]=1; p[N+1]=n+1;	 /** set ends of sequences **/
              	    sum_score=AlnSeqSMatrix(LenSeq(E), XSeqPtr(E), N, sM,p);
            	    /***********
            	      if(r==2) PutSeqInfo(stderr,E); 
            	      fprintf(stderr,
            	        "r = %d; score = %d; sum_score=%d; sum_score-last_sum=%d\n", 
            		r,score,sum_score,sum_score - last_sum_score);
            	    /***********/
            	    prob = SMatrixProb(sum_score - last_sum_score, smx);
            	    last_sum_score = sum_score;
            	    prob = log10(prob*factor);
            	    if(prob <= F->log10ecut) {
                         for(m = 0; m < N; ){
            		      m0 = (m % N0) + 1; m++;
            		      field = p[m+1] - (p[m-1] + w[m0-1]) - w[m0] + 1;
                              prob = PvalWModel(SeqPtr(E),p[m],F->M[m0]);
                              pv[m] = (float)-log10(prob*(double)field);
                         }
            	         I = ReMakeScanInfo(I,pval,N,r,pv,p);
            	    } else break; 
            	 }
               } /********* END CHECK FOR REPEATS *********************/ 

	        List = AppendScanInfo(I, List);
		E = NULL;
	     } 
	   } /*** END DYNAMIC PROGRAMMING **/
/**************** SPOUGE IId ************************/
   if(E != NULL && F->gapfunct){
	print_error("Spouge gap functions deactivated");
#if 0
	switch (F->mode) {
	 case 'g': gapscore=GapFuncAlnSeqSMatrix(LenSeq(E), XSeqPtr(E), N0,
			sM,p, GapScoreSpouge(Spg)); break;
	 default: print_error("input error in OScanScan");
	 break;
	}
	if(gapscore >= SpougeThreshold(Spg,F->total,LenSeq(E))){
	          p[0]=1; p[F->N+1]=LenSeq(E)+1;
		  for(m = 1; m<= F->N; m++){
		      field = p[m+1] - (p[m-1] + w[m-1]) - w[m] + 1;
                      prob = PvalWModel(SeqPtr(E),p[m],F->M[m]);
		      pv[m] = (float)-log10(prob*(double)field);
		  }
		  UnXSeq(E); /** remove X'ed out regions **/
		  I=MakeScanInfo(E,log10(SpougePvalue(Spg,F->total,LenSeq(E),gapscore)),
				N0,1,pv,p,s,F->mode);
	          List = AppendScanInfo(I, List);
#if 0 /*** DEBUG ***/
	          score=AlnSeqSMatrix(LenSeq(E), XSeqPtr(E), N0, sM,p);
		  prob = SMatrixProb(score, smx); 
		  pval = prob*factor; 
		  fprintf(stderr,"\n************************\n");
		  PutSeqInfo(stderr,E);
		  fprintf(stderr,"\n (%d residues)\n",LenSeq(E));
	   	  fprintf(stderr,"score=%d,Eval=%f\n",score,pval);
		  fprintf(stderr,"gapped score=%d; Eval = %g\n",
			gapscore,SpougePvalue(Spg,F->total,LenSeq(E),gapscore)); 
	  	  fprintf(stderr,"spouge threshold score = %d\n",
			SpougeThreshold(Spg,F->total,LenSeq(E)));
#endif /*** DEBUG ***/
		  E = NULL; /** don't destroy; saved in I **/
	}
#endif
   }
/**************** SPOUGE IId ************************/

	   } /** end if (model length <= seq_length) ***/
	   if(E!= NULL) {
	      if(F->dfp != NULL) PutSeq(F->dfp,E,A); 
	      NilSeq(E); E = NULL;
	   }
	}
	if(F->gaps) printf("gapped = %d; ungapped = %d; expected score = %g\n",
		gap,ungap,ExpScoreSMatrix(smx)); 
	hpsz = LengScanInfo(List);
	/** printf("hpsz = %d\n",hpsz); /***/
	sH = MakeScanHeap(hpsz+1,F->N,w,F->singleEval);
	if(List != NULL){
	  while((I = RmLastScanInfo(List)) != List){
		if(InsertScanHeap(I, sH)==NULL) NilScanInfo(I);
	  }
	  if(InsertScanHeap(I, sH)==NULL) NilScanInfo(I);
	}
/**************** SPOUGE III ************************/
   if(F->mode != 'G'){ 
   	PutHist(stderr,60,H[0]); 
   	PutHist(stderr,60,H[1]); 
	free(scoreGap); 
	// NilSpouge(Spg); 
   }
   NilHist(H[0]); NilHist(H[1]);
/**************** SPOUGE III ************************/
	
	free(pv); free(p); free(w); free(sM); NilSMatrix(smx);
	return sH;
}

snh_typ	OScanScan1(FILE *fptr,Int4 number, unsigned short *nsize, osn_typ F,
	BooLean combine)
/*** impose colinearity constraint and return heap with results. ***/
{
	e_type	E;
	wm_type	*M=F->M,*M2=F->M2;
	a_type	A=F->A;
	unsigned char	*seq;
	char	c,mode;
	Int4	s,i,j,start,end,*w,tot_len,len;
	Int4	score,score1,score2,sum_score,last_sum_score;
	Int4	field,hpsz,N,n,r,rmax,N0,m,m0,*p,*p0,rpts,repeats;
	Int4	gap=0,ungap=0;
	double	pval,pval2,k,factor,prob,eval,mask_drop=1.0,target;
	float	*pv,singleEval;
	BooLean	masked,okay;
	smx_typ	smx,*sM,sm,*fsM;
	snh_typ	sH;
	sni_typ	I,I0,List=NULL;
	Int4    **scoreGap;
	Int4	max_gap_score=0;
	// js_type	Spg = NULL;
	h_type	H[2];
	/*** PDB ***/
	FILE	*fp;
	Int4	*pdb=NULL;
	/*** PDB ***/

	H[0] = Histogram("scores for sequences discarded", 0,800,10); 
	H[1] = Histogram("scores for sequences checked by Spouge functions",
			0,800,10); 
	singleEval = (float) -log10(F->singleEval);

	/**** ALLOCATE MEMORY & OBJECTS ****/
	rpts = F->maxrpts;
	N0 = F->N; N = rpts*N0; 
	NEW(p,N+3,Int4); NEW(pv,N+3,float); 
	NEW(sM,N+3,smx_typ); NEW(fsM,N+3,smx_typ); 
	NEW(w,N0+2,Int4); 
	for(w[0]=0, tot_len=0,m=1; m<=N0; m++) { 
		w[m] = LenWModel(M[m]); tot_len += w[m];
	}
	/**** CONSTRUCT FULL SCORING MATRIX. ****/
	smx = MkSMatrix(2.0,tot_len,F->freq,A);
	for(i=1,m=1; m<=N0; m++) { 
	    for(j=1; j<=w[m]; j++,i++) { 
		for(c=0; c<=nAlpha(A); c++) { 
			score = CellScoreWModel(c, j, M[m]); 
			SetSMatrix(c,i,score,smx);
		}
	    }
	    sM[m] = GetSMatrixWModel(M[m]);
	    if(M2[m] == NULL) fsM[m] = sM[m];
	    else fsM[m] = GetSMatrixWModel(M2[m]);
	    for(r=1; r < rpts; r++){ j=m+r*N0; sM[j]=sM[m]; fsM[j]=fsM[m]; }
	}
	if(F->mode != 'G'){
	   print_error("Spouge gap functions deactivated");
#if 0
           for(n=0,s=1; s<=number;s++) {
		n = MAXIMUM(Int4,n,nsize[s]);
	   }
	   Spg = MkSpouge(nAlpha(A)+1,F->freq,N0,sM,n+1,F->scoreGap,F->mode);
	   if(islower(F->mode)){
	   	NEWP(scoreGap,N+3,Int4); 
	   	scoreGap[0] = GapScoreSpouge1(0,Spg);
	   	for(m=1; m<=N0; m++) { 
	      	   scoreGap[m] = GapScoreSpouge1(m,Spg);
	      	   for(r=1; r < rpts; r++){ 
			j = m + r*N0; 
			scoreGap[j] = GapScoreSpouge1(m,Spg); 
		   }
		}
	   } else scoreGap=NULL;
	   max_gap_score = MaxGapScoreSpouge(Spg);
	   fprintf(stderr,"n = %d; max_gap_score = %d\n",n,max_gap_score);
#endif
	}
/**************** PDB ************************/
if(F->pdbfa != NULL){
	if(rpts != 1) print_error("OScanScan( ) pdb hom option input error");
        E = ReadSeqFA(F->pdbfa, 1, F->A);
	NEW(pdb,N0+3,Int4); 
	pdb[0]=1; pdb[N0+1]=LenSeq(E)+1; /** set ends **/
        AlnSeqSMatrix(LenSeq(E), SeqPtr(E),N0,sM,pdb);
	fp = open_file(F->pdbfa,".hom","w");
	PutSeqMtfMask(fp,N0, w, pdb, pdb,E,A);
	NilSeq(E);
} else pdb=NULL;
/**************** PDB ************************/

	/**** SEARCH DATABASE FOR MATCHES. ****/
        for(s=1; s<=number;s++) {

	  /**** SCAN NEXT SEQUENCE FOR MOTIFS ****/ 
	  mode = F->mode;
	  E = ReadSeq(fptr,s,nsize[s],A); len = LenSeq(E);
	  if(F->shuffle) ShuffleSeq(E);
	  if(tot_len <= len){
	   seq = SeqPtr(E);
	   n = len - tot_len; m = N0;
	   target = F->maxEval/((double)F->total/(double) len);
	   factor = bico(n+m,m)*F->total/(double) len; /* Altschul's way */
	   /***  factor = bico(n+m,m)*(double) number; /** OLD **/

	   /**** DO AN INITIAL CHECK ****/ 
	   switch (F->mode) {
	     case 'D': case 'd':
	   	score=LocalAlnSeqSMatrix(len,SeqPtr(E),N0,sM,NULL);
		break;
	     case 'G': case 'g':
	   	score=AlnSeqSMatrix(len,SeqPtr(E),N0,sM,NULL);
		break;
	     case 'O': case 'o':
	     	print_error("OScanScan( ) mode not yet implemented");
		break;
	     default: print_error("OScanScan( ) mode input error");
	   } 
	   if(F->mode == 'G'){
	   	pval = log10(SMatrixProb(score,smx)*factor); 
		if(pval > F->log10ecut) okay = FALSE;
		else okay = TRUE;
	   } else {
	        print_error("Spouge gap functions deactivated");
#if 0
	        if((score + max_gap_score) < SpougeThreshold(Spg,target,len)){
		     if(combine && F->mode == 'g'){
	   		pval = log10(SMatrixProb(score,smx)*factor); 
			if(pval > F->log10ecut) okay = FALSE;
			else okay = TRUE;
		     } else okay = FALSE; 
		     IncdHist(score,H[0]); 
	        } else { okay = TRUE; IncdHist(score,H[1]); }
#endif
	   } 

	   /**** IF INITIAL CHECK LOOKS PROMISING THEN DO RIGOROUS CHECK ****/ 
	   if(okay) {	
	     rmax = len/tot_len; rmax = MINIMUM(Int4,rmax,rpts);
	     I = NULL; 
	     if(F->segmask) {
	        masked=ProcessSeqPSeg(14,2.2,2.5,100,E,A); /* compositional bias? */
	     } else {
		masked=FALSE;
	     }
	     switch (F->mode) {
	        case 'D':   /** score already computed **/
	            print_error("Spouge gap functions deactivated");
#if 0
		    pval = log10(SpougePvalue(Spg,F->total,len,score)); 
		    if(masked){
			score2=LocalAlnSeqSMatrix(len,XSeqPtr(E),N0,sM,NULL);
		        pval2 = log10(SpougePvalue(Spg,F->total,len,score2)); 
		    } else pval2=pval; 
#endif
		    break;
		case 'd': 
	            print_error("Spouge gap functions deactivated");
#if 0
	            score=LocalGapFuncAlnSMatrix(len,SeqPtr(E),
				N0,sM,NULL,GapScoreSpouge(Spg));
		    pval = log10(SpougePvalue(Spg,F->total,len,score)); 
		    if(masked){
	                score2=LocalGapFuncAlnSMatrix(len,XSeqPtr(E),
				N0, sM,NULL,GapScoreSpouge(Spg));
		        pval2 = log10(SpougePvalue(Spg,F->total,len,score2)); 
		    } else pval2=pval; 
#endif
		    break;
	        case 'G':   /** score and positions already computed **/
		    pval = log10(SMatrixProb(score, smx)*factor); 
		    if(masked){
		       score2=AlnSeqSMatrix(len,XSeqPtr(E),N0,sM,NULL); 
		       pval2 = log10(SMatrixProb(score2, smx)*factor); 
		    } else pval2=pval; 
		    break;
		case 'g':
	            print_error("Spouge gap functions deactivated");
#if 0
		    if(combine) prob=log10(SMatrixProb(score, smx)*factor); 
	            score=GapFuncAlnSeqSMatrix(len, SeqPtr(E),
				N0, sM,NULL,GapScoreSpouge(Spg));
		    pval = log10(SpougePvalue(Spg,F->total,len,score)); 
		    if(combine && pval < SpougeThreshold(Spg,F->total,len) && prob < pval){
			score=AlnSeqSMatrix(len,SeqPtr(E),N0,sM,NULL);
			pval=prob; mode = 'G';
		    } 
		    if(masked){
	               score2=GapFuncAlnSeqSMatrix(len, XSeqPtr(E),
				N0, sM,NULL,GapScoreSpouge(Spg));
		       pval2 = log10(SpougePvalue(Spg,F->total,len,score2)); 
		       if(combine) {
			   score2=AlnSeqSMatrix(len,XSeqPtr(E),N0,sM,NULL); 
			   prob = log10(SMatrixProb(score2, smx)*factor); 
			   pval2 = MINIMUM(double,prob,pval2); 
		       }
		    } else pval2=pval; 
#endif
		    break;
	        default: print_error("OScanScan( ) mode input error");
	     }
	     if(pval <= F->log10ecut && (pval2 <= (F->log10ecut+mask_drop))){ 

               /****************** CHECK FOR REPEATS *********************/
               N0 = F->N; repeats=1;
               if(rpts > 1){
                 rmax = MINIMUM(Int4,len/tot_len,rpts);
                 factor = bico((len + N0 - tot_len),N0);  /** single model & seq **/
                 for(last_sum_score = score, r=2; r <= rmax; r++){
                    N = r*N0;
	            switch (F->mode) {
	        	case 'D': 
	            	  print_error("Spouge gap functions deactivated");
#if 0
	   	    	  sum_score=LocalAlnSeqSMatrix(len,SeqPtr(E),N,sM,NULL);
		          prob=SpougePvalue(Spg,F->total,len,sum_score-last_sum_score); 
#endif
			  break;
	        	case 'd': 
	            	  print_error("Spouge gap functions deactivated");
#if 0
			  sum_score=LocalGapFuncAlnSMatrix(len,SeqPtr(E),N,sM,
						NULL,scoreGap);
		          prob=SpougePvalue(Spg,F->total,len,sum_score-last_sum_score); 
#endif
			  break;
	        	case 'G': 
                    	  sum_score=AlnSeqSMatrix(len,SeqPtr(E),N,sM,NULL);
                          prob=factor*SMatrixProb(sum_score-last_sum_score,smx);
			  break;
	        	case 'g': 
	            	  print_error("Spouge gap functions deactivated");
#if 0
	            	  sum_score=GapFuncAlnSeqSMatrix(len,SeqPtr(E),N,sM,
						NULL,scoreGap);
		          prob=SpougePvalue(Spg,F->total,len,sum_score-last_sum_score); 
#endif
			  break;
	        	default: print_error("OScanScan( ) mode input error");
		    }
                    last_sum_score = sum_score;
                    if(log10(prob) <= F->log10ecut) repeats = r;
                    else break;
                 }
               }
	       N = repeats*N0;
	       UnXSeq(E); /** remove X'ed out regions **/
/*** THIS IS STILL NOT WORKING CORRECTLY; I'D LIKE TO COMPUTE Spouge BOTH WITH
     AND WITHOUT IGNORED COLUMNS BUT DON"T WANT TO USE MORE COMPUTE TIME.  ASK 
     SPOUGE IF CAN BE DONE AS A SIDE JOB (I WOULD THINK IT CERTAINLY COULD). ***/
	       switch (F->mode) {
	          case 'D': LocalAlnSeqSMatrix(len,SeqPtr(E),N,fsM,p); break;
	          case 'd': LocalGapFuncAlnSMatrix(len,SeqPtr(E),N,fsM,p,scoreGap);
			    break;
	          case 'G': AlnSeqSMatrix(len,SeqPtr(E),N,fsM,p); break;
	          case 'g': GapFuncAlnSeqSMatrix(len,SeqPtr(E),N,fsM,p,scoreGap);
			    break;
	          default: print_error("OScanScan( ) mode input error");
	       }
	       p[0]=1; p[N+1]=len+1; /** set ends **/
/**************** PDB ************************/
if(pdb != NULL){ PutSeqMtfMask(fp,N0, w, pdb, p,E,A); }
/**************** PDB ************************/
               for(m = 0; m < N; ){
		   m0 = (m % N0) + 1; m++;
                   field = p[m+1] - (p[m-1] + w[m0-1]) - w[m0] + 1;
                   prob = PvalWModel(SeqPtr(E),p[m],F->M[m0]);
                   pv[m] = (float)-log10(prob*(double)field);
               }
		I=MakeScanInfo(E,pval,N,repeats,pv,p,s,mode);
	        List = AppendScanInfo(I, List);
	        E = NULL;
	     } /** end of if(pval... **/
	   } /*** end if(okay)... **/
	  } /** end if (model length <= seq_length) ***/
	  if(E!= NULL) {  /** Sequences that failed scan **/
	      if(F->dfp != NULL) PutSeq(F->dfp,E,A); 
	      NilSeq(E); E = NULL;
	  }
	}
	hpsz = LengScanInfo(List);
	sH = MakeScanHeap(hpsz+1,F->N,w,F->singleEval);
	if(List != NULL){
	  while((I = RmLastScanInfo(List)) != List){
		if(InsertScanHeap(I, sH)==NULL) NilScanInfo(I);
	  }
	  if(InsertScanHeap(I, sH)==NULL) NilScanInfo(I);
	}

/**************** PDB ************************/
if(F->pdbfa != NULL){ free(pdb); fclose(fp); }
/**************** PDB ************************/
/**************** TEST ************************/
if(F->mode != 'G' && F->mode != 'D'){
	PutHist(stderr,60,H[0]); 
	PutHist(stderr,60,H[1]); 
}
	NilHist(H[0]); NilHist(H[1]);
/**************** TEST ************************/
	
        // if(Spg != NULL) { NilSpouge(Spg); free(scoreGap); }
	free(pv); free(p); free(w); free(sM); free(fsM); NilSMatrix(smx);
	return sH;
}

