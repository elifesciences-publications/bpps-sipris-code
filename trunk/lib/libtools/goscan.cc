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

#include "goscan.h"

gsn_typ MakeGOScan(char *snfile, a_type A, Int4 maxrpts, char method, char mode, 
	double ecut,double Ecut,Int4 maxLength, BooLean weight,float minmap,
	double pseudo,double *freq)
{
	gsn_typ F;

	NEW(F,1,goscan_type);
	F->pssm=0; F->pssm_cma=0;
	F->num_models = 1;
	NEW(F->pm,F->num_models +3, ptm_typ);
	F->shuffle=FALSE;
	F->recomb = NULL;
	F->snfile = AllocString(snfile);
	F->dfp = NULL;
	if(ecut <= 0.0 || Ecut <= 0.0)
		print_error("maximum E-value must be >= 0.0");
	F->maxEval = ecut; F->log10ecut = log10(ecut);
	F->singleEval = Ecut;
	F->repeatEval = 0.01;
	F->PM = MakePrtnModel(snfile, A, maxrpts, method, mode,
        	maxLength, weight,minmap,pseudo,0.0,freq);
	if(!F->PM){
	    assert(F->PM);
	}
	F->pm[1] = F->PM;
	return F;
}

gsn_typ MakeMGOScan(char *snfile, a_type A, Int4 maxrpts, char method, char mode, 
	double ecut,double Ecut,Int4 maxLength, BooLean weight,float minmap,
	double pseudo,double *freq)
{
	gsn_typ F;
	Int4	i;
	ptm_typ	tmp_pm[1002];

	NEW(F,1,goscan_type);
	F->pssm=0; F->pssm_cma=0;
	FILE *fp=open_file(snfile,"","r");
	i=0;
	do {
		i++;
		if(i > 1000) print_error("MakeMGOScan(): too many input models");
		tmp_pm[i]=MakePrtnModel(fp,A,maxrpts,method,mode,
				maxLength,weight,minmap,pseudo,0.0,freq);
	} while(tmp_pm[i]);
	// for(i=1; snfile[i] != NULL; ) i++; 
	F->num_models = i-1;
	NEW(F->pm,F->num_models +3, ptm_typ);
	F->shuffle=FALSE;
	F->recomb = NULL;
	F->snfile = AllocString(snfile);
	F->dfp = NULL;
	if(ecut <= 0.0 || Ecut <= 0.0)
		print_error("maximum E-value must be >= 0.0");
	F->maxEval = ecut; F->log10ecut = log10(ecut);
	F->singleEval = Ecut;
	F->repeatEval = 0.01;
#if 0	// old version...
	F->PM = MakePrtnModel(snfile[1], A, maxrpts, method, mode,
        	maxLength, weight,minmap,pseudo,0.0,freq);
	F->pm[1] = F->PM;
	for(i=2; i<=F->num_models; i++){
	    F->pm[i] = MakePrtnModel(snfile[i], A, maxrpts, method, mode,
        	maxLength, weight,minmap,pseudo,0.0,freq);
	}
#else
	F->PM = tmp_pm[1];
	for(i=1; i <= F->num_models; i++) F->pm[i] = tmp_pm[i];
#endif
	return F;
}

char	*DescriptionGOScan(Int4 i,gsn_typ F)
{
	if(i > 0 && i <= F->num_models) return DescriptionsmaPrtnModel(F->pm[i]);
	else return 0;
}

void    PutSmxGOScan(FILE *fptr,gsn_typ F)
{
	Int4	i,N = NumModelsPrtnModel(F->PM);
	
	for(i=1; i<=N; i++){
		fprintf(fptr,"============= SMatrix %d ===========\n",i);
		PutSmxPrtnModel(i, F->PM);
		fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n");
}

void	NoMaskGOScan(gsn_typ F) { NoMaskPrtnModel(F->PM); }

void	MaskNonGlobularGOScan(gsn_typ F) { MaskNonGlobularPrtnModel(F->PM); }

void	SetEvalGOScan(double ecut, double Ecut, gsn_typ F)
{
	if(ecut <= 0.0 || Ecut <= 0.0)
        	print_error("maximum E-value must be >= 0.0");
        F->maxEval = ecut; F->log10ecut = log10(ecut);
	F->singleEval = Ecut;
}

FILE	*OpenDatabaseGOScan(gsn_typ F)
{ F->dfp=open_file(F->snfile,".dbs","w"); return F->dfp; }

void	NilGOScan(gsn_typ F)
{
	free(F->snfile);
	if(F->dfp != NULL) fclose(F->dfp);
	if(F->recomb != NULL) fclose(F->recomb);
	NilPrtnModel(F->PM);
	if(F->num_models > 1){
	   for(Int4 i=2; i<= F->num_models; i++) NilPrtnModel(F->pm[i]);
	}
	free(F->pm);
	if(F->pssm) delete F->pssm;
	if(F->pssm_cma) TotalNilCMSA(F->pssm_cma);
	free(F);
}

Int4	gap_seq_goscan(char *operation, Int4 start0, gsn_typ F, Int4 *left, Int4 *right)
/****************************************************************************
  create 'fake' subsequence with gaps and flanking 
  trace = "EDdmmmmmmmmmmmiiiMmmmmmmmmmmmmmmmmmmmmiiMmmmmmmmmmmE"
 ****************************************************************************/
{
        Int4            o,j,t,n,offsetX=0,start=0,end=0;
        char            state;

        assert(operation[0]=='E'); 
        assert(operation[strlen(operation)-1]=='E');

	Int4 numM=NumModelsPrtnModel(PrtnModelGOScan(F));
        // 2. Add aligned-region sequence.
        for(j=start0,end=t=0,n=o=1,state='E'; operation[o] != 'E'; o++){
            switch(operation[o]){
               case 'M': 
		  if(t==0) start=j; 
		  if(t == numM){ // end of profile...
			if(state != 'i') end =j-1; 
			left[n] = start; right[n] = end; n++;
			start =j; t=1;
		  } else t++;
               case 'm': j++; break;
               case 'D': // starting deletion in sequence relative to profile.
		  if(t==0) start=j; 
		  if(t == numM){ // end of sequence...
			if(state != 'i') end =j-1; 
			left[n] = start; right[n] = end; n++;
			start=j; t=1;
		  } else t++;
               case 'd': // internal deletion in sequence relative to profile.
		  break;
               case 'i': // insert is between profile blocks;
                  assert(state != 'i'); 
		  if(t == numM) end=j-1; 
                  while(operation[o]=='i'){ j++; o++; } o--; 
		  break;
               case 'I': // Insert within a profile block; delete from seq.
                  assert(state != 'I'); 
                  while(operation[o]=='I'){ j++; o++; } o--; 
		  break;
               default:
                        fprintf(stderr,"operations = %s\n",operation);
                        print_error("gap_seq_goscan( ): input error"); break;
            }  state=operation[o];
        }
	if(state != 'i') end =j-1;
	left[n] = start; right[n] = end;
// fprintf(stderr,"DEBUG1: s=%d; e=%d; %d ; %d\n",s,e,start,end);
	return n;
}

BooLean contain_segs_goscan(Int4 *s, Int4 *e, Int4 R, Int4 *gs, 
	Int4 *ge, Int4 gR,double cutoff)
// see if blocked segment set is contained within gapped set
/**********************************************************************

        if(gs[g] <= e[r] && s[r] <= ge[g])
        TRUE:                          FALSE:
               s         e                           s        e
        =======-----------=========    ==============----------=====
        ===========-----------=====    ====---------=========
                   gs        ge            gs      ge
 **********************************************************************/
{
	BooLean	found;
	Int4	r,g,S,E,minlen,us,ue;
	double	overlap;

	assert(cutoff > 0.0 && cutoff <= 1.0);
	for(r=1; r <= R; r++){
	   us=s[r]; ue=e[r];
	   for(found=FALSE,g=1; g <= gR; g++){
                if(gs[g] <= ue && us <= ge[g]){  // i.e., some overlap
		    S = MAXIMUM(Int4,us,gs[g]);
		    E = MINIMUM(Int4,ue,ge[g]);
		    minlen = MINIMUM(Int4,ue-us+1,ge[g]-gs[g]+1);
		    overlap = (double)(E-S+1)/(double)minlen;
		    if(overlap >= cutoff){ found = TRUE; break; }
		}
	   } if(!found) return FALSE;
	} return TRUE;
}

e_type	**GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 left, Int4 right,
	Int4 min_rpt, Int4 number,UInt8 total, unsigned short *nsize, 
	gsn_typ F)
{ return GapSeqGOScanScan(stderr,fptr,a,b,left,right,min_rpt,number,total,nsize,F); }

e_type	**GapSeqGOScanScan(FILE *ofp, FILE *fptr,Int4 a, Int4 b, Int4 left, Int4 right,
	Int4 min_rpt, Int4 number,UInt8 total, unsigned short *nsize, gsn_typ F)
// find optimum gapped sequence repeats 
{
	e_type	E;
	a_type	A;
	Int4	r,s,st,e,totN,*p,repeats,N;
	char	mode;
	double	pval;
	ptm_typ PM=PrtnModelGOScan(F);
	char	*operation,*best_operation;
	Int4	start,score,lastscore,this_score;
	float	*pv;
	Int4	*gstrt,*gend,max_rpts,best_rpts;

	assert(!F->shuffle); assert(min_rpt > 0);
	totN=NmaxPrtnModel(F->PM); N=NumModelsPrtnModel(F->PM);
	NEW(p,totN+3,Int4); NEW(pv,totN+3,float); 
	mode = ModePrtnModel(F->PM); A=PrtnModelA(F->PM);
	max_rpts=MaxRptsPrtnModel(F->PM);
	NEW(gstrt,max_rpts+2,Int4); NEW(gend,max_rpts+2,Int4);
	Int4 *lenM=PrtnModelLengths(F->PM);

	e_type	**ListE;
	Int4	hits;
	NEWP(ListE, number +2, e_type);
        for(hits=0, s=1; s<=number;s++) {
	    E = ReadSeq(fptr,s,nsize[s],A); assert(E);
	    repeats=ComparePrtnModelRpts(E,total,F->maxEval,p,pv,&pval,
			F->repeatEval,F->PM);
	    if(repeats < 0) repeats = -repeats;
	    if(repeats >= min_rpt){

	      // 1. Find optimum number of gapped repeats.
	      best_operation=0; lastscore=0;
	      for(r = repeats; r <= max_rpts; r++){
        	operation=gapped_aln_seq_smatrixSW(a,b,LenSeq(E),SeqPtr(E),
		   r*N,SMatricesPrtnModel(PM),GapScoresPrtnModel(PM),
		   &start,&score);
		if(r == repeats){  // guarranteed to be best.
		  this_score = score/r;  
		} else { this_score = score-lastscore; }
		lastscore = score;
		pval=EvaluePrtnModel(this_score, LenSeq(E), PM);
#if 0
		fprintf(stderr,"rpts = %d; score = %d; this_score = %d; Eval = %g\n",
			r,score,this_score,pval);
#endif
		if(pval <= F->repeatEval){
			if(best_operation) free(best_operation);
			best_operation=operation;  best_rpts=r;
		} else { free(operation);  break; }
	      }
	      if(best_operation){
#if 1
		PutSeqInfo(ofp,E);
                // put_seqaln_smatrixSW(ofp,best_operation,LenSeq(E),SeqPtr(E),0,
                put_seqaln_smatrixSW(ofp,best_operation,LenSeq(E),SeqPtr(E),OffSetSeq(E),
		   start,best_rpts*N,SMatricesPrtnModel(PM));
		fprintf(ofp,"\n  (%d repeats)\n\n",best_rpts);
                // fprintf(ofp,"operations = %s\n",best_operation);
#endif
	        gap_seq_goscan(best_operation,1,F,gstrt,gend);
		NEW(ListE[s], best_rpts+2, e_type); hits++;
	        for(Int4 j=1; j <= best_rpts; j++){
		   st = gstrt[j]-left; st = MAXIMUM(Int4, st, 1);
		   e = gend[j]+right; e = MINIMUM(Int4, e, LenSeq(E));
                   ListE[s][j] = MkSubSeq(st,e,E);
		   // PutSubSeq(stdout,gstrt[j]-left,gend[j]+right,E,PrtnModelA(PM));
		   // NOTE PutSubSeq( ) adjusts boundary extensions.
	        } free(best_operation);
	      } NilSeq(E);	// destroy for now... May want to return...
	    } else { NilSeq(E); } 
	}
	free(pv); free(p); free(gstrt); free(gend); 
	e_type **rtnE;
	NEWP(rtnE, hits + 2, e_type);
	for(r=1,s=1; s<=number; s++){
	   if(ListE[s] != NULL){ rtnE[r] = ListE[s]; r++; }
	} free(ListE);
	return rtnE;
}

snh_typ	GOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	gsn_typ F)
{ return GOScanScan(fptr,number,total,nsize,1,0,F); }

snh_typ	GOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	Int4 min_rpt, gsn_typ F)
{ return GOScanScan(fptr,number,total,nsize, min_rpt, 0,F); }

snh_typ	GOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	Int4 min_rpt, Int4 MaxBadBlks, gsn_typ F)
/*** impose colinearity constraint and return heap with results. ***/
{
	e_type	E;
	a_type	A;
	Int4	s,hpsz,totN,*p,repeats,N,*lengths;
	char	mode;
	double	pval;
	float	*pv;
	snh_typ	sH;
	sni_typ	I,List=NULL;

	totN=NmaxPrtnModel(F->PM); N=NumModelsPrtnModel(F->PM);
	NEW(p,totN+3,Int4); NEW(pv,totN+3,float); 
	mode = ModePrtnModel(F->PM); A = PrtnModelA(F->PM);
        for(s=1; s<=number;s++) {
	    E = ReadSeq(fptr,s,nsize[s],A);
	    if(F->shuffle) ShuffleSeq(E);
	    repeats=ComparePrtnModelRpts(E,total,F->maxEval,p,pv,&pval,
			F->repeatEval,F->PM);
	    if(repeats >= min_rpt){
		I=MakeScanInfo(E,pval,repeats*N,repeats,pv,p,s,mode);
	        List = AppendScanInfo(I, List);
	    } else if(repeats <= -min_rpt){
		repeats = -repeats;
		I=MakeScanInfo(E,pval,repeats*N,repeats,pv,p,s,toupper(mode));
	        List = AppendScanInfo(I, List);
	    } else {
		if(F->dfp != NULL) PutSeq(F->dfp,E,A); 
		NilSeq(E); 
	    } 
	}
	hpsz = LengScanInfo(List);
	lengths = PrtnModelLengths(F->PM);
	sH = MakeScanHeap(hpsz+1,N,lengths,F->singleEval,MaxBadBlks);
	if(List != NULL){
	  while((I = RmLastScanInfo(List)) != List){
		if(InsertScanHeap(I, sH)==0) NilScanInfo(I);
	  }
	  if(InsertScanHeap(I, sH)==0) NilScanInfo(I);
	}
	free(pv); free(p); 
	return sH;
}

void	GOScanSWScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	gsn_typ F)
{ GOScanSWScan(fptr,number,total,nsize, 10,5,F); }

void	GOScanSWScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	Int4 open, Int4 extend,gsn_typ F)
{ GOScanSWScan(fptr,number,total,nsize, open, extend,F,' '); }

void	GOScanSWScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	Int4 open, Int4 extend,gsn_typ F,char mode)
// search a database for full Smith-Waterman alignment score against models
{
	e_type	E;
	a_type	A;
	Int4	s,totN,*p,N;
	double	pval;
	float	*pv;

	totN = NmaxPrtnModel(F->PM); 
	NEW(p,totN+3,Int4); NEW(pv,totN+3,float); 
	N = NumModelsPrtnModel(F->PM);
	A = PrtnModelA(F->PM);
	if(mode == 'C'){
	   Int4 End=TotLenPrtnModel(F->PM);
	   fprintf(stdout,
		"[0_(1)=goscan(%d){go=10000,gx=2000,pn=1000.0,lf=0,rf=0}:",number);
	   fprintf(stdout,"\n(%d)",End);
	   for(s=1; s <= End; s++){
		fprintf(stdout,"*");
	   } fprintf(stdout,"\n\n");
	}
        for(s=1; s<=number;s++) {
	  E = ReadSeq(fptr,s,nsize[s],A);
	  if(F->shuffle) ShuffleSeq(E);
	  if(mode == 'C'){ fprintf(stdout,"$%d=",s); }
	  ComparePrtnModelSW(E,total,F->maxEval,&pval,open,extend,F->PM,mode);
	  if(pval <= F->maxEval && F->dfp != NULL) PutSeq(F->dfp,E,A); 
	  NilSeq(E); 
	}
	if(mode == 'C') fprintf(stdout,"_0].\n\n");
}


snh_typ	*MGOScanScan(FILE *fptr,Int4 number,UInt8 total,unsigned short *nsize, 
	gsn_typ F)
// impose colinearity constraint and return heap with results for multiple 
// gapped models.
{
	e_type	E;
	a_type	A;
	Int4	i,s,hpsz,n,n_max,repeats,N,*lengths;
	Int4	*p,*p_max,*p_tmp,rpts_max,i_max;
	char	mode;
	double	pval,pval_max;
	float	*pv,*pv_max,*pv_tmp;
	snh_typ	*sH;
	sni_typ	I,*List=NULL;

	NEW(sH, F->num_models+3, snh_typ);
	NEW(List, F->num_models+3, sni_typ);
        for(n_max=0,i=1; i<=F->num_models;i++) {
		n=NmaxPrtnModel(F->pm[i]); 
		if(n > n_max) n_max=n;
	}
	NEW(p,n_max+3,Int4); NEW(pv,n_max+3,float); 
	NEW(p_max,n_max+3,Int4); NEW(pv_max,n_max+3,float); 
	A = PrtnModelA(F->PM);
        for(s=1; s<=number;s++) {
	  E = ReadSeq(fptr,s,nsize[s],A);
	  if(F->shuffle) ShuffleSeq(E);
	  for(pval_max=DBL_MAX,rpts_max=0,i=1; i <= F->num_models; i++){
	    N=NumModelsPrtnModel(F->pm[i]);
	    mode = ModePrtnModel(F->pm[i]);
	    repeats=ComparePrtnModelRpts(E,total,F->maxEval,p,pv,&pval,F->repeatEval,F->pm[i]);
	    if(repeats > 0 && pval < pval_max){
		pval_max=pval;
		pv_tmp = pv_max; pv_max = pv; pv = pv_tmp;
		p_tmp = p_max; p_max = p; p = p_tmp;
		i_max = i; rpts_max = repeats;
	    }
	  }
	  if(rpts_max > 0){
	  	N=NumModelsPrtnModel(F->pm[i_max]);
		I=MakeScanInfo(E,pval_max,rpts_max*N,rpts_max,pv_max,p_max,s,mode);
	        List[i_max] = AppendScanInfo(I, List[i_max]);
	  } else if(rpts_max < 0){
		rpts_max= -rpts_max;
	  	N=NumModelsPrtnModel(F->pm[i_max]);
		I=MakeScanInfo(E,pval_max,rpts_max*N,rpts_max,pv_max,p_max,s,toupper(mode));
	        List[i_max] = AppendScanInfo(I, List[i_max]);
	  } else { if(F->dfp != NULL) PutSeq(F->dfp,E,A); NilSeq(E); } 

	}
	for(i =1; i <= F->num_models; i++){
	  N=NumModelsPrtnModel(F->pm[i]);
	  hpsz = LengScanInfo(List[i]);
	  lengths = PrtnModelLengths(F->pm[i]);
	  sH[i] = MakeScanHeap(hpsz+1,N,lengths,F->singleEval);
	  if(List[i] != NULL){
		while((I = RmLastScanInfo(List[i])) != List[i]){
		  if(InsertScanHeap(I, sH[i])==0) NilScanInfo(I);
	  	}
	  	if(InsertScanHeap(I, sH[i])==0) NilScanInfo(I);
	  }
	}
	free(List);
	free(pv); free(p); free(pv_max); free(p_max); 
	return sH;
}

// NEW GAPPED ROUTINE FOR GENERATING CMA_TYP FILES

Int4    gap_seq_goscan(char *operation,Int4 start0,gsn_typ F,Int4 *left,
	Int4 *right, char **SubOper)
/****************************************************************************
  create 'fake' subsequence with gaps and flanking
  trace = "EDdmmmmmmmmmmmiiiMmmmmmmmmmmmmmmmmmmmmiiMmmmmmmmmmmE"
 ****************************************************************************/
{
        Int4     o,i,j,t,n,len,offsetX=0,start=0,end=0,begin,endss,b;
        char     state;

        assert(operation[0]=='E');
        assert(operation[strlen(operation)-1]=='E');

        Int4 numM=NumModelsPrtnModel(PrtnModelGOScan(F));
        // 2. Add aligned-region sequence.
        for(j=start0,end=t=0,n=o=1,state='E'; operation[o] != 'E'; o++){
            switch(operation[o]){
               case 'M':
                  if(t==0){ start=j; begin=o; }
                  if(t == numM){ // end of profile...
                        if(state != 'i'){ endss=o-1; end =j-1; }
                        left[n]=start; right[n]=end; 

		        len = endss-begin+1;
			NEW(SubOper[n],len+4,char); SubOper[n][0]='E';
			for(b=1,i=begin; b <= len; b++,i++)
				SubOper[n][b]=operation[i]; 
			SubOper[n][b++]='E'; SubOper[n][b]=0;
	
                        n++; start=j; begin=o; t=1;
                  } else t++;
               case 'm': j++; break;
               case 'D': // starting deletion in sequence relative to profile.
                  if(t==0){ start=j; begin = o; }
                  if(t == numM){ // end of sequence...
                        if(state != 'i') { endss=o-1; end =j-1; }
                        left[n] = start; right[n] = end; 

			len = endss-begin+1;
			NEW(SubOper[n], len+4,char); SubOper[n][0]='E';
			for(b=1,i=begin; b <= len; b++,i++)
				SubOper[n][b]=operation[i]; 
			SubOper[n][b++]='E'; SubOper[n][b]=0;
	
                        n++; start =j; begin=o; t=1;
                  } else t++;
               case 'd': // internal deletion in sequence relative to profile.
                  break;
               case 'i': // insert is between profile blocks;
                  assert(state != 'i');
                  if(t == numM) { endss=o-1; end=j-1; }
                  while(operation[o]=='i'){ j++; o++; } o--;
                  break;
               case 'I': // Insert within a profile block; delete from seq.
                  assert(state != 'I');
                  while(operation[o]=='I'){ j++; o++; } o--;
                  break;
               default:
                        fprintf(stderr,"operations = %s\n",operation);
                        print_error("gap_seq_goscan( ): input error"); break;
            }  state=operation[o];
        }
        if(state != 'i') { endss=o-1; end =j-1; }
        left[n] = start; right[n] = end;

	len = endss-begin+1;
	NEW(SubOper[n], len+4,char); SubOper[n][0]='E';
	for(b=1,i=begin; b <= len; b++,i++) SubOper[n][b]=operation[i]; 
	SubOper[n][b++]='E'; SubOper[n][b]=0;
	
// fprintf(stderr,"DEBUG1: s=%d; e=%d; %d ; %d\n",s,e,start,end);
        return n;
}

Int4	GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 gapo, Int4 gapx,
	Int4 left, Int4 right, Int4 min_rpt, Int4 number,UInt8 total, 
	unsigned short *nsize, char ***Operation, e_type **RtnE, e_type **FullE, 
	unsigned short **FullR, Int4 **Start, char Mode, gsn_typ F)
{
	return GapSeqGOScanScan(fptr,a,b,gapo,gapx,left,right,min_rpt,number,
		total,nsize,Operation,RtnE,FullE,FullR,Start,Mode,F,0,UINT4_MAX);
}

Int4	GapSeqGOScanScan(FILE *fptr,Int4 a, Int4 b, Int4 gapo, Int4 gapx,
	Int4 left, Int4 right, Int4 min_rpt, Int4 number,UInt8 total, 
	unsigned short *nsize, char ***Operation, e_type **RtnE, e_type **FullE, 
	unsigned short **FullR, Int4 **Start, char Mode, gsn_typ F, 
	UInt4 minlen, UInt4 maxlen)
// find optimum gapped sequence repeats 
{
	e_type	E,*fullE;
	a_type	A;
	Int4	r,s,st,e,totN,*p,repeats,N;
	char	mode;
	double	pval;
	ptm_typ PM=PrtnModelGOScan(F);
	char	*operation,*best_operation,***gstr;
	Int4	start,score,lastscore,this_score;
	float	*pv;
	Int4	**gstrt,*gend,max_rpts,best_rpts;
	unsigned short	*fullR,*fullR0;
// Int4 key_hit=0;

	assert(!F->shuffle); assert(min_rpt > 0);
	totN=NmaxPrtnModel(F->PM); N=NumModelsPrtnModel(F->PM);
	NEW(p,totN+3,Int4); NEW(pv,totN+3,float); 
	mode = ModePrtnModel(F->PM); A=PrtnModelA(F->PM);
	max_rpts=MaxRptsPrtnModel(F->PM);
	NEW(gend,max_rpts+2,Int4);
	Int4 *lenM=PrtnModelLengths(F->PM);

	e_type	**ListE;
	Int4	fullhits,hits;
	NEWP(ListE, number +2, e_type);
	NEW(ListE[0], number +2, e_type); 
	NEWPP(gstr, number +2, char);
	NEWP(gstrt, number +2, Int4);
	NEW(fullR0, number +2, unsigned short);
//************************** HEAP FOR psgaln scan **********************
// Use best single repeat score for key right now...
	pah_typ *pah=0;
	cls_typ	*cls=0;
	if(Mode == 'H'){	 // Test Heap for scan
	   Int4 maxlength;
           for(maxlength=0, s=1; s<=number;s++) { 
		if(nsize[s] > maxlength) maxlength = nsize[s];
	   }
	   cls = new cls_typ(A,maxlength);
	   pah = new pah_typ(1000);
	}
//************************** HEAP FOR psgaln scan **********************
        for(fullhits=hits=0, s=1; s<=number;s++) {
	    E = ReadSeq(fptr,s,nsize[s],A); assert(E);
	    if(LenSeq(E) < minlen || LenSeq(E) > maxlen) repeats=0;
	    else {
	        repeats=ComparePrtnModelRpts(E,total,F->maxEval,p,pv,&pval,
			F->repeatEval,F->PM);
	       if(repeats < 0) repeats = -repeats;
	    }
	    // if(repeats >= min_rpt){  // TOO EARLY TO TELL!!!
	    if(repeats > 0){ 
	      // 1. Find optimum number of gapped repeats.
	      best_operation=0; lastscore=0;
	      for(r = repeats; r <= max_rpts; r++){
        	operation=gapped_aln_seq_smatrixSW(gapo,gapx,LenSeq(E),SeqPtr(E),
		   r*N,SMatricesPrtnModel(PM),GapScoresPrtnModel(PM),
		   &start,&score);
		if(r == repeats){  // guarranteed to be best.
		  this_score = score/r;  
		} else { this_score = score-lastscore; }
		lastscore = score;
		pval=EvaluePrtnModel(this_score, LenSeq(E), PM);
		if(pval <= F->repeatEval){
			if(best_operation) free(best_operation);
			best_operation=operation;  best_rpts=r;
		} else { free(operation);  break; }
	      }
	      if(best_operation){
		// fprintf(stderr,"operations(%d,%d) = %s\n",gapo,gapx,best_operation);
		free(best_operation); best_operation=0;
		if(F->pssm == 0){
		  if(best_rpts >= min_rpt){
        	    best_operation=gapped_aln_seq_smatrixSW(a,b,LenSeq(E),
		       SeqPtr(E),best_rpts*N,SMatricesPrtnModel(PM),
		       GapScoresPrtnModel(PM), &start,&score);
		    // fprintf(stderr,"start = %d\n",start);
		    start=1;  // this method keeps 'i' insertion residues.
		  }
		} else if(F->hmm){	
		  if(Mode == 'S'){	
		    Int4 tmp_rpts,oper_len;
		    best_operation=F->hmm->FindBestRpts(E,best_rpts,&tmp_rpts,
					&score,&start,&oper_len);
		    best_rpts=tmp_rpts;
		    // best_operation=F->hmm->Align(stderr,E,best_rpts,&score,&start,&oper_len);
		    fprintf(stderr,"total score = %d\n",score);
	    	    if(best_rpts < min_rpt){ free(best_operation); best_operation=0; }
                  } else if(best_rpts >= min_rpt){
                    // assert(N == 1);
                    Int4 oper_len;
                    best_operation=F->hmm->Align(E,best_rpts,&score,&start,&oper_len);
		    F->hmm->PutAlign(stdout,E,best_operation, oper_len,start,best_rpts);
                    fprintf(stderr,"score = %d; repeats = %d\n",score,best_rpts);
                  }
		} else if(Mode == 'S'){	
		  // Test simulation...
		  // use sampling procedure...then add all significant repeats 
		  // assert(N == 1);	// Only for 1 block models for now...
		  Int4 tmp_rpts,oper_len;
		  best_operation=F->pssm->FindBestRpts(E,best_rpts,&tmp_rpts,&score,&start,&oper_len);
		  best_rpts=tmp_rpts;
		  // best_operation=F->pssm->Align(stderr,E,best_rpts,&score,&start,&oper_len);
		  fprintf(stderr,"total score = %d\n",score);
	    	  if(best_rpts < min_rpt){ free(best_operation); best_operation=0; }
//************************** HEAP FOR psgaln scan **********************
		} else if(Mode == 'H'){	 // Test Heap for scan
	    	  if(best_rpts >= min_rpt){
		    Int4 oper_len;
		    best_operation=F->pssm->Align(E,best_rpts,&score,&start,&oper_len);
		    // best_operation=F->pssm->Align(stderr,E,best_rpts,&score,&start,&oper_len);
		    // fprintf(stderr,"score = %d\n",score);
		    //****************************************************
		    // Get score for heap (probably want to use seg and/or coils here)
		    //****************************************************
		    Int4 tmp_score,tmp_start;
		    unsigned char tmp_rpts=1;
		    if(best_rpts > 4 && FALSE) tmp_rpts=6; else tmp_rpts=4;
		    e_type tmpE=CopySeq(E);
		    if(cls->MaskCoils(tmpE,14,0.85)){ // mask coiled coils
			cls->Put(stderr,E,0.50);
		    }
		    BooLean masked=ProcessSeqPSeg(17,2.2,2.5,100,tmpE,A);
		    char *tmp_operation=F->pssm->Align(tmpE,tmp_rpts,
				&tmp_score,&tmp_start,&oper_len);
		    if(!pah->Insert(tmpE, tmp_operation,tmp_score,tmp_rpts,tmp_start)){
			free(tmp_operation); NilSeq(tmpE);
		    }
		  }
//************************** HEAP FOR psgaln scan **********************
		} else {	// use Alex's PSSM.
	    	  if(best_rpts >= min_rpt){
		    // assert(N == 1);
		    Int4 oper_len;
		    best_operation=F->pssm->Align(E,best_rpts,&score,&start,&oper_len);
		    // best_operation=F->pssm->Align(stderr,E,best_rpts,&score,&start,&oper_len);
		    // fprintf(stderr,"score = %d\n",score);
		  }
		}
	        if(best_operation){
		 NEWP(gstr[s],best_rpts+2,char);
		 NEW(gstrt[s],best_rpts+2,Int4);
		 // fprintf(stderr,"operations(%d,%d) = %s\n",a,b,best_operation);
	         gap_seq_goscan(best_operation,start,F,gstrt[s],gend,gstr[s]);
		 NEW(ListE[s], best_rpts+2, e_type); 
		 hits+=best_rpts; fullhits++; fullR0[s]=best_rpts;
	         for(Int4 j=1; j <= best_rpts; j++){
		   st = gstrt[s][j]-left; st = MAXIMUM(Int4, st, 1);
		   e = gend[j]+right; e = MINIMUM(Int4, e, LenSeq(E));
                   ListE[s][j] = MkSubSeq(st,e,E);
		   // reset start to correspond to subseq.
		   gstrt[s][j] = MINIMUM(Int4,left+1,gstrt[s][j]);
		   // PutSubSeq(stdout,gstrt[s][j]-left,gend[j]+right,E,PrtnModelA(PM));
		   // NOTE PutSubSeq( ) adjusts boundary extensions.
	         }
		 free(best_operation);
		} else { NilSeq(E); E=0; }
	      } else { NilSeq(E); E=0; }  // end if(best_operation) [above] 
	      if(E) ListE[0][s] = E;
	    } else { NilSeq(E); E=0; }  // end if(repeats > 0)...
	} free(pv); free(p); free(gend); 
//************************** HEAP FOR psgaln scan **********************
	if(Mode == 'H'){	 // Test Heap for scan
		if(pah){
		    Int4 tmp_score,tmp_start;
		    unsigned char tmprpts;
		    char *tmp_operation;
		    e_type tmpE;
		    pah->Put(stdout);
		    while(pah->DelMax(&tmpE,&tmp_operation, &tmprpts,
						&tmp_start,&tmp_score)){
		     fprintf(stdout,"======================\n");
		     PutSeqInfo(stdout,tmpE);
		     F->pssm->PutAlign(stdout,tmpE,tmp_operation,
				strlen(tmp_operation),tmp_start,tmprpts);
		     fprintf(stdout,"Score = %d\n",tmp_score);
		     NilSeq(tmpE); free(tmp_operation);
		   }
		    delete pah;
		    delete cls;
		}
	}
//************************** HEAP FOR psgaln scan **********************
	e_type	*rtnE;
	char	**rtnStr;
	Int4	*rtnStrt;
	NEW(rtnE, hits + 2, e_type);
	NEW(fullR, fullhits + 2, unsigned short);
	NEW(fullE, fullhits + 2, e_type);
	NEW(rtnStrt, hits + 2, Int4);
	NEWP(rtnStr, hits + 2, char);
	Int4 f;
	for(f=r=0,s=1; s<=number; s++){
	    if(ListE[s] != NULL){ 
		f++; fullE[f] = ListE[0][s]; fullR[f] = fullR0[s];
		for(Int4 t=1; ListE[s][t]; t++){
		   r++;
		   rtnE[r] = ListE[s][t]; 
		   rtnStrt[r] = gstrt[s][t];
		   rtnStr[r] = gstr[s][t];
		}
		free(ListE[s]); free(gstrt[s]); free(gstr[s]);
	    }
	} free(ListE[0]); free(ListE); free(gstrt); free(gstr);
	free(fullR0);
	*RtnE = rtnE; *Operation = rtnStr; *Start = rtnStrt;
	*FullR = fullR; *FullE = fullE;
	return r;
}

#define	USAGE_GSN2CMA	"USAGE: scan2cma database msafile [options]\n\
   msafile = file with aligned segments\n\
   options:\n\
     -B<char>   - define residue frequencies to use for background model\n\
                  m = use multiple alignment background counts.\n\
                  d = use database background counts.\n\
                  s = use standard background counts.\n\
                  (default = use standard background counts).\n\
     -C<float>  - minimum field map Cutoff for individual blocks (0.0)\n\
     -D         - Don't mask low complexity regions (will mask by default)\n\
     -d<float>  - maximum E-value for repeat domain detection\n\
     -e<float>  - maximum E-value for printing alignments\n\
     -E<float>  - maximum E-value for single motif block\n\
     -F<file>   - add full_counts to cma file & domain information from <file>\n\
     -f         - add full_counts to cma file\n\
     -G=<int>,<int>,<int> - gapo, gapx, pernats (default: 25,4,5)\n\
     -g<int>,<int> - gapped repeat opening and extension penalties (default: 18,2)\n\
     -h         - use hmm type search\n\
     -I<int>:<int> - left & right flank lengths for domain sampling\n\
     -M<int>    - maximum number of input sequences (default: 1,000,000)\n\
     -m<char>   - method = <char> (lowercase = use informative columns only)\n\
                  M or m = modified Gribskov method.\n\
                  F or f = modified Gribskov method with off columns=blosum45.\n\
                  D or d = Dirichlet Mixtures priors.\n\
                  R or r = product multinomial.\n\
                  H or h = Henikoff's method\n\
                  B or b = Henikoff's method with motif residue background\n\
     -N<int>    - maximum number of repeats to look for\n\
     -n         - mask potential nonglobular regions using 'seg x 45 3.4 3.75'\n\
     -P         - use complex PSSM protein profile with default settings or...\n\
     -P<int>..<int>,<int>..<int>:<int>,<int>..<int>\n\
     -P<int>..<int>,<int>..<int>/<int>..<int>:<int>,<int>..<int>\n\
                - with specified insertion & deletion opening & extension penalties\n\
     -p<float>  - pseudo counts for product multinomial model\n\
     -r<int>    - minimum number of repeats to output a hit (1-1000:1)\n\
     -S<int>:<int> - size range of sequences to be searched (default: all seqs)\n\
     -T<char>   - test mode\n\
                  'S' = Sample to determine repeat cutoff\n\
     -t<real>   - set target p-value for repeats with -P option (default: 0.05)\n\
     -u<char>   - scan method with or without gap function (default 'O')\n\
                  'g' = global with gap function\n\
                  'G' = global without gap function\n\
                  'd' = local with gap function\n\
                  'D' = local without gap function\n\
                  'o' = global overlap with gaps\n\
                  'O' = global overlap without gaps\n\
                  'c' = local core with gaps\n\
                  'C' = local core without gaps\n\
     -s<int>    - random seed\n\
     -v         - verbose mode\n\
     -w         - DON'T apply weights to sequences in scan file\n\n"

cma_typ	GOScanToCMSA(Int4 argc,char *argv[], a_type A)
{ 
	Int4	m,n,i,arg,number,maxrpts=1,right=5,left=5;
	Int4	*counts,time1,MAX_IN_SEQS=10000000,*mtfcnts=NULL;
	UInt8	total;
	Int4	gapo=25,gapx=4,pernats=5,a=18,b=2;
	UInt4 min_rpt=1;
	unsigned short	*nsize;
	char	aafreq='s',mode='O',method='H';
	float	minmap=0.0,domEval=0.05,expect=0.01,pseudo=0.5,singleEval=1.0;
	float	target_pvalue=0.05;
	double	*freq;
        FILE    *fptr;
	gsn_typ F;
	BooLean	segmask=TRUE, mask_nonglobular=FALSE,weights=TRUE,AddFull=FALSE;
	BooLean	verbose=FALSE;
	char	Test=' ';
	BooLean	use_hmm=FALSE;
	char	*dom_file=0;
	static char default_arg[ ] = "-P350..1700,20..50:400..500,40..700";
	char	*pssm_arg=0;
	UInt4 seed=18364592;
	Int4	minlen=0,maxlen=INT4_MAX;

// mode='G'; // use gapped not overlap...

	time1=time(NULL);
	if(argc < 3) print_error(USAGE_GSN2CMA);
	// A=MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	for(arg = 3; arg < argc; arg++){
	   if(argv[arg][0] != '-') print_error(USAGE_GSN2CMA);
	   switch(argv[arg][1]) {
	     case 'B': if(!isalpha(aafreq=argv[arg][2])) 
			print_error(USAGE_GSN2CMA); break;
	     case 'C': minmap=RealOption(argv[arg],'C',-1000,+5000,USAGE_GSN2CMA);
		break;
	     case 'D': segmask= FALSE; break;
	     case 'd': domEval=RealOption(argv[arg],'d',0,500000,USAGE_GSN2CMA); 
		break;
	     case 'E': singleEval=RealOption(argv[arg],'E',0,10000,USAGE_GSN2CMA); 
		break;
	     case 'e': expect=RealOption(argv[arg],'e',0,500000,USAGE_GSN2CMA); 
		break;
	     case 'F': if(!isprint(argv[arg][2])) print_error(USAGE_GSN2CMA); 
		       else dom_file=AllocString(argv[arg]+2); AddFull=TRUE; break;
	     case 'f': AddFull=TRUE; break;
	     case 'G': if(sscanf(argv[arg],"-G=%d,%d,%d",&gapo,&gapx,&pernats) != 3) 
			print_error(USAGE_GSN2CMA); 
		     if(gapo < 1 || gapx < 1 || pernats < 1) print_error(USAGE_GSN2CMA);
		     break;
	     case 'g': if(sscanf(argv[arg],"-g%d,%d",&a,&b) != 2) 
			print_error(USAGE_GSN2CMA); 
		     if(a < 0 || b < 0) print_error(USAGE_GSN2CMA); break;
	     case 'h': use_hmm=TRUE; break;
	     case 'I': if(sscanf(argv[arg],"-I%d:%d",&left,&right) != 2)
			print_error(USAGE_GSN2CMA); break;
	     case 'M': MAX_IN_SEQS=IntOption(argv[arg],'M',
			1000,2000000000,USAGE_GSN2CMA); break;
	     case 'm': if(!isalpha(method=argv[arg][2])) 
			print_error(USAGE_GSN2CMA); break;
	     case 'N': maxrpts=IntOption(argv[arg],'N',1,1000,USAGE_GSN2CMA); break;
	     case 'n': mask_nonglobular=TRUE; break;
	     case 'P': pssm_arg=argv[arg]; break;
	     case 'p': pseudo=RealOption(argv[arg],'p',0,500,USAGE_GSN2CMA); break;
             case 'r': min_rpt=IntOption(argv[arg],'r',1,1000,USAGE_GSN2CMA); break;
	     case 'S': if(sscanf(argv[arg],"-S%d:%d",&minlen,&maxlen) != 2) 
			print_error(USAGE_GSN2CMA); 
		     if(minlen < 0 || maxlen < minlen) print_error(USAGE_GSN2CMA); break;
	     case 'T': Test=argv[arg][2]; if(Test==0) print_error(USAGE_GSN2CMA); break;
	     case 't': target_pvalue=RealOption(argv[arg],'t',0,1.0,USAGE_GSN2CMA); break;
	     case 's': seed=atoi(argv[arg]+2); break;
	     case 'u': if(!isalpha(mode=argv[arg][2])) 
			print_error(USAGE_GSN2CMA); break;
	     case 'w': weights = FALSE; break;
	     case 'v': verbose= TRUE; break;
	     case ' ': break; // ignore these 
	     default: print_error(USAGE_GSN2CMA);
	   }
	}

	number = GetFastaInfo(argv[1], MAX_IN_SEQS, &counts, &nsize, A);
	for(total=0, i=0; i<=nAlpha(A); i++) total += counts[i];
	for(arg = 0; arg < argc; arg++) {
		if(argv[arg][1] != ' ') fprintf(stderr,"%s ",argv[arg]);
	}
	if(seed == 18364592) {  // not provided by user
          seed = (UInt4) time(NULL);
          fprintf(stderr,"-s%d\n",seed);
   	} else fprintf(stderr,"\n");
	sRandom(seed);

	fptr = open_file(argv[1],"","r");
	NEW(freq,nAlpha(A)+2,double);
	if(aafreq == 'm') {
	  NEW(mtfcnts, nAlpha(A)+2, Int4);
	  sma_typ MA=ReadSMA(argv[2]); CountsSMA(mtfcnts, MA); NilSMA(MA);
	  for(n=0, i=0; i<=nAlpha(A); i++) n += mtfcnts[i];
    	  for(i=0; i<=nAlpha(A); i++) freq[i]= (double)mtfcnts[i]/(double)n;
	  free(mtfcnts);
	} else if(aafreq == 'd') {
    	  for(i=0; i<=nAlpha(A); i++) freq[i]= (double)counts[i]/(double)total;
	} else for(i=0; i<=nAlpha(A); i++) freq[i] = blosum62freq[i];
	for(m=0,i=1; i<= number; i++) m = MAXIMUM(Int4,nsize[i],m);
	F=MakeGOScan(argv[2],A,maxrpts,method,mode,expect,singleEval,m,weights,
		minmap,pseudo,freq);
	F->hmm=0;
	SetRptEvalGOScan(domEval,F);
	if(!segmask) NoMaskGOScan(F);
	if(mask_nonglobular){ MaskNonGlobularGOScan(F); }

        char    	**operation;
        e_type 		*ListE,*FullE;
	unsigned short	*FullR;
        Int4    	*start;
	cma_typ		cma=0;
        ptm_typ 	PM=PrtnModelGOScan(F);
        Int4    	nblks=NumModelsPrtnModel(PM);
	BooLean		*skip;

	if(use_hmm && pssm_arg==0) pssm_arg=default_arg;
	if(pssm_arg){ // default = "-P100..1100,20..120:500,40..400"
	  char	str[200];
	  Int4 len=strlen(argv[2]); len -= 4;  // delete .msa;
	  for(i=0; i < len; i++) str[i]=argv[2][i]; str[i]=0;
          strcat(str,".cma");
	  std::cerr << str; std::cerr << "\n";
          F->pssm_cma = ReadCMSA2(str,PrtnModelA(PM));
#if 1	// NEED to REMOVE POOR BLOCKS JUST AS FOR PM...
	  BooLean *remove=0;
	  for(i=nBlksCMSA(F->pssm_cma); i > 0; i--){
		double map = FieldRelMapCMSA(F->pssm_cma,i);
		if(map < minmap){
		   if(!remove) NEW(remove,nBlksCMSA(F->pssm_cma) + 2, BooLean);
		   remove[i]=TRUE;
		}
	  }
	  if(remove){
	   for(i=nBlksCMSA(F->pssm_cma); i > 0; i--){
	    if(remove[i]){
		cma = RmBlkCMSA(i,F->pssm_cma);
		NilCMSA(F->pssm_cma); F->pssm_cma=cma;
	    }
	   } free(remove); cma=0;
	  }
#endif
	  if(use_hmm){
		F->pssm=new psm_typ(MaxRptsPrtnModel(PM),0,F->pssm_cma);
	  	F->hmm= new HMM_typ(MaxRptsPrtnModel(PM),pssm_arg,F->pssm_cma,200,0);
	  	F->hmm->SetRptPval(target_pvalue);
                if(verbose) F->hmm->Put(stderr);
          } else {
		F->pssm = new psm_typ(MaxRptsPrtnModel(PM),pssm_arg,F->pssm_cma);
	  	if(verbose) F->pssm->Put(stderr);
	  	F->pssm->SetRptPval(target_pvalue);
		F->hmm=0;
	  }
	}
        Int4 NumHits=GapSeqGOScanScan(fptr,a,b,gapo,gapx,left,right,min_rpt,number,
		total,nsize,&operation,&ListE,&FullE,&FullR,&start,Test,F,minlen,maxlen);
	if(NumHits > 0){ 
          cma=MakeCMSA(ListE,NumHits,operation,start,nblks,PrtnModelLengths(PM),
			gapo,gapx,pernats,left,right,argv[1],PrtnModelA(PM),FullE,FullR);
			// 10000,2000,1000,left,right,argv[1],PrtnModelA(PM),FullE,FullR);
	  ss_type FullSeq = FullSeqCMSA(cma);
	  if(dom_file){
		FILE *fp = open_file(dom_file,"","r");
		dom_typ *dom;  dom=new dom_typ(fp); 
		fclose(fp);
		NEW(skip,number+2,BooLean);
		Int4 s,S;
		for(s=1; s <= number; s++) skip[s]=TRUE;
		for(s=1; s <= NSeqsSeqSet(FullSeq); s++){
		    e_type E = SeqSetE(s,FullSeq);
		    S = SeqI(E); skip[S] = FALSE;
		}
		fp=tmpfile(); dom->Put(fp,skip); rewind(fp);
		delete dom; dom = new dom_typ(fp); fclose(fp);
		AddDomainsCMSA(dom,cma);
	  }
	  if(!AddFull) RmFullCountsCMSA(cma);
	  if(verbose) { PutAlnCMSA(stdout,cma); }
	}  // WARNING: NOT DEALLOCATING ALL MEMORY WHEN NumHits == 0! Fix here later...
	free(operation); free(start);
	if(use_hmm) delete F->hmm;
	NilGOScan(F); fclose(fptr); 
	free(counts); free(nsize); free(freq);
	return cma;
}

