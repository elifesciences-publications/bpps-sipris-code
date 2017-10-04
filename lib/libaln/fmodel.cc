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

#include "fmodel.h"
#include "blosum62.h"

fm_type	MkFModel(BooLean *null, Int4 length, Int4 maxlen, double npseudo, 
	UInt4 *observedNS,UInt4 *counts, double *freq, a_type A)
/* create and return a fmodel with totsites = 0 and unfragmented columns. */
{
	fm_type	M;
	Int4	i,j,r,x;
	
	if(maxlen < length) {
		maxlen = length;
		// fmodel_error("maxlen must be >= length");
	}
	if(length*3 > maxlen) maxlen = length*3;
// fprintf(stderr,"MkFModel: %d max = %d\n",length,maxlen);
	NEW(M,1,fmodel_type);
	M->maxlen = maxlen;
	M->A = A; M->totsites = 0; M->length = length; 
	if(null != NULL){
	   for(M->ncols=length,j=1; j<=length; j++) if(null[j]) M->ncols--;
	   if(M->ncols < 3){
		fprintf(stderr,"length = %d; ncols = %d\n",length,M->ncols);
	        for(j=1; j<=length; j++){
			if(null[j]) fprintf(stderr,".");
			else fprintf(stderr,"*");
		}
		fprintf(stderr,"\n"); // assert(M->ncols >= 3);
		fmodel_error("model must have at least 3 columns");
	   }
	} else { M->ncols = length; }
	M->observedNS=observedNS;	// Borrowed from above...
	M->freq=freq;			// Borrowed from above...
	M->counts=counts;			// Borrowed from above...
	M->start = MAXIMUM(Int4,1,(M->maxlen - length)/2);
	M->end = M->start + length - 1;
	M->npseudo = MAXIMUM(double,npseudo,0.0001);
	NEWP(M->observed, M->maxlen+1, Int4);
	NEWP(M->likelihood, M->maxlen+1, double);
	NEW(M->sumlngamma, M->maxlen+1, double);
	NEW(M->tmp_val, M->maxlen+1, double);
	for(j = 0; j <= M->maxlen; j++){ M->observed[j] = NULL; }
	for(M->ncols = 0, i=1,j = M->start; j <= M->end; j++,i++){	
	  if(null == NULL || !null[i]){
	    NEW(M->observed[j],nAlpha(M->A) +2, Int4);
	    NEW(M->likelihood[j], nAlpha(M->A)+2, double);
            M->likelihood[j][0] = 0.05; /*** 'X' residues ***/
	    M->ncols++;
	  }
	}
	NEW(M->Ps,nAlpha(M->A)+2,double);
	NEW(M->targetNS,nAlpha(M->A)+2,double);
	NEWP(M->blsm62p,nAlpha(M->A) +2, double);
	for(r=0; r<= nAlpha(M->A); r++){
	     NEW(M->blsm62p[r],nAlpha(M->A) +2, double);
	     for(x=0; x<= nAlpha(M->A); x++) {
		M->blsm62p[r][x]=blosum62P[r][x];
	     }
	}
	// OLD InitFModel():
	Int4	b;
	M->totsites = 0; 
	for(b=1; b<= nAlpha(M->A); b++) { M->Ps[b]= M->npseudo * M->freq[b]; }
	for( j = M->start; j <= M->end; j++){	
	  if(M->observed[j]!=NULL){
	   for(b=0; b <= nAlpha(M->A); b++) M->observed[j][b] = 0;
	  }
	}
	M->recalc = M->update = TRUE;
	return M;
}

BooLean	EnlargeFModel(Int4 maxlen, fm_type M)
// if maxlen > M->maxlen then reallocate arrays.
{
        Int4    oldmax,i,**observed,**lpp;
	double	**likelihood,*sumlngamma,**dpp,*dp;

#if 0
std::cerr << "enlarging fmodel...\n";
PutFModel(stderr, M); 
#endif

        if(maxlen <= M->maxlen) return FALSE;
	oldmax=M->maxlen;
        M->maxlen = maxlen;
        NEWP(observed, M->maxlen+1, Int4);  
	lpp=M->observed; M->observed=observed; observed=lpp;
        NEWP(likelihood, M->maxlen+1, double);
	dpp=M->likelihood; M->likelihood=likelihood; likelihood=dpp;
        NEW(sumlngamma, M->maxlen+1, double);
	dp=M->sumlngamma; M->sumlngamma=sumlngamma; sumlngamma=dp;
	for(i=0; i<=oldmax; i++){
		M->observed[i] = observed[i];
		M->likelihood[i] = likelihood[i];
		M->sumlngamma[i] = sumlngamma[i];
	}
	free(likelihood); free(observed); free(sumlngamma); 
	free(M->tmp_val); NEW(M->tmp_val, M->maxlen+1, double);
#if 0
std::cerr << "done enlarging fmodel\n";
PutFModel(stderr, M); 
#endif
	return TRUE;
}

fm_type	CopyFModel(fm_type M)
/* return a fmodel with totsites = 0 and unfragmented columns. */
{
	fm_type	M2;
	BooLean	*null;
	Int4	i,j;
	
	NEW(null,M->length+2,BooLean);
	for(i=1,j = M->start; j <= M->end; j++,i++){	
		if(M->observed[j] == NULL) null[i] = TRUE;
		else null[i] = FALSE;
	}
	M2 = MkFModel(null, M->length, M->maxlen, M->npseudo, M->observedNS,
		M->counts, M->freq, M->A);
	free(null);
	return M2;
}

fm_type NilFModel(fm_type M)
/* Destroy model M */
{
   Int4 j;

   if(M!=NULL){ 
	free(M->Ps); 
	free(M->targetNS); 
	free(M->sumlngamma);
	for(j=M->start; j<=M->end; j++){
	  if(M->observed[j]!=NULL){
	    free(M->observed[j]); free(M->likelihood[j]);
	  }
	}
	free(M->observed); 
	free(M->likelihood); free(M->tmp_val);
	for(j=0; j<= nAlpha(M->A); j++) free(M->blsm62p[j]);
	free(M->blsm62p);
	free(M);
   }
   return (fm_type) NULL;
}


void	RmFModel(unsigned char *seq, Int4 site, fm_type M)
/* remove the segment at site in seq from model */
{
        Int4 j;
        for(j=M->start; j<=M->end; j++,site++){
	    if(M->observed[j] != NULL){
		M->observed[j][seq[site]]--;
#if 0
if(M->observed[j][seq[site]] < 0){
std::cerr << "M->observed[j][seq[site]] = "; std::cerr << M->observed[j][seq[site]]; std::cerr << std::endl;
std::cerr << "site = "; std::cerr << site; std::cerr << std::endl;
std::cerr << "j = "; std::cerr << j-M->start+1; std::cerr << std::endl;
PutFModel(stderr, M); return FALSE;
	}
#endif
	    }
        }
	M->totsites--; M->recalc = M->update = TRUE;
}

void	Add2FModel(unsigned char *seq, Int4 site, fm_type M)
/* Add the segment at site in seq to model */
{
	Int4 j;
	for(j=M->start; j<=M->end; j++,site++){
	    if(M->observed[j] != NULL) {
	    	M->observed[j][seq[site]]++;
	    }
	}
	M->totsites++; M->recalc = M->update = TRUE;
}

BooLean	NullSiteFModel(Int4 s,fm_type M)
/* if site s == null in model M return TRUE; else return FALSE. */
{
	Int4	site;
	if(s < 1 || s > M->length) return TRUE;
	site = s + M->start - 1;
	if(M->observed[site] == NULL) return TRUE;
	else return FALSE;
}

Int4	NullSitesFModel(BooLean *null, fm_type M)
/* modifies null array so that null sites in model are TRUE and
   nonnull sites are FALSE; returns the length of model */
{
	Int4	i,j;

	for(i=1,j=M->start; j <= M->end; j++,i++){
		if(M->observed[j]==NULL) null[i] = TRUE;
		else null[i] = FALSE;
	}
	return M->length;
}

/************************ Column Move Routines ****************************/

Int4	WorstOnFModel(fm_type M)
/** Return the worst 'on' column in the model. */
{
	Int4	i,k,worst;
	double	r,min;

	min = RatioFModel(M->observed[M->start],1,M);
	for(worst=1,k=2,i=M->start+1; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		r = RatioFModel(M->observed[M->start],k,M);
		if(r < min){
			min=r; worst=k;
		}
	     }
	}
	return worst;
}

Int4	ChoiceLemonFModel(fm_type M)
/** Sample a column in model proportional to how BAD it is. */
{
	Int4	i,k;
	double	r,rand,total;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		r = RatioFModel(M->observed[M->start],k,M);
		M->tmp_val[k] = r; total += r;
	     }
	}
	rand = ((double) Random()/(double) RANDOM_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		rand -= M->tmp_val[k];
		if(rand <= 0.0) return k;
	     }
	}
	fprintf(stderr,"start=%d; end=%d; k=%d; i=%d; rand=%g; total=%g\n",
        	M->start,M->end, k,i,rand,total);
	for(k=1,i=M->start; i<=M->end; i++,k++){
		fprintf(stderr,"RatioFModel(%d) = %g\n",k,M->tmp_val[k]);
	}
	return INT4_MIN;	// temporary fix for overflow.
	fmodel_error("ChoiceLemonFModel( )... this should not happen.");
}

Int4	ChoiceOrangeFModel(fm_type M)
/** Sample a column in model at random. */
{
	Int4	i,k;
	double	rand,total;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL) total += 1.0;
	}
	rand = ((double) Random()/(double) RANDOM_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		rand -= 1.0;
		if(rand <= 0.0) { return k; }
	     }
	}
	fmodel_error("ChoiceOrangeFModel( )... this should not happen.");
}

Int4	LemonFModel(fm_type M)
/** Sample a column in model proportional to how BAD it is. */
// if overflow occurs return INT4_MIN.
{
	double	r,rand,total;
	Int4	i,k;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		r = RatioFModel(M->observed[M->start],k,M);
		// Fix overflow problem
		// M->tmp_val[k] = r; total += r;  // OLD..
		if(r > 1000000000.0){   // ratio > 1 billion
			M->tmp_val[k] = 1000000000.0; total +=1000000000.0;
		} else {
			M->tmp_val[k] = r; total += r;
		}
	     }
	}
	rand = ((double) Random()/(double) RANDOM_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL){
		rand -= M->tmp_val[k];
		if(rand <= 0.0) return k;
	     }
	}
	// overflow error below...
	fprintf(stderr,"k=%d; rand=%g; total=%g\n",k,rand,total);
	for(Int4 b=1;b<= nAlpha(M->A); b++) {
		fprintf(stderr,"observed[%c] = %d\n",
			AlphaChar(b,M->A),M->observed[M->start][b]);
	}
	for(k=1,i=M->start; i<=M->end; i++,k++){
		fprintf(stderr,"RatioFModel(%d) = %g\n",k,M->tmp_val[k]);
		if(M->tmp_val[k] >= DBL_MAX){
		  for(Int4 b=1;b<= nAlpha(M->A); b++) {
		     fprintf(stderr,"observed[%c] = %d\n",
				AlphaChar(b,M->A),M->observed[i][b]);
		  }
		}
	}
	PutFModel(stderr,M);
	return INT4_MIN;	// temporary fix for overflow.
	fmodel_error("LemonFModel( )... this should not happen.");
}

Int4	OrangeFModel(fm_type M)
/** Sample a random column in model. */
{
	double	rand,total;
	Int4	i,k;

	for(total=0.0,k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL) total += 1.0;
	}
	/** total = (double) M->ncols; /***/
	rand = ((double) Random()/(double) RANDOM_MAX)*total;
	for(k=1,i=M->start; i<=M->end; i++,k++){
	     if(M->observed[i] != NULL) if((rand-= 1.0) <= 0.0) return k;
	}
	fmodel_error(" RandColFModel( )... this should not happen.");
}

double	RatioFModel(Int4 *observed, Int4 d, fm_type M)
{ return exp(LnRatioFModel(observed, d, M)); }

double	LnRatioFModel(Int4 *observed, Int4 d, fm_type M)
// Return the ratio of new to old observed.
{
	double	sum;
	Int4	b,j,*mobserved;

	if(observed == NULL || d < 1 || d > M->length) return 0.0;
	if(M->recalc){
	    for(j=0; j<=M->maxlen; j++) M->sumlngamma[j] = FMODEL_UNDEF;
	    M->recalc = FALSE;
	}
	j = d+M->start-1;
	if(M->sumlngamma[j] == FMODEL_UNDEF){
	  if((mobserved = M->observed[j])==NULL) return 0.0;
	  for(M->sumlngamma[j]=0.0,b=1;b<= nAlpha(M->A); b++) {
	     if(M->Ps[b] != 0.0){
		M->sumlngamma[j] += lngamma(((double)mobserved[b]+M->Ps[b]));
	     }
	  }
	} 
	for(sum=0.0,b=1;b<= nAlpha(M->A); b++) {
	   if(M->Ps[b]!=0.0) sum+=lngamma(((double)observed[b]+M->Ps[b]));
	}
#if 0
	double D = (sum-M->sumlngamma[j]);
	if(D > 27.0){
		 fprintf(stderr,"ratio[%d] = %g - %g = %g; j=%d\n",
			j,sum,M->sumlngamma[j],D);
	}
	return D;
#endif
	return (sum-M->sumlngamma[j]);
}

Int4	MvColumnFModel(Int4 *observed, Int4 lemon, Int4 pos, fm_type M)
/*  Remove column at position lemon and add observed to pos in M */
{
	Int4	oms,ms,oldstart;
	
	oldstart = M->start;
	if(!RmColumnFModel(lemon, M)) fmodel_error("can't remove column!?");
	if(oldstart != M->start) pos += oldstart - M->start;
	oms = M->start; 	/* center_model may alter oms */
	ms = add_column_fmodel(observed, pos, M);
	if(ms!=oms) oldstart += ms - oms; /* re-adjust oldstart */
	return  M->start - oldstart;
}

Int4	RmColumnFModel2(Int4 lemon, fm_type M)
/** remove column and return change in start site **/
{
	Int4	oldstart;
	oldstart = M->start;
	if(!RmColumnFModel(lemon, M)) fmodel_error("can't remove column!?");
	return  M->start - oldstart;
}

BooLean RmColumnFModel(Int4 pos, fm_type M)
/*  if possible position "pos" is removed from model and TRUE is returned;
    if position "pos" cannot be removed FALSE is returned */
{
	Int4	site;

	if(M->length == 1) fmodel_error("zero length model not allowed.");
	site = pos + M->start - 1;
	if(site < M->start || site > M->end)  return FALSE;
	else if(M->observed[site]!=NULL){
		free(M->observed[site]); free(M->likelihood[site]); 
		M->observed[site] = NULL; M->likelihood[site] = NULL; 
		if(site == M->start){		/* shrink from left */
			while(M->observed[site]==NULL) site++;
			M->start=site; 
			M->length = M->end - M->start + 1;
		}else if(site == M->end){	/* shrink from right */
			while(M->observed[site]==NULL) site--;
			M->end=site; 
			M->length = M->end - M->start + 1;
		}
		M->ncols--;
		M->update = TRUE;
		return TRUE;
	} else return FALSE;
}

Int4	AddColumnFModel(Int4 *observed, Int4 pos, fm_type M)
/** add column and return change in start site **/
{
	Int4	oldstart;
	oldstart=add_column_fmodel(observed, pos, M);
	return  M->start - oldstart;
}

void	ShiftFModel(Int4 *observed, BooLean left, fm_type M)
{
	if(left){
		RmColumnFModel(1, M);
		add_column_fmodel(observed, M->length+1, M);
	} else {
		RmColumnFModel(M->length, M);
		add_column_fmodel(observed, 0, M);
	}
	M->update = TRUE;
}

/****************************** Output ********************************/

void	PutFModelShort(FILE *fptr, fm_type M)
/* Report the current frequency model. */
{
	Int4	j,b,r,pos,v;
	double	i,total,info,p,q;

	if(M->update) update_fmodel(M);
	/*** freq model ****/
	fprintf(fptr," pos ");
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += (double) M->observed[j][b] + M->Ps[b];
	      if(j==M->start) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(j==M->start) fprintf(fptr,"  Info\n");
	    fprintf(fptr,"%4d ",pos);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double) M->observed[j][b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			i = p*log(p/q)/log(2.0);
			info += i;
		}
		v = (Int4)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
	    	else fprintf(fptr,"%3d", (Int4)(100*p+0.5));
	    }
	    fprintf(fptr,"   %1.1f\n", info);
	  }
	}
	fprintf(fptr," ave.");	
	for(total=0.0,b = 1; b <= nAlpha(M->A); b++){
		total += (double) M->observedNS[b];
	}
	for(b = 1; b <= nAlpha(M->A); b++){
	    	p = (double) M->observedNS[b]/total;
		if(p > 0.0){ q = M->freq[b]; i = p*log(p/q)/log(2.0); }
		else i=0.0;
		v = (Int4)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else {
			r = (Int4) (100.0 * M->freq[b]);
	    		fprintf(fptr,"%3d", r);
		}
	} fprintf(fptr,"\n");
}

void	PutFModel(FILE *fptr, fm_type M)
/* Report the current frequency model. */
// use PutFModelShort eventually...but need to test it first...
{
	Int4	j,b,r,pos,v;
	double	i,total,info,p,q;
	h_type	H;

	if(M->update) update_fmodel(M); /*** freq model ****/
#if 0	// use this instead eventually...
	fprintf(fptr,"Motif model (residue frequency x 100):\n");
	PutFModelShort(fptr,M);
#else
	fprintf(fptr,"Motif model (residue frequency x 100):\n");
	fprintf(fptr,"POS  ");
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += (double) M->observed[j][b] + M->Ps[b];
	      if(j==M->start) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(j==M->start) fprintf(fptr,"  Info\n");
	    fprintf(fptr,"%4d ",pos);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double) M->observed[j][b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			i = p*log(p/q)/log(2.0);
			info += i;
			v = (Int4)floor(10*i+0.5);
		} else v = 0;
		if(v<=0) fprintf(fptr,"  .");
	    	else fprintf(fptr,"%3d", (Int4)(100*p+0.5));
	    }
	    fprintf(fptr,"   %1.1f\n", info);
	  }
	}
	fprintf(fptr,"non-\nsite ");	
	for(total=0.0,b = 1; b <= nAlpha(M->A); b++){
		total += (double) M->observedNS[b];
	}
	for(b = 1; b <= nAlpha(M->A); b++){
	    	p = (double) M->observedNS[b]/total;
		if(p > 0.0){ q = M->freq[b]; i = p*log(p/q)/log(2.0); }
		else i=0.0;
		v = (Int4)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else {
			r = (Int4) (100.0 * M->freq[b]);
	    		fprintf(fptr,"%3d", r);
		}
	}
	fprintf(fptr,"\n%d columns\n\n",M->ncols);
#endif
	H =Histogram("model information",0,LenFModel(M)+1,1.0);
	fprintf(fptr,
	  "Information (relative entropy) contribution in 1/100th bits:\n");
	fprintf(fptr,"POS  ");
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += (double) M->observed[j][b] + M->Ps[b];
	      if(j==M->start) fprintf(fptr,"%3c", AlphaChar(b, M->A));
	    }
	    if(j==M->start) fprintf(fptr,"  Info\n");
	    fprintf(fptr,"%4d ",pos);
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double) M->observed[j][b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			i = p*log(p/q)/log(2.0);
			info += i;
		} else i = 0.0;
		v = (Int4)floor(10*i+0.5);
		if(v<=0) fprintf(fptr,"  .");
		else fprintf(fptr,"%3d", v);
	    }
	    v = (Int4)floor(100*info+0.5);
	    fprintf(fptr,"  %3d\n", v);
	    IncdMHist(pos, v, H);
	  }
	}
	PutHist(fptr,60,H); NilHist(H);
}

/****************************** Statistics ********************************/
/** see macro LikelihoodFModel(seq, pos, M) & likelihood_fmodel below. **/

double	NormLikelihoodFModel(fm_type M)
{
	Int4	s,r;
        double	L=1.0,d,sum,n,**likelihood,min;

	if(M->update) update_fmodel(M);
	min=(0.2/(double)nAlpha(M->A));
	likelihood = M->likelihood + M->start;
	s = M->end - M->start;
        while(s >= 0){
                if(likelihood[s]!=NULL){
		   for(sum=n=0.0,r=nAlpha(M->A); r > 0; r--){
			 d=likelihood[s][r];
			 if(d >= min) { sum+=d; n+=1.0; }
		   } L*=sum/n;
		} s--; 
        } return L;
}

double	ProbFModel(register unsigned char *seq, register Int4 pos, 
	register double p, register fm_type M)
/* Return the relative probability of site at pos in seq being in model M */
{
        register double L = LikelihoodFModel(seq, pos, M);
	return (L*p/((1.0-p)+L*p));
}

Int4	*SeeColumnFModel(Int4 c, fm_type M)
{if(c<1||c>LenFModel(M))return 0;else return M->observed[M->start+c-1];}

Int4	ObservedFModel(Int4 **array, fm_type M)
{
	register Int4	i,j;

	if(M->update) update_fmodel(M);
	array[0] = 0;
        for(i=0,j=M->start; j<=M->end; j++){
		i++; array[i] = M->observed[j];
	} return i;
}

void	SetPseudoFModel(double npseudo, fm_type M)
{ 	char b;

	M->npseudo = npseudo;
	for(b=1; b<= nAlpha(M->A); b++) {
		M->Ps[b]= M->npseudo * M->freq[b];
	}
	M->update=TRUE;
}

#if 0	// modify LikelihoodFModel( ) so that 'X' residues count as deletions.
#define LikelihoodFModel(s,p,M) ( ((M)->update) ? update_fmodel(M), \
        likelihood_fmodel((s+p),((M)->end-(M)->start),((M)->likelihood+(M)->start)):\
        likelihood_fmodel((s+p),((M)->end-(M)->start),((M)->likelihood+(M)->start)))
#endif

/********************************* private ********************************/
/*****************************************************************
  Return the likelihood of site at pos in seq being in model M 
  WARNING: variables s and likelihood are assumed to have private 
	information:
	i.e., s = M->start - M->end.
	i.e., likelihood = M->likelihood + M->start.
/*****************************************************************/
double	likelihood_fmodel(register unsigned char *seq, register Int4 s, 
	register double **likelihood)
{
        register double L=1.0;

        while(s >= 0){
                if(likelihood[s]!=NULL) L*=likelihood[s][(seq[s])];
 		s--; 
        }
	return L;
}

double	score_fmodel(register unsigned char r, register Int4 s, 
	register double **likelihood) 
// return Nonsite if == NULL ????
{ if(likelihood[s]!=NULL) return log(likelihood[s][r]); else return 0.0; }

void    fmodel_error(const char *s) { fprintf(stderr,"fmodel: %s\n",s); exit(1); }

Int4	add_column_fmodel(Int4 *observed, Int4 pos, fm_type M)
/********************************************************************
 Add the column observed at pos in model M.  The center_model( ) 
 function may change the value of M->start therefore add_column_fmodel( )
 returns the value of the (possibly changed) old M->start to the 
 calling environment which may need this value.

 Note: 0 = one position to left of start!
*********************************************************************/
{
	Int4	site,oldstart=M->start;

	site = pos + M->start - 1;
	if(site < 1 || site > M->maxlen){
		fprintf(stderr,"\ncentering model\n");
		center_model(pos, M);
		site = pos + M->start - 1;
		oldstart=M->start;
	}
	if(site < M->start){ 
		if(site==0) fprintf(stderr,"site=0\n");
		M->start = site; 
		M->length = M->end - M->start + 1;
	} else if(site > M->end){	
		M->end = site; 
		M->length = M->end - M->start + 1;
	} else if(M->observed[site]!=NULL) 
		fmodel_error("attempt to add column where one exists!?");
	NEW(M->likelihood[site],nAlpha(M->A) +2, double);
        M->likelihood[site][0] = 0.05; /*** 'X' residues ***/
	M->observed[site] = observed; M->ncols++;
	M->update = TRUE;
	return oldstart;
}

void	center_model(Int4 pos, fm_type M)
/************************************************************************
 Recenters the model so that don't run off end of array. site is location
 of site to be added that is beyond one of the ends. 
 CAUTION: This operation clears all sumlngammas.
 ************************************************************************/
{
	Int4	length,flank,start,end,New,old,site;

	site = pos + M->start - 1;
	if(site < 1){			/** case I: adding column to left **/
	   length = M->end - site +1;
	   flank = M->maxlen - length;
	   start = MAXIMUM(Int4,1,flank/2); /*** try to put about in center ***/
	   if((start+pos)<2) start=-pos+2; /** pos < 0 -> start > 0 **/
#if 1
	   fprintf(stderr,
  	     "(L) start: %d->%d; leng: %d->%d (max=%d); site: %d->%d\n",
	     M->start,start,M->length,length,M->maxlen,site,pos+start-1);
#endif

	} else if(site > M->maxlen) {	/** case II: adding column to right **/
	   length = site - M->start +1;
	   flank = M->maxlen - length;
	   start = MAXIMUM(Int4,1,flank/2); /*** try to put about in center ***/

#if 1
	   fprintf(stderr,
  	     "(R) start: %d->%d; leng: %d->%d (max=%d); site: %d->%d\n",
	     M->start,start,M->length,length,M->maxlen,site,pos+start-1);
#endif

	} else return;
	M->recalc = TRUE;
	if(length > M->maxlen){ 
		fprintf(stderr,"site=%d;len=%d;maxlen=%d;ncols=%d\n",
			site,M->length,M->maxlen,M->ncols);
		// assert(length <= M->maxlen);
		fmodel_error("added site creates model > maxleng.");
	}
	end = start + M->length - 1;
	if(end > M->maxlen) fmodel_error("centering error; end > maxlen");
	if(start < M->start){		/*** :...New<-old: ***/
	    for(New = start,old= M->start; New <= end; New++,old++){
		M->observed[New] = M->observed[old];
		M->likelihood[New] = M->likelihood[old];
		M->observed[old] = NULL; M->likelihood[old] = NULL;
	    }
	} else if(start > M->start){  /*** :old->New...: ***/
	    for(New = end,old= M->end; New >= start; New--,old--){
		M->observed[New] = M->observed[old];
		M->likelihood[New] = M->likelihood[old];
		M->observed[old] = NULL; M->likelihood[old] = NULL;
	    }
	}
	M->start = start; M->end = end;
}

void	update_fmodel2(fm_type M)
/* Normalize frequencies to avoid overflow */
// THIS SITE ONLY versus ALL ELSE
{
	Int4		j;
	register Int4	b,*observed;
	register double	*Ps,*likelihood,totalS,*targetNS,totalNS;

	Ps = M->Ps; 
	for(totalNS=0.0,b=nAlpha(M->A);b>0;b--) totalNS+=(double)M->counts[b];
	totalNS += M->npseudo - (double) M->totsites;
	for(b=nAlpha(M->A); b > 0; b--){	// nonsite model
	   if(Ps[b] > 0.0) M->targetNS[b]=((double)M->counts[b]+Ps[b]); 
	}
	targetNS = M->targetNS;
	totalS = (double) (M->totsites + M->npseudo);
        for(j=M->start; j<=M->end; j++){
	    if(M->observed[j] != NULL){
		likelihood=M->likelihood[j]; observed=M->observed[j];
		for(b=nAlpha(M->A); b > 0; b--){
               	    likelihood[b]=((observed[b]+Ps[b])/totalS)/
			((targetNS[b]-(double)observed[b])/totalNS);
           	} // 'X' residue likelihoods are set once for all (above) 
	    }
        }
	M->update = FALSE; 
}

void	update_fmodel3(fm_type M)
// Sample using Jun's LPR ...
{
	Int4		j;
	register Int4	b,*observed;
	register double	*Ps,*likelihood,totalS,*targetNS;
	double		totalNS;

	Ps = M->Ps; 

	for(totalNS=0.0,b=nAlpha(M->A);b > 0; b--){totalNS+=M->observedNS[b];}
	totalNS+=M->npseudo;
	for(b=nAlpha(M->A); b > 0; b--){
	   if(Ps[b] > 0.0) M->targetNS[b]=
		( (lngamma(M->observedNS[b] - 1.0 +Ps[b]) - lngamma(M->observedNS[b]+Ps[b])) 
			- (lngamma(totalNS - 1.0 ) - lngamma(totalNS)) );
	   else M->targetNS[b]=0.0;
	} targetNS = M->targetNS;

	totalS = (double) (M->totsites + M->npseudo);
        for(j=M->start; j<=M->end; j++){
	    if(M->observed[j] != NULL){
		likelihood=M->likelihood[j]; observed=M->observed[j];
		for(b=nAlpha(M->A); b > 0; b--){
		   likelihood[b]=exp(
			( (lngamma(observed[b]+1.0+Ps[b]) - lngamma(observed[b]+Ps[b]))
				- (lngamma(totalS+1.0) - lngamma(totalS))  ) + targetNS[b]);
           	} // 'X' residue likelihoods are set once for all (above) 
	    }
        }
	M->update = FALSE; 
}

void	update_fmodel(fm_type M)
/* Normalize frequencies to avoid overflow */
{
	Int4		j;
	register Int4	b,*observed;
	register double	*Ps,*likelihood,totalS,*targetNS;
	double		totalNS;

	Ps = M->Ps; 
	for(totalNS=0.0,b=nAlpha(M->A);b > 0; b--){totalNS+=M->observedNS[b];}
	totalNS+=M->npseudo;
	for(b=nAlpha(M->A); b > 0; b--){
	   if(Ps[b] > 0.0) M->targetNS[b]=(M->observedNS[b]+Ps[b])/totalNS;  
	   else M->targetNS[b]=1.0;  // this should not be used...
	} targetNS = M->targetNS;
	totalS = (double) (M->totsites + M->npseudo);
        for(j=M->start; j<=M->end; j++){
	    if(M->observed[j] != NULL){
		likelihood=M->likelihood[j]; observed=M->observed[j];
		for(b=nAlpha(M->A); b > 0; b--){
               	    likelihood[b]=((observed[b]+Ps[b])/totalS)/targetNS[b];
		    //assert(likelihood[b] < DBL_MAX && likelihood[b] > -DBL_MAX);
           	} // 'X' residue likelihoods are set once for all (above) 
	    }
        }
	M->update = FALSE; 
}

smx_typ	SampleSmatrixFModel(double pernats, Int4 wt, fm_type M)
// return an smatrix with parameters sampled from the beta distribution ...
{
        Int4	i,j,c,o;
	Int4	*observed,total;
	double	s,*likelihood,*Ps,*targetNS,totalS,obs,totalNS;
	static Int4 Seed=0;

	if(Seed==0) { 
	   Seed=-Random();  // need to initialize with a negative number.
	   fprintf(stderr,"INITIALIZING SEED (%d)\n",Seed);
	}
	assert(wt > 0 && wt <= 20);
	Ps = M->Ps; 
	smx_typ	smx = MkSMatrix(2.0,M->length,M->freq,M->A);
	for(totalNS=0.0,c=nAlpha(M->A);c > 0; c--){totalNS+=M->observedNS[c];}
	totalNS+=M->npseudo;
	for(c=nAlpha(M->A); c > 0; c--){
	   if(Ps[c] > 0.0) M->targetNS[c]=(M->observedNS[c]+Ps[c])/totalNS;  
	   else M->targetNS[c]=1.0;  // this should not be used...
	} targetNS = M->targetNS;
	total=M->totsites;
	totalS = (double) (M->totsites + M->npseudo);
        for(i=1,j=M->start; j<=M->end; j++,i++){
            if(M->observed[j]!=NULL){
	       observed=M->observed[j];
               for(c=nAlpha(M->A); c>0; c--){
		if(Ps[c] > 0.0){
		  o=observed[c];
		  if(o > 0 && o < total){
		    obs=(double)total*betadev(o*wt,(total-o)*wt,&Seed);
		    if(obs > 0.0){
		      s=((obs+Ps[c])/(totalS+(obs-o)))/targetNS[c];
		    } else s=((o+Ps[c])/totalS)/targetNS[c];
		  } else if(o) s=1.0/targetNS[c]; else s=(Ps[c]/totalS)/targetNS[c];
		  SetSMatrix(c,i,(Int4)floor((pernats*log(s)+0.5)),smx);
                } else SetSMatrix(c,i,0,smx);
	       }
	       if(Ps[0] == 0.0) SetSMatrix(0,i,0,smx);
	       else SetSMatrix(0,i,(Int4)floor((pernats*log(M->likelihood[j][0]))+0.5),smx);
             } else for(c=0; c<=nAlpha(M->A); c++) SetSMatrix(c,i,0,smx);
        }
	// PutFModel(stderr, M); PutSMatrix(stderr, smx); if(TRUE) exit(1);
	return smx;
}

smx_typ	SampleDirichletSmatrixFModel(double pernats, Int4 wt, fm_type M)
/***********************************************************************
 Returns an smatrix with parameters sampled from the Dirichlet distribution
 corresponding to the observed counts with uninformed priors (i.e., with
 one pseudocount per residue type).
 ***********************************************************************/
{
        Int4	i,j,c,*observed;
	static Int4 Seed=0;
	UInt4	*observedNS,total;
	double	*p,*q;
	Int4	*a;

	if(Seed==0) { 
	   Seed=-Random();  // need to initialize with a negative number.
	   fprintf(stderr,"INITIALIZING SEED (%d)\n",Seed);
	}
	assert(wt <= 200);
	smx_typ	smx = MkSMatrix(2.0,M->length,M->freq,M->A);

	NEW(p,nAlpha(M->A)+3,double);
	NEW(q,nAlpha(M->A)+3,double);
	NEW(a,nAlpha(M->A)+3,Int4);
	observedNS=M->observedNS;
        for(total=0,c=nAlpha(M->A); c>0; c--) total+=a[c]=(observedNS[c]*wt + 1); 
	if(wt <= 0){	// Then take maximum likelihood...
	   for(c=nAlpha(M->A); c>0; c--) q[c] = (double)a[c]/(double)total;
	} else {
	  DirichletDev(a,q,nAlpha(M->A),&Seed);
	  for(c=nAlpha(M->A); c>0; c--) if(q[c] < 0.0001) q[c] = 0.0001;
	}
        for(i=1,j=M->start; j<=M->end; j++,i++){
            if(M->observed[j]!=NULL){
	       observed=M->observed[j];
	       if(wt <= 0){	// Then take maximum likelihood...
                  for(total=0,c=nAlpha(M->A); c>0; c--) total+=a[c]=(observed[c] + 1); 
		  for(c=nAlpha(M->A); c>0; c--) p[c]=(double)a[c]/(double)total;
	       } else {		// Sample parameters from Dirichlet distribution.
                  for(c=nAlpha(M->A); c>0; c--) a[c]=(observed[c]*wt + 1); 
	          DirichletDev(a,p,nAlpha(M->A),&Seed);
	       }
               for(c=nAlpha(M->A); c>0; c--){
		 if(p[c] < 0.0001) p[c] = 0.0001;
		 SetSMatrix(c,i,(Int4)floor((pernats*log(p[c]/q[c])+0.5)),smx);
#if 0
fprintf(stderr,"%d %c: p = %.5f q = %.5f; %.5f %.5f\n", i,AlphaChar(c,M->A),p[c],q[c],
		(double)(observed[c]+1)/(M->totsites+20),
		(double)(observedNS[c]+1)/total);
#endif
	       }
	       if(M->Ps[0] == 0.0) SetSMatrix(0,i,0,smx);
	       else SetSMatrix(0,i,(Int4)floor((pernats*log(M->likelihood[j][0]))+0.5),smx);
            } else for(c=0; c<=nAlpha(M->A); c++) SetSMatrix(c,i,0,smx);
        }
	free(p); free(a); free(q);
	// PutFModel(stderr, M); PutSMatrix(stderr, smx); if(TRUE) exit(1);
	return smx;
}

smx_typ	GetSmatrixFModel(double pernats, fm_type M)
{
        Int4	i,j,score,c;
	double	s;
	smx_typ	smx;

	if(M->update) update_fmodel(M);
	smx = MkSMatrix(2.0,M->length,M->freq,M->A);
        for(i=1,j=M->start; j<=M->end; j++,i++){
            if(M->observed[j]!=NULL){
               for(c=0; c<=nAlpha(M->A); c++){
		if(M->Ps[c] > 0.0){
		    s = M->likelihood[j][c];
                    score = (Int4) floor((pernats*log(s) + 0.5));
		    SetSMatrix(c,i,score,smx);
#if 0
if(score < -999){ 
std::cerr << "M->likelihood[j][c] = "; std::cerr << M->likelihood[j][c]; 
std::cerr << "\ns = "; std::cerr << s; std::cerr << std::endl;
PutFModel(stderr, M); PutSMatrix(stderr,smx); 
}
#endif
                } else SetSMatrix(c,i,0,smx);
	       }
             } else for(c=0; c<=nAlpha(M->A); c++) SetSMatrix(c,i,0,smx);
        }
	// PutSMatrix(stderr, smx);
	return smx;
}

float	*InfoFModel(fm_type M)
/* Return the information content in tenth bits. */
{
	Int4	j,b,r,pos,v;
	double	total,info,p,q;
	float	*Info;

	if(M->update) update_fmodel(M);
	NEW(Info, LenFModel(M) +3,float);
	for(pos=1,j=M->start; j<= M->end; j++,pos++){
	  if(M->observed[j]!=NULL){
	    for(total=0.0, b = 1; b <= nAlpha(M->A); b++){
	      total += (double) M->observed[j][b] + M->Ps[b];
	    }
	    for(info=0.0,b = 1; b <= nAlpha(M->A); b++){
	    	p = ((double)M->observed[j][b]+M->Ps[b])/total;
		if(p > 0.0){
			q = M->freq[b];
			info += p*log(p/q);
		}
	    }
	    Info[pos] = info*1.442695; // in bits
	  } else Info[pos] = -1.0;
	} return Info;
}

