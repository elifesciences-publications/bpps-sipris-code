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

#include "profile.h"

pfl_typ	MkProfile(Int4 N, Int4 K, double *freq, a_type A)
{
	pfl_typ	M;
	Int4	i;
	char	c;
	
	NEW(M,1,profile_type);
	M->A = A;
	M->nlet = nAlpha(A); M->K = K;
	M->smx = MkSMatrixN(N,K,freq,A);
	M->calc_prof = TRUE; M->nseqs = 0;
	NEWP(M->cnts,K+2,Int4);
	NEW(M->maxseg,K+2,char);
	for(i=1; i<=K;i++) {
		NEW(M->cnts[i],M->nlet+2,Int4);
		for(c=0; c<=M->nlet; c++) M->cnts[i][c]=0;
	}
	NEW(M->max,K+2,Int4); NEW(M->min,K+2,Int4); 
	return M;
}

void	NilProfile(pfl_typ M)
{
	Int4	i;

	for(i=1; i <= M->K;i++) { free(M->cnts[i]);}
	free(M->cnts); free(M->min); free(M->max); free(M->maxseg);
	NilSMatrix(M->smx);
	free(M);
}

void	RmProfile(unsigned char *seq, Int4 start, pfl_typ M)
{
	Int4	j,s;

	for(s=start,j=1; j <= M->K; s++,j++) M->cnts[j][seq[s]]--;
	M->nseqs--; M->calc_prof = TRUE;
}

void	Add2Profile(unsigned char *seq, Int4 start, pfl_typ M)
{
	Int4	j,s;

	for(s=start,j=1; j <= M->K; s++,j++) M->cnts[j][seq[s]]++;
	M->nseqs++; M->calc_prof = TRUE;
}

void	calc_profile(pfl_typ M)
/** Chip()? method *****************************
	We want to find log[q(r)/p(r)] - the odds of r being in the
	model (i.e., q(r)) versus being in a nonsite (i.e., p(r). 
		q(r) = sum q(x)*p(r|x)		parental residue is x
			x

		     = sum q(x)*p(r,x)/p(x)	by def.
			x

	-> q(r)/p(r) = sum q(x)*p(r,x)/p(x)*p(r)
			x

		     = sum q(x)*blosum62L(x,r)	by def.
		        x
*******************************************************************/
{
	Int4	i,s,r,x;
	double	sum,q;
	
	M->max[0]=M->min[0]=0;
	for(i = 1; i<= M->K; i++){
	    M->max[i]=INT_MIN; M->min[i]=INT_MAX;
	    for(r = 1; r<= M->nlet; r++){
		/**** GRIBSKOV'S PROFILE METHOD ****/
		for(sum = 0.0,x = 1; x<= M->nlet; x++){
		    sum  += (double)M->cnts[i][x]*blosum62[r][x];
		}
		s = (Int4) floor(2.0*sum + 0.5);
		/**** CHIP'S()? METHOD ****
		for(sum = 0.0,x = 1; x<= M->nlet; x++){
		    q  = ((double)M->cnts[i][x]/(double)M->nseqs);
		    sum += q * (double) blosum62L[r][x];
		}
		s = (Int4) floor((10.0*log(sum) + 0.5));
		/***************************/
		M->min[i] = MINIMUM(Int4,M->min[i],s);
		if(M->max[i] < s){ M->max[i]=s; M->maxseg[i]=r; }
		SetSMatrix(r,i,s,M->smx);
	    }
	    if(M->min[i] < 0) {
		SetSMatrix(0,i,2*M->min[i],M->smx);
	    } else {
		SetSMatrix(0,i, M->min[i], M->smx);
	    }
	}
	M->calc_prof = FALSE;
}

void	PutProfile(FILE *fptr, pfl_typ M)
{ if(M->calc_prof) calc_profile(M); PutSMatrix(fptr,M->smx); }

Int4	ScoreProfile(register unsigned char *seq, register Int4 start, register pfl_typ M)
{ if(M->calc_prof) calc_profile(M); return ScoreSMatrix(seq,start,M->smx); }

double	ProfileProb(Int4 score, pfl_typ M)
{ if(M->calc_prof) calc_profile(M); return SMatrixProbFast(score,M->smx); }

/****************************** private **********************************/

void	profile_error(char *s){fprintf(stderr,"Profile: %s\n",s);exit(1);}

