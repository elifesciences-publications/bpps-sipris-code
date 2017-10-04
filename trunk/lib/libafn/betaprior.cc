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

#include "betaprior.h"

bp_type MkBPrior(Int4 expect, double weight, double N)
/*******************************************************************
  N = total sites; 	- Pseudo and Real sites.
  tot_sites  = number of sites;
  alpha = number site pseudo counts;
  beta = number site pseudo counts;
 *******************************************************************/
{
	bp_type	B;

	if(N <=0) bprior_error("total sites must be > 0");
	NEW(B,1,beta_prior_type);
        B->expect = expect;
        B->weight = weight; 
	SetBPriorN(N, B);
	B->success = 0;
        B->p = (B->alpha + (double)B->expect)/B->T;
	B->calc = FALSE;
	return B;
}

bp_type CopyBPrior(bp_type B)
/*******************************************************************
 Make and return an exact copy of B.
 *******************************************************************/
{
	bp_type	B2;

	NEW(B2,1,beta_prior_type);
        B2->A = B->A; B2->N = B->N; B2->T = B->T;
        B2->alpha = B->alpha; B2->beta = B->beta;
        B2->success = B->success; B2->expect = B->expect;
        B2->p = B->p;
        B2->weight = B->weight;
        B2->calc = B->calc;
	return B2;
}

double	BPriorNullMAP(register bp_type B)
{ return LnBeta(B->alpha,B->N + B->beta) - LnBeta(B->alpha,B->beta); }

double	RatioBPrior(register bp_type B1, register bp_type B2)
/** Return the ratio of new to old. **/
{
	register double      v1,v2;

        v1 = (double) B1->success;
        v2 = (double) B2->success;
        v1 = LnBeta(v1+B1->alpha,B1->N - v1 + B1->beta)
		- LnBeta(B1->alpha,B1->beta);
        v2 = LnBeta(v2+B2->alpha,B2->N - v2 + B2->beta)
		- LnBeta(B2->alpha,B2->beta);
	return exp(v1 - v2);
}

double	BPriorMAP(register bp_type B)
{
	register double      n;

        n = (double) B->success;
        return LnBeta(n+B->alpha,B->N - n + B->beta) - LnBeta(B->alpha,B->beta);
}

double	LogLikeBPrior(bp_type B)
/** in bits of information **/
{
	double	n,p,L,N;

        n = (double) B->success;
	N = (double)B->N;
	p = n/N;
	L = n*log(p) + (N-n)*log(1.0 - p);
        return 1.4427*L;
}

void    SetBPriorN(double N, bp_type B)
{ 
        double  ratio;

        ratio = (B->weight/(1.0 - B->weight));
	B->N = N;
        B->A =N*ratio;		/* A = N*(w/(1-w)) */
        B->alpha = (double) B->expect*ratio;
        B->beta= B->A - B->alpha;
        B->T=(B->A + B->N);
	B->calc = TRUE; 
}

void    SetBPriorS(Int4 success, bp_type B)
{ B->success = success; B->calc = TRUE; }

void    ClearBPrior(bp_type B)
{
	B->success = 0;
        B->p = (B->alpha + (double)B->expect)/B->T;
        B->calc = FALSE;
}

void	PutBPrior(FILE *fptr, bp_type B)
{
        double  variance,stdev;

        variance = (B->alpha + (double)B->expect);
        variance *= (B->beta + B->N - (double)B->expect);
        variance /= (B->A+B->N)*(B->A+B->N)*(B->A+B->N+1.0);
        stdev = sqrt(variance);
        stdev *= 2*B->N;
        fprintf(fptr, "%d (+/-%4.1f) out of %.0f ",
                        B->expect, stdev,B->N);
        fprintf(fptr,"  a = %g; b = %g; p = %g\n",
                B->alpha, B->beta, B->p);
}

void    bprior_error(char *s) {fprintf(stderr,"betaprior: %s\n",s); exit(1);}

